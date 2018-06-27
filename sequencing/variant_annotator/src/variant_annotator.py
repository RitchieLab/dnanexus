import os
import subprocess
import re
import multiprocessing
from itertools import izip
from collections import defaultdict

import dx_utils
import dxpy


@dxpy.entry_point('run_anno')
def run_anno(**job_inputs):
    # Download input files
    gatk_jar = dx_utils.download_and_gunzip_file(job_inputs['gatk_jar'])
    vcf_fn = dx_utils.download_and_gunzip_file(job_inputs['vcf'], skip_decompress=True)
    dx_utils.download_and_gunzip_file(job_inputs['tbi'])
    fasta_fn = dx_utils.download_and_gunzip_file(job_inputs['genome_fastagz'])
    dx_utils.download_and_gunzip_file(job_inputs['fasta_fai'])
    dx_utils.download_and_gunzip_file(job_inputs['fasta_dict'])

    # Calculate some local variables
    ofn = '{0}.annotated.vcf.gz'.format(re.sub('.vcf.gz', '', vcf_fn))
    mem = int(dx_utils.get_memory('M') * 0.9)

    # Form and run our GATK VariatAnnotator command
    cmd = ['java', '-Xmx{0}m'.format(mem), '-jar', gatk_jar, '-T', 'VariantAnnotator', 
        '-nt', str(multiprocessing.cpu_count()), '-R', fasta_fn, '-V', vcf_fn, '-o', ofn]
    if job_inputs['hdr_only'] and job_inputs['no_geno']:
        cmd += ['--sites_only']
    if 'cmd_params' in job_inputs:
        cmd += job_inputs['cmd_params'].split()
    dx_utils.run_cmd(cmd)

    # Upload the output
    output = {}
    if job_inputs['hdr_only']:
        if job_inputs['no_geno']:
            dx_utils.run_cmd(['mv', ofn, 'header.{0}'.format(ofn)])
            dx_utils.run_cmd(['mv', ofn + '.tbi', 'header.{0}.tbi'.format(ofn)])
        else:
            cmd = ['pigz', '-dc', ofn]
            pigz_proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            cmd = ['cut', '-f1-8']
            cut_proc = subprocess.Popen(cmd, stdin=pigz_proc.stdout, stdout=subprocess.PIPE)
            cmd = ['bgzip', '-c']
            with open('header.{0}'.format(ofn), 'w') as fh:
                bgzip_proc = subprocess.Popen(cmd, stdin=cut_proc.stdout, stdout=fh)
                bgzip_proc.communicate()
            cmd = ['tabix', '-p', 'vcf', 'header.{0}'.format(ofn)]
            dx_utils.run_cmd(cmd)
        output['vcf_hdr_out'] = dxpy.dxlink(dxpy.upload_local_file('header.{0}'.format(ofn)))
        output['vcfidx_hdr_out'] = dxpy.dxlink(dxpy.upload_local_file('header.{0}.tbi'.format(ofn)))

    if not job_inputs['no_geno']:
        output['vcf_out'] = dxpy.dxlink(dxpy.upload_local_file(ofn))
        output['vcfidx_out'] = dxpy.dxlink(dxpy.upload_local_file('{0}.tbi'.format(ofn)))

    return output


@dxpy.entry_point('main')
def main(**job_inputs):
    # Do some sanity checks
    if job_inputs['no_geno'] and not job_inputs['hdr_only']:
        raise dxpy.AppInternalError('ERROR: No output requested, nothing to do!')

    if len(job_inputs['variants_vcfgztbis']) != len(job_inputs['variants_vcfgzs']):
        raise dxpy.AppinternalError('ERROR: Number of VCFs and VCF indexes do NOT match!')

    # Download reference fasta
    genome_fasta_name = dx_utils.download_and_gunzip_file(job_inputs['genome_fastagz'])

    # Create reference genome indices
    pool = dx_utils.ExceptionAwarePool()
    pool.apply_async(dx_utils.run_cmd, (['samtools', 'faidx', genome_fasta_name],))
    pool.apply_async(dx_utils.run_cmd, (['java', '-jar', '/picard.jar', 'CreateSequenceDictionary',
        'R={0}'.format(genome_fasta_name), 'O={0}.dict'.format(os.path.splitext(genome_fasta_name)[0])],))
    pool.close()
    pool.join()
    fasta_fai = dx_utils.gzip_and_upload(genome_fasta_name + '.fai')
    fasta_dict = dx_utils.gzip_and_upload('{0}.dict'.format(os.path.splitext(genome_fasta_name)[0]))

    # Sort the vcf's and tbis by filename
    variants_vcfgzs = sorted(job_inputs['variants_vcfgzs'], key=lambda x: dxpy.describe(x)['name'])
    variants_vcfgztbis = sorted(job_inputs['variants_vcfgztbis'], key=lambda x: dxpy.describe(x)['name'])

    # Now launch a subjob per vcf/tbi pair.
    output = defaultdict(list)
    subjob_inputs = job_inputs.copy()
    subjob_inputs['fasta_fai'] = fasta_fai
    subjob_inputs['fasta_dict'] = fasta_dict
    for vcfgz, tbi in izip(variants_vcfgzs, variants_vcfgztbis):
        subjob_inputs['vcf'] = vcfgz
        subjob_inputs['tbi'] = tbi
        job = dxpy.new_dxjob(subjob_inputs, 'run_anno')
        if job_inputs['hdr_only']:
            output['vcf_hdr_out'].append(job.get_output_ref('vcf_hdr_out'))
            output['vcfidx_hdr_out'].append(job.get_output_ref('vcfidx_hdr_out'))
        if not job_inputs['no_geno']:
            output['vcf_out'].append(job.get_output_ref('vcf_out'))
            output['vcfidx_out'].append(job.get_output_ref('vcfidx_out'))

    return dict(output)