#!/usr/bin/env python
from __future__ import print_function
import unittest
import yaml
import json
import subprocess
import os
import dxpy


TEST_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
TEST_DATA_YAML = os.path.join(TEST_DIRECTORY, '../../../test_data.yml')
MAIN_APP_DIR = os.path.join(TEST_DIRECTORY, '../')


def remove_output(output_value):
    if type(output_value) is list:
        map(remove_output, output_value)
    elif dxpy.is_dxlink(output_value):
        dxpy.remove(output_value)
    elif type(output_value) is dict:
        map(remove_output, output_value.values())


class VariantAnnotatorTests(unittest.TestCase):
    def setUp(self):
        print('Building applets')
        with open(TEST_DATA_YAML) as fh:
            # Load the test data.
            self.test_data = yaml.safe_load(fh)['production']

            # Check if we are on staging or production
            if dxpy.APISERVER_HOST.startswith('staging') and 'staging' in self.test_data:
                self.test_data = self.test_data['staging']
            elif 'production' in self.test_data:
                self.test_data = self.test_data['production']

            # Get the region of interest
            region = dxpy.describe(dxpy.PROJECT_CONTEXT_ID)['region']
            if region in self.test_data:
                self.test_data = self.test_data[region]

        self.applet = json.loads(subprocess.check_output(['dx', 'build', '-f', MAIN_APP_DIR]).strip())['id']


    def tearDown(self):
        print('Removing applet')
        subprocess.check_call(['dx', 'rm', self.applet])
        print('Removing job outputs.')
        map(lambda job: remove_output(job.describe()['output']), self.jobs)


    def test_applet(self):
        self.assertIsNotNone(self.test_data)

        applet = dxpy.DXApplet(self.applet)
        self.jobs = []
        job_input = {'variants_vcfgzs': [dxpy.dxlink(self.test_data['hg001_wes_vcfgz']), dxpy.dxlink(self.test_data['hg002_wes_vcfgz'])],
                     'variants_vcfgztbis': [dxpy.dxlink(self.test_data['hg001_wes_vcfgz_tbi']), dxpy.dxlink(self.test_data['hg002_wes_vcfgz_tbi'])],
                     'genome_fastagz': dxpy.dxlink(self.test_data['hs37d5_fasta']),
                     'gatk_jar': dxpy.dxlink(self.test_data['gatk_jar']),
                     'no_geno': False,
                     'hdr_only': False}
        self.jobs.append(applet.run(job_input))
        print('Running full test in job: {0}.'.format(self.jobs[-1].get_id()))

        job_input['hdr_only'] = True
        self.jobs.append(applet.run(job_input))
        print('Running header only in job: {0}.'.format(self.jobs[-1].get_id()))

        job_input['no_geno'] = True
        self.jobs.append(applet.run(job_input))
        print('Running header only, no genotyping test in job: {0}.'.format(self.jobs[-1].get_id()))

        # Note that if we specify hdr_only = False and no_geno= True (the only combination we
        # aren't testing here) then no output is produced and the job fails.

        for job in self.jobs:
            job.wait_on_done()
