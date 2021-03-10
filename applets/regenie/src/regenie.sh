#!/bin/bash
# regenie 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://documentation.dnanexus.com/developer for tutorials on how
# to modify this file.

main() {

    ### validate and download input files
    # original filenames in `dx describe "variable" --name`

    mkdir input
    PLINK_ARGS="--threads $(nproc)"
    REGENIE_ARGS="--threads $(nproc)"
    REGENIE_ARGS_STEP1="--step 1"
    REGENIE_ARGS_STEP2="--step 2"

    if [[ -n "${input_bed}" ]]; then
        if [[ -z "${input_bim}" ]]; then
            dx-jobutil-report-error "Input variants data not found; the .bim input file must be specified when using .bed genotypes."
        elif [[ -z "${input_fam}" ]]; then
            dx-jobutil-report-error "Input samples data not found; the .fam input file must be specified when using .bed genotypes."
        fi
        dx download "${input_bed}" -o input/geno.bed
        dx download "${input_bim}" -o input/geno.bim
        dx download "${input_fam}" -o input/geno.fam
        PLINK_ARGS="${PLINK_ARGS} --bfile input/geno"
        REGENIE_ARGS="${REGENIE_ARGS} --bed input/geno"
    elif [[ -n "${input_bgen}" ]]; then
        if [[ -z "${input_sample}" ]]; then
            dx-jobutil-report-error "Input samples data not found; the .sample input file must be specified when using .bgen genotypes."
        fi
        dx download "${input_bgen}" -o input/geno.bgen
        dx download "${input_sample}" -o input/geno.sample
        PLINK_ARGS="${PLINK_ARGS} --bgen input/geno.bgen --sample input/geno.sample"
        REGENIE_ARGS="${REGENIE_ARGS} --bgen input/geno.bgen --sample input/geno.sample"
    else
        dx-jobutil-report-error "Input genotype data not found; either the .bed/.bim/.fam or .bgen/.sample input files must be specified."
    fi

    if [[ -n "${input_keep_samples}" ]]; then
        dx download "${input_keep_samples}" -o input/keep.samples
        PLINK_ARGS="${PLINK_ARGS} --keep input/keep.samples"
        REGENIE_ARGS="${REGENIE_ARGS} --keep input/keep.samples"
    fi
    if [[ -n "${input_remove_samples}" ]]; then
        dx download "${input_remove_samples}" -o input/remove.samples
        PLINK_ARGS="${PLINK_ARGS} --remove input/remove.samples"
        REGENIE_ARGS="${REGENIE_ARGS} --remove input/remove.samples"
    fi

    if [[ -n "${input_keep_variants}" ]]; then
        dx download "${input_keep_variants}" -o input/keep.variants
        PLINK_ARGS="${PLINK_ARGS} --extract input/keep.variants"
        REGENIE_ARGS_STEP1="${REGENIE_ARGS_STEP1} --extract input/keep.variants"
    fi
    if [[ -n "${input_remove_variants}" ]]; then
        dx download "${input_remove_variants}" -o input/remove.variants
        PLINK_ARGS="${PLINK_ARGS} --exclude input/remove.variants"
        REGENIE_ARGS_STEP1="${REGENIE_ARGS_STEP1} --exclude input/remove.variants"
    fi

    if [[ -n "${input_pred}" ]]; then
        if [[ "${run_step_1}" == "true" ]]; then
            echo "WARNING: ignoring input predictions file because step 1 will be run"
        else
            dx download "${input_pred}" -o input/pred.list
            REGENIE_ARGS_STEP2="${REGENIE_ARGS_STEP2} --pred input/remove.variants"
        fi
    elif [[ "${run_step_1}" != "true" && "${run_step_2}" == "true" ]]; then
        dx-jobutil-report-error "Input predictions file not found; step 1 predictions must be supplied when running step 2 only"
    fi

    dx download "${input_pheno}" -o input/pheno.txt
    REGENIE_ARGS="${REGENIE_ARGS} --phenoFile input/pheno.txt"
    if [[ -n "${pheno_columns}" ]]; then
        REGENIE_ARGS="${REGENIE_ARGS} --phenoColList ${pheno_columns}"
    fi

    dx download "${input_covar}" -o input/covar.txt
    REGENIE_ARGS="${REGENIE_ARGS} --covarFile input/covar.txt"
    if [[ -n "${covar_columns}" ]]; then
        REGENIE_ARGS="${REGENIE_ARGS} --covarColList ${covar_columns}"
    fi
    if [[ "${drop_missing_covars}" == "true" ]]; then
        grep -P "(^|\s)NA(\s|$)" input/covar.txt | awk '{print $1,$2}' > input/covar.drop.samples
        if [[ -s input/covar.drop.samples ]]; then
            echo "WARNING: $(cat input/covar.drop.samples | wc -l) samples have missing covariate(s) and will be dropped"
            PLINK_ARGS="${PLINK_ARGS} --remove input/covar.drop.samples"
            REGENIE_ARGS="${REGENIE_ARGS} --remove input/covar.drop.samples"
        fi
    fi

    REGENIE_ARGS="${REGENIE_ARGS} --bsize ${block_size}"
    if [[ "${flag_loocv}" == "true" ]]; then
        REGENIE_ARGS="${REGENIE_ARGS} --loocv"
    fi
    if [[ "${flag_bt}" == "true" ]]; then
        REGENIE_ARGS="${REGENIE_ARGS} --bt"
    fi
    if [[ "${flag_lowmem}" == "true" ]]; then
        mkdir lowmem
        REGENIE_ARGS="${REGENIE_ARGS} --lowmem --lowmem-prefix lowmem/tmp"
    fi
    if [[ -n "${max_iter}" ]]; then
        REGENIE_ARGS_STEP1="${REGENIE_ARGS_STEP1} --niter ${max_iter}"
    fi
    if [[ -n "${firth_limit}" ]]; then
        REGENIE_ARGS_STEP2="${REGENIE_ARGS_STEP2} --firth --pThresh ${firth_limit}"
        if [[ "${flag_firth_approx}" == "true" ]]; then
            REGENIE_ARGS_STEP2="${REGENIE_ARGS_STEP2} --approx"
        fi
    fi
    if [[ -n "${min_mac}" ]]; then
        REGENIE_ARGS_STEP2="${REGENIE_ARGS_STEP2} --minMAC ${min_mac}"
    fi

    echo "===== INPUTS ====="
    ls -la input


    ### run regenie step 1?
    mkdir output1
    if [[ "${run_step_1}" == "true" ]]; then
        if [[ "${auto_qc}" == "true" ]]; then
            ### apply MAF/LD QC for step 1
            echo "===== RUNNING PLINK QC: ${PLINK_ARGS} ====="
            plink ${PLINK_ARGS} --out input/qc
            REGENIE_ARGS_STEP1="${REGENIE_ARGS_STEP1} --extract input/qc.prune.in"
        fi

        REGENIE_ARGS_STEP1="${REGENIE_ARGS_STEP1} ${extra_options_1}"
        echo "===== RUNNING REGENIE STEP 1: ${REGENIE_ARGS} ${REGENIE_ARGS_STEP1} ====="
        for P in $(seq 1 $(head -n 1 input/pheno.txt | wc -w)) ; do
            touch output1/step1_${P}.loco
        done
        regenie ${REGENIE_ARGS} ${REGENIE_ARGS_STEP1} --out output1/step1
        
        echo "===== STEP 1 OUTPUTS ====="
        ls -la output1
        
        REGENIE_ARGS_STEP2="${REGENIE_ARGS_STEP2} --pred output1/step1_pred.list"
    fi


    ### run regenie step 2?
    mkdir output2
    if [[ "${run_step_2}" == "true" ]]; then
        REGENIE_ARGS_STEP2="${REGENIE_ARGS_STEP2} ${extra_options_2}"
        echo "===== RUNNING REGENIE STEP 2: ${REGENIE_ARGS} ${REGENIE_ARGS_STEP2} ====="
        regenie ${REGENIE_ARGS} ${REGENIE_ARGS_STEP2}  --out output2/step2

        echo "===== STEP 2 OUTPUTS ====="
        ls -la output2
    fi


    ### upload outputs

    if [[ -s input/covar.drop.samples ]]; then
        dropped_samples="$(dx upload --brief input/covar.drop.samples)"
        dx-jobutil-add-output dropped_samples "${dropped_samples}" --class=file
    fi
    if [[ -s output1/step1.log ]]; then
        output_log_1="$(dx upload --brief output1/step1.log)"
        dx-jobutil-add-output output_log_1 "${output_log_1}" --class=file
    fi
    if [[ -s output1/step1_pred.list ]]; then
        output_pred="$(dx upload --brief output1/step1_pred.list)"
        dx-jobutil-add-output output_pred "${output_pred}" --class=file
    fi
    if [[ -s output2/step2.log ]]; then
        output_log_2="$(dx upload --brief output2/step2.log)"
        dx-jobutil-add-output output_log_2 "${output_log_2}" --class=file
    fi
    if [[ -s output2/step2.regenie ]]; then
        output_regenie="$(dx upload --brief output2/step2.regenie)"
        dx-jobutil-add-output output_regenie "${output_regenie}" --class=file
    fi
}
