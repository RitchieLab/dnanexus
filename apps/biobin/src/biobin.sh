#!/bin/bash
set -e -x -o pipefail

main() {

  ulimit -c unlimited
  # switch to a temp directory and download all user input files

  NUM_CORES=$(nproc)
  TMPDIR=$(mktemp -d)
  cd "$TMPDIR"
  mkdir input
  BIOBIN_ROLE_ARG=""
  if [[ -n "$role" ]]; then
    dx download "$role_rol" -o input/input.role
    BIOBIN_ROLE_ARG="--role-file input/input.role"
  fi
  BIOBIN_PHENO_ARG=""
  BIOBIN_TEST_ARG=""
  if [[ -n "$phenotype_phe" ]]; then
    dx download "$phenotype_phe" -o input/input.phenotype
    BIOBIN_PHENO_ARG="--phenotype-file input/input.phenotype"
    if [[ ${#regression_types[*]} -gt 0 ]]; then
      BIOBIN_TEST_ARG="--test $(IFS="," ; echo "${regression_types[*]}")"
    fi
  fi
  BIOBIN_COVAR_ARG=""
  if [[ -n "$covariate_cov" ]]; then
    dx download "$covariate_cov" -o input/input.covariate
    BIOBIN_COVAR_ARG="--covariates input/input.covariate"
  fi
  BIOBIN_WEIGHT_ARG=""
  if [[ -n "$weight_wgt" ]]; then
    dx download "$weight_wgt" -o input/input.weight
    BIOBIN_WEIGHT_ARG="--weight-file input/input.weight"
  fi
  BIOBIN_REGION_ARG=""
  if [[ -n "$region_rgn" ]]; then
    dx download "$region_rgn" -o input/input.region
    BIOBIN_REGION_ARG="--region-file input/input.region"
  fi
  BIOBIN_INCLUDE_REGION_ARG=""
  if [[ -n "$include_region_gen" ]]; then
    dx download "$include_region_gen" -o input/input.gene
    BIOBIN_INCLUDE_REGION_ARG="--include-region-file input/input.gene"
  fi
  VCF_FILE="input/input.vcf.gz"
  TBI_FILE="$VCF_FILE.tbi"
  dx download "$variants_vcfgz" -o "$VCF_FILE"
  tabix -p vcf "$VCF_FILE" 2>&1 | tee -a output.log

  # fetch the executable(s) and shared resource file(s)
  mkdir -p shared
  dx download "$loki_db" -o shared/loki.db

  # run biobin
  mkdir biobin
  biobin \
    --threads "$NUM_CORES" \
    --settings-db shared/loki.db \
    --vcf-file "$VCF_FILE" \
    $BIOBIN_ROLE_ARG \
    $BIOBIN_PHENO_ARG \
    $BIOBIN_COVAR_ARG \
    $BIOBIN_TEST_ARG \
    $BIOBIN_WEIGHT_ARG \
    $BIOBIN_REGION_ARG \
    $BIOBIN_INCLUDE_REGION_ARG \
    --report-prefix "biobin/$output_prefix" \
    $biobin_args \
    2>&1 | tee -a output.log
  ls -laR biobin

  RGX="--force-all-control[[:space:]]+(y|Y)"
  ALL_CONTROL=""
  if [[ $biobin_args =~ $RGX ]]; then
    ALL_CONTROL="--all-control"
  fi

  # run summary script
  python /usr/bin/biobin-summary.py \
    --prefix="biobin/$output_prefix" \
    $ALL_CONTROL \
    > "biobin/${output_prefix}-summary.tsv"

  # permute, if enabled
  if [[ "$permu_count" -gt 0 ]]; then
    mkdir permu
    for p in $(seq 1 $permu_count) ; do
      mkdir permu/$p
      biobin-permute-pheno.py input/input.phenotype $p > permu/$p/input.phenotype
      biobin \
        --threads "$NUM_CORES" \
        --settings-db shared/loki.db \
        --vcf-file "$VCF_FILE" \
        $BIOBIN_ROLE_ARG \
        --phenotype-file permu/$p/input.phenotype \
        $BIOBIN_COVAR_ARG \
        $BIOBIN_TEST_ARG \
        $BIOBIN_WEIGHT_ARG \
        $BIOBIN_REGION_ARG \
        $BIOBIN_INCLUDE_REGION_ARG \
        --report-prefix "permu/$p/output" \
        $biobin_args \
        2>&1 | tee -a permu/$p/output.log
      python /usr/bin/biobin-summary.py \
        --prefix=permu/$p/output \
        $ALL_CONTROL \
        > permu/$p/output-summary.tsv
    done
    ls -laR permu
    biobin-permute-collate.py \
      "biobin/${output_prefix}-summary.tsv" \
      "permu/%d/output-summary.tsv" \
      $p \
      > "permu/${output_prefix}-permute-summary.tsv"
    permu_file=$(dx upload permu/*-permute-summary.tsv --brief)
    dx-jobutil-add-output permu_file --class="file" "$permu_file"
  fi

  # upload output files
  for f in biobin/*-bins.csv ; do
    bin_file=$(dx upload "$f" --brief)
    dx-jobutil-add-output bins_files --class="array:file" "$bin_file"
  done

  summary_file=$(dx upload biobin/*-summary.tsv --brief)
  dx-jobutil-add-output summary_file --class="file" "$summary_file"

  locus_file=$(dx upload biobin/*-locus.csv --brief)
  dx-jobutil-add-output locus_file --class="file" "$locus_file"

  mv output.log "${output_prefix}.log"
  log_file=$(dx upload "${output_prefix}.log" --brief)
  dx-jobutil-add-output log_file --class="file" "$log_file"

  for f in biobin/*-unlifted.csv ; do
    if [[ -f "$f" ]]; then
      unlifted_file=$(dx upload "$f" --brief)
      dx-jobutil-add-output unlifted_files --class="array:file" "$unlifted_file"
    fi
  done

}
