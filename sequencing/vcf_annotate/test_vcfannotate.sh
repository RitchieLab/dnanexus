#!/usr/bin/env bash
export input_vcf_gz=
export input_vcf_tbi=
export DX_RESOURCES_ID=
src/vcf_annotate.sh || exit 1
