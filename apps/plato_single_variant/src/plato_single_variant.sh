#!/bin/bash
# plato_single_variant 0.0.1
# Generated by dx-app-wizard.
#

set -e -x -o pipefail
main() {
  echo "Value of input_plink_binary: '${input_plink_binary[@]}'"
  echo "Value of input_phenotype: '$input_phenotype'"
  echo "Value of plato_analysis_string: '$plato_analysis_string'"
  echo "Value of missingness: '$missingness'"
  echo "Value of input_continuous_covariate: '$input_continuous_covariate'"
  echo "Value of input_categorical_covariate: '$input_categorical_covariate'"
  echo "Value of input_samples: '$input_samples'"
  echo "Value of input_markers: '$input_markers'"
  echo "Value of maf_threshold: '$maf_threshold'"
  echo "Value of Association Type: '${association_type[@]}'"
  echo "Value of Phenotype per job: '$split_phenotype'"
  echo "Value of Phenotype per job: '$correction'"

  # Process Plink input files
  for i in ${!input_plink_binary[@]}
  do
    name=$(dx describe "${input_plink_binary[$i]}" --name)
    if [[ "${name}" =~ \.bed$ ]]; then
      bed_file="${input_plink_binary[$i]}"
    elif [[ "${name}" =~ \.bim$ ]]; then
      bim_file="${input_plink_binary[$i]}"
    elif [[ "${name}" =~ \.fam$ ]]; then
      fam_file="${input_plink_binary[$i]}"
    fi
  done

  if [ ${split_phenotype} -gt 0 ]; then
    dx download "${input_phenotype}" -o input_phenotype
    # Get column numbers from the phenotype file
    head -n 1 input_phenotype | sed 's/ /\t/g' | tr '\t' '\n' | awk '{ print FNR "\t" $0 }' | cut -f1 | tail -n+3  > pheno_col_index
    # Split the phenotypes
    split -l ${split_phenotype} -a 3 -d pheno_col_index pheno_job

    postprocess_arg=""
    for i in pheno_job*; do
      process_jobs[$i]=$(dx-jobutil-new-job plato_reg -ibed_file="${bed_file}" \
        -ibim_file="${bim_file}" \
        -ifam_file="${fam_file}" \
        -iinput_phenotype="${input_phenotype}" \
        -iinput_continuous_covariate="${input_continuous_covariate}" \
        -iinput_categorical_covariate="${input_categorical_covariate}" \
        -iregression="${regression}" \
        -ioutcome="${outcome}" \
        -imissingness="${missingness}" \
        -iinput_samples="${input_samples}" \
        -iinput_markers="${input_markers}" \
        -imaf_threshold="${maf_threshold}" \
        -iassociation_type="${association_type}" \
        -ipheno_col=$(cat $i | tr '\n' ',' | sed 's/,$//g') \
        -icovariates="${covariates}" \
        -ioutput_filename="${output_filename}" \
        -imem="${mem}" \
        -iplato_analysis_string="${plato_analysis_string}" \
        -icase_threshold="${case_threshold}" \
        -icorrection="${correction}")
      postprocess_arg="$postprocess_arg -iplato_out_files=${process_jobs[$i]}:plato_out -iplato_log_files=${process_jobs[$i]}:plato_log -iplato_other_files=${process_jobs[$i]}:plato_other"
    done
  else
    process_jobs[$i]=$(dx-jobutil-new-job plato_reg -ibed_file="${bed_file}" \
      -ibim_file="${bim_file}" \
      -ifam_file="${fam_file}" \
      -iinput_phenotype="${input_phenotype}" \
      -iinput_continuous_covariate="${input_continuous_covariate}" \
      -iinput_categorical_covariate="${input_categorical_covariate}" \
      -iregression="${regression}" \
      -ioutcome="${outcome}" \
      -imissingness="${missingness}" \
      -iinput_samples="${input_samples}" \
      -iinput_markers="${input_markers}" \
      -imaf_threshold="${maf_threshold}" \
      -iassociation_type="${association_type}" \
      -ipheno_col=$(cat $i | tr '\n' ',' | sed 's/,$//g') \
      -icovariates="${covariates}" \
      -ioutput_filename="${output_filename}" \
      -imem="${mem}" \
      -iplato_analysis_string="${plato_analysis_string}" \
      -icase_threshold="${case_threshold}" \
      -icorrection="${correction}")
    postprocess_arg="$postprocess_arg -iplato_out_files=${process_jobs[$i]}:plato_out -iplato_log_files=${process_jobs[$i]}:plato_log -iplato_other_files=${process_jobs[$i]}:plato_other"
  fi

  postprocess=$(dx-jobutil-new-job postprocess $postprocess_arg -ioutfile:string=$output_filename -icorrection=${correction} --depends-on ${process_jobs[@]})

  dx-jobutil-add-output output_files "$postprocess:output_files" --class=jobref
}

plato_reg() {

  export LD_LIBRARY_PATH=/usr/lib:/usr/local/lib:$LD_LIBRARY_PATH

  INPUTDIR=$(mktemp -d)
  OUTPUTDIR=$(mktemp -d)

  cd "${INPUTDIR}"

  dx download "$bed_file" -o input_plink.bed
  dx download "$bim_file" -o input_plink.bim
  dx download "$fam_file" -o input_plink.fam

  # Download Phenotype files
  dx download "$input_phenotype" -o input_phenotype

  #Parallelize by phenotype
  if [[ $(echo ${pheno_col} | sed 's/,/ /g' | wc -w) -gt 0 ]];then
    cut -f1,2,"${pheno_col}" input_phenotype > input_phenotype_subset
    pheno_string="$INPUTDIR/input_phenotype_subset"
  else
    pheno_string="$INPUTDIR/input_phenotype"
  fi

  # Download Continous Covariate files
  if [ -n "$input_continuous_covariate" ];then
    dx download "$input_continuous_covariate" -o input_continuous_covariate
    load_cont="--file $INPUTDIR/input_continuous_covariate"
  fi

  # Download Categorical Covariate Files
  if [ -n "$input_categorical_covariate" ];then
    dx download "$input_categorical_covariate" -o input_categorical_covariate
    load_cat="load-categorical --file $INPUTDIR/input_categorical_covariate --missing $missingness --extra-samples"
  fi

  # Download sample file if provided
  if [ -n "$input_samples" ]; then
    dx download "$input_samples" -o input_samples
    plinkargs=" --keep $INPUTDIR/input_samples"
  fi

  #Download marker file if provided
  if [ -n "$input_markers" ];then
    dx download "$input_markers" -o input_markers
    NF=$(head -n1 input_markers | wc -w)

    # Check the format marker files, If its RSID format then use --extract else for range format use --extract range
    if [[ "$NF" == 1 ]];then
      plinkargs="$plinkargs --extract $INPUTDIR/input_markers"
    else
      plinkargs="$plinkargs --extract range $INPUTDIR/input_markers"
    fi
  fi

  # MAF Threshold
  if [ -n "$maf_threshold" ];then
    plinkargs="$plinkargs --maf $maf_threshold"
  fi

  # perform any filtering with plink as specified by sample, marker lists
  # or maf threshold
  cd

  PLINK_CMD="plink2 "
  if [ -n "$plinkargs" ]; then
    mv $INPUTDIR/input_plink.bed orig.bed
    mv $INPUTDIR/input_plink.bim orig.bim
    mv $INPUTDIR/input_plink.fam orig.fam
    plinkargs="$plinkargs --make-bed --out $INPUTDIR/input_plink"
    $PLINK_CMD --bfile orig $plinkargs
  fi
  # create output directory.  Everything in this directory should
  # be returned after job is done.
  cd $OUTPUTDIR

  plato_analysis_string=${plato_analysis_string}

  #Check if command-line string provided
  if [[ -z "$plato_analysis_string" ]]; then
    # Regression Type
    if [[ -n "$regression" ]]; then
      if [[ "$regression" == "firth" ]];then
        regression="logistic --firth"
      fi
      # Plato memory option. Use --lowmem by default
      if [ "$mem" == "true" ];then
        plato_analysis_string=" $regression --lowmem"
      else
        plato_analysis_string=" $regression"
      fi
    fi

    # Any covariates
    if [[ -n "$covariates" ]];then
      plato_analysis_string="$plato_analysis_string --covariates $covariates"
    fi

    # Type of analysis in PLATO
    for i in ${association_type[@]}
    do
      if [[ "$i" == "PheWAS" ]];then
        plato_analysis_string="$plato_analysis_string --phewas"
      elif [[ "$i" == "GWAS" ]];then
        plato_analysis_string=" $plato_analysis_string --outcome $outcome"
      fi
    done

    # Bonferoni or FDR correction
    if [[ -n "$correction" ]];then
      correction_val=$(echo ${correction[*]} | tr ' ' ',')
      plato_analysis_string="$plato_analysis_string --correction $correction_val"
    else
      plato_analysis_string="$plato_analysis_string"
    fi

    # If output file name provided
    if [[ -n $output_filename ]];then
      plato_analysis_string=" $plato_analysis_string --output $output_filename"
      outfile="$output_filename"
    else
      plato_analysis_string=" $plato_analysis_string --output output.txt"
      outfile="output.txt"
    fi
  else
    plato_analysis_string="$plato_analysis_string"
    outfile=$(echo $plato_analysis_string | sed 's/^.*--output //g' | sed 's/ .*$//g')
  fi

  # Use all the core in an AWS instance
  if ! [[ "$plato_analysis_string" =~ --threads ]];then
    threads="--threads $(nproc)"
    analysis2=${plato_analysis_string/linear /linear $threads }
    analysis3=${analysis2/regress-auto /regress-auto $threads }
    analysis4=${analysis3/logistic /logistic $threads }
    plato_analysis_string=$analysis4
  fi

  # Plato command
  plato load-data \
    --bed $INPUTDIR/input_plink.bed \
    --bim $INPUTDIR/input_plink.bim \
    --fam $INPUTDIR/input_plink.fam \
    recode-alleles --auto \
    load-trait \
    --extra-samples \
    --missing $missingness \
    --file $pheno_string \
    $load_cont \
    $load_cat \
    $plato_analysis_string || touch ${outfile}


  # Filter out the results with MAF or Case filter provided.

  mkdir temp

  if [ -n "$maf_threshold" ] && [ -n "$case_threshold" ]; then
    awk 'BEGIN{OFS=FS="\t"}{gsub(":","\t",$4)}1' ${outfile}  | sed 's/Var1_MAF/Var1_Allele\tMAF/g' | awk -v var1="$maf_threshold" 'BEGIN{OFS=FS="\t"}{if($5>var1) print $0}' > temp/${outfile}_maf_threshold
    #	cat temp/${outfile}_maf_threshold
    case_col=$(head -n1 temp/${outfile}_maf_threshold | tr '\t' '\n' | grep -n "Num_Cases" | cut -d":" -f1)
    awk -v var1="$case_col" -v var2="$case_threshold" 'BEGIN{OFS=FS="\t"}{if($var1>=var2)print $0}' temp/${outfile}_maf_threshold > $outfile
  elif [ -n "$maf_threshold" ]; then
    awk 'BEGIN{OFS=FS="\t"}{gsub(":","\t",$4)}1' $outfile | sed 's/Var1_MAF/Var1_Allele\tMAF/g' | awk -v var1="$maf_threshold" 'BEGIN{OFS=FS="\t"}{if($5>var1) print $0}' > temp/${outfile}_maf_threshold
    cp temp/${outfile}_maf_threshold ${outfile}
  elif [ -n "$case_threshold" ]; then
    awk 'BEGIN{OFS=FS="\t"}{gsub(":","\t",$4)}1' ${outfile} | sed 's/Var1_MAF/Var1_Allele\tMAF/g' > temp/${outfile}_maf
    case_col=$(head -n1 temp/${outfile}_maf | tr '\t' '\n' | grep -n "Num_Cases" | cut -d":" -f1)
    awk  -v var1="$case_col" -v var2="$case_threshold" 'BEGIN{OFS=FS="\t"}{if($var1>=var2)print $0}' temp/${outfile}_maf > ${outfile}
  fi
  rm -rf temp

  mkdir $OUTPUTDIR/out
  mv $outfile $OUTPUTDIR/out/

  mkdir $OUTPUTDIR/log
  mv *.log $OUTPUTDIR/log/

  if [ $(find . -maxdepth 1 -type f | wc -l) -gt 0 ]; then
    mkdir $OUTPUTDIR/other
    mv $(find . -maxdepth 1 -type f) $OUTPUTDIR/other/
  fi

  cd

  # The following line(s) use the utility dx-jobutil-add-output to format and
  # add output variables to your job's output as appropriate for the output
  # class.  Run "dx-jobutil-add-output -h" for more information on what it
  # does.
  out=$(dx upload --brief $OUTPUTDIR/out/*)
  dx-jobutil-add-output plato_out "${out}" --class=file

  log=$(dx upload --brief $OUTPUTDIR/log/plato.log)
  dx-jobutil-add-output plato_log "${log}" --class=file

  if [ -d "$OUTPUTDIR/other/" ]; then
    for i in $OUTPUTDIR/other/*
    do
      other=$(dx upload --brief $i)
      dx-jobutil-add-output plato_other "${other}" --class=array:file
    done
  else
    touch other.txt
    other=$(dx upload --brief other.txt)
    dx-jobutil-add-output plato_other "${other}" --class=file
  fi
  #sleep 36000
}

postprocess(){

  OUTPUTDIR=$(mktemp -d)
  mkdir $OUTPUTDIR/output_files

  for i in "${!plato_out_files[@]}"; do
    out_name=$(dx describe "${plato_out_files[$i]}" --name)
    dx download "${plato_out_files[$i]}" -o plato_out-$i
    mv plato_out-$i $OUTPUTDIR/
  done

  for i in "${!plato_log_files[@]}"; do
    dx download "${plato_log_files[$i]}" -o plato_log-$i
    mv plato_log-$i $OUTPUTDIR/
  done

  for i in "${!plato_other_files[@]}"; do
    name=$(dx describe "${plato_other_files[$i]}" --name)
    dx download "${plato_other_files[$i]}" -o $name
    if [ -s "$name" ]; then
      echo "File is not empty! Moving to output_files/"
      mv $name $OUTPUTDIR/output_files/
    else
      echo "File is empty! Removing empty file"
      rm $name
    fi
  done


  cd $OUTPUTDIR

  out_header=$(head -n1 plato_out-0 | sed 's/\t/\\t/g')

  pval_col=$(head -n1 plato_out-0 | tr '\t' '\n' | grep -n "Overall_Pval_(LRT)" | cut -d":" -f1)

  if [ "$correction" == "Bonferoni" ]; then
    sort -s -m -gk$pval_col,$pval_col $(ls plato_out-*) | grep -v "Var1" \
      | cut -f1-$pval_col \
      | awk 'BEGIN{FS=OFS="\t"}{print $0, $var1*var2}' \
      var1=$pval_col \
      var2=$(ls plato_out-* | xargs -n 1 tail -n+2  | wc -l) \
      | awk 'BEGIN{FS=OFS="\t"}{if($var3<1) print $0,$var4=$var3; else print $0,$var4=1}'  \
      var3=$(echo "$pval_col" + 1 | bc -l) \
      var4=$(echo "$pval_col" + 2 | bc -l) \
      | cut --complement -f$(echo "$pval_col" + 1 | bc -l) \
      | sed "1i $(echo ${out_header} | sed -E 's/Overall_Pval_adj_Bonferroni|Overall_Pval_adj_FDR//g')Overall_Pval_adj_Bonferroni" \
      > output_files/${out_name}
  elif [ "$correction" == "FDR" ]; then
    pmax=1
    iter=0
    num_test=$(ls plato_out-* | xargs -n 1 tail -n+2  | wc -l)

    paste -d "\t" <(sort -s -m -gk$pval_col,$pval_col $(ls plato_out-*) \
      | grep -v "Var1" | cut -f1-$pval_col) \
      <(for i in $(sort -s -m -gk$pval_col,$pval_col $(ls plato_out-*) \
      | grep -v "Var1" | cut -f$pval_col | tac); \
    do pcal=$(echo "$(echo $i | sed -e 's/[eE]+*/\*10\^/')*${num_test}/(${num_test}-${iter})" | bc -l ); \
      pmax=$(if [ $(echo "${pcal} < $(echo ${pmax} | sed -e 's/[eE]+*/\*10\^/' | bc -l)" | bc -l) == 1 ]; \
    then echo $pcal ; \
    else echo ${pmax}; \
    fi);iter=$(($iter+1)); \
    pmax=$(printf %e ${pmax}); echo $pmax ; \
  done | tac) | sed "1i $(echo ${out_header} \
  | sed -E 's/Overall_Pval_adj_Bonferroni|Overall_Pval_adj_FDR//g')Overall_Pval_adj_FDR" > output_files/${out_name}
  else
    sort -s -m -gk$pval_col,$pval_col $(ls plato_out-*) | grep -v "Var1" | sed "1i ${out_header}" > output_files/${out_name}
  fi

  cat plato_log-* > output_files/${out_name}.log

  for i in $OUTPUTDIR/output_files/*; do
    out=$(dx upload --brief $i)
    dx-jobutil-add-output output_files "$out" --class=array:file
  done
}
