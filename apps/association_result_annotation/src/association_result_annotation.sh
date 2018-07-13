#!/bin/bash
# association_result_annotation 1.0.0
set -e -x -o pipefail

#########################################
# Performs initial checks on input file
# Arguments: None
# Returns: None
#########################################
check_inputs() {

    # check if the input file is tab separated
    if [ -z "$(sed -n '/\t/p;q' input_file )" ]
    then
        dx-jobutil-report-error "ERROR: Incorrect input file format. App only accepts tab separated input file."
    else
        input_filename=$(dx describe --name "$associations")
    fi

    # Seperate Var1_Pos from PLATO output to "chr" and "pos" columns
    if [ "$chr_col" == "Var1_Pos" ] && [ "$pos_col" == "Var1_Pos" ]
    then
        if [ "$(head -n 1 input_file | tr '\t' '\n' | grep "Var1_Pos")" == "Var1_Pos" ]
        then
            awk 'BEGIN{OFS=FS="\t"}{gsub(":","\t",$3)}1' input_file | sed 's/Var1_Pos/chr\tpos/g' > up_input_file
            mv up_input_file input_file
            chr_col="chr"
            pos_col="pos"
        else
            dx-jobutil-report-error "ERROR: Something went wrong! Please check if column names for chromosome and position are provided correctly."
        fi
    fi

}

##############################################################
# Establishes the context id where resource files are located
# Arguments: None
# Returns: None
##############################################################
establish_resource_context() {

    DX_ASSETS_ID=""
    if [[ "$DX_RESOURCES_ID" != "" ]]; then
        
        DX_ASSETS_ID="$DX_RESOURCES_ID"
    
    else

        DX_ASSETS_ID="$DX_PROJECT_CONTEXT_ID"

    fi

}

##########################
# Runs main code pipeline
# Arguments: None
# Returns: None
##########################
main() {
       
    # get the context id where to find resources 
    establish_resource_context
     
    # download input TSV file containing GWAS data
    dx download "$associations" -o input_file
    out_suffix=""

    # parse through input options
    check_inputs

    #################################################################################
    # Check Gene and GWAS annotation provided 
    # Download GRASP and GWAS database when requested
    # If Gene annotations requested, run biofilter and import the results in SQLite 
    #################################################################################

    # check if all the values are provided if any Gene and GWAS annotation options are true
    if ([ "$ebi_gwas" == "true" ] || [ "$grasp" == "true" ] || [ "$gene" == "true" ] || [ "$up_gene" == "true" ] || [ "$down_gene" == "true" ]) && ([ "$chr_col" == "" ] || [ "$pos_col" == "" ])
    then
        dx-jobutil-report-error "ERROR: One or more column names under Gene and GWAS annotation option not provided!"
    elif ([ "$ebi_gwas" == "true" ] || [ "$grasp" == "true" ] || [ "$gene" == "true" ] || [ "$up_gene" == "true" ] || [ "$down_gene" == "true" ]) && ([ "$chr_col" != "" ] || [ "$pos_col" != "" ])
    then
        # check if each column provided by user in the Gene and GWAS annotation exist within the input file.
        if [ "$(head -n 1 input_file | tr '\t' '\n' | grep "$chr_col")" == "$chr_col" ]
        then
            chr_col_num=$(head -n 1 input_file | tr '\t' '\n' | grep -n "$chr_col" | cut -d":" -f1)
        else
            dx-jobutil-report-error "ERROR: Column name provided for chromosome was not found in input file"
        fi

        if [ "$(head -n 1 input_file | tr '\t' '\n' | grep "$pos_col")" == "$pos_col" ]
        then
            pos_col_num=$(head -n 1 input_file | tr '\t' '\n' | grep -n "$pos_col" | cut -d":" -f1)
        else
            dx-jobutil-report-error "ERROR: Column name provided for basepair position was not found in input file"
        fi


        if $ebi_gwas || $grasp
        then
            dx download "$sql_database" -o anno.db

            if [ "$grasp_pval" != "" ]
				then
sqlite3 anno.db <<!
.separator \t
create table grasp_subset as select * from grasp where Pvalue < $grasp_pval;
create index "grasp_subset_chr" on "grasp_subset" ("chr","pos");
create index "grasp_subset_chr_2" on "grasp_subset" ("chr");
create index "grasp_subset_pos" on "grasp_subset" ("pos"); 
!
            fi

        fi

		if $ebi_gwas;
		then
		    gwas_join="left join ebi_gwas_pos egp on a."${chr_col}"=egp.chr and a."${pos_col}"=egp.pos left join ebi_gwas eg on eg.chr_id=egp.chr_hg38 and eg.chr_pos=egp.pos_hg38"
			gwas_query=",IFNULL(group_concat(distinct eg.DISEASE_TRAIT), 'NA') as GWAS_trait"
			out_suffix="${out_suffix}_gwas"
		fi

		if $grasp;
		then
            grasp_join="left join grasp_subset as g on a."$chr_col"=g.chr and a."$pos_col"=g.pos"
            grasp_query=",IFNULL(group_concat(distinct Phenotype), 'NA') as GRASP_Trait"
            out_suffix="${out_suffix}_grasp"
        fi

        # Create position file of SNPs from the association result file
        cut -f$chr_col_num,$pos_col_num input_file |  sort -u  > input_file_snp_position.txt

        # Run biofilter to get gene, upstream gene, downstream gene and GWAS
        sudo mkdir /usr/share/biofilter
        sudo chmod a+rwx /usr/share/biofilter
        dx download -r "$DX_ASSETS_ID:/Biofilter/2.4/*" -o /usr/share/biofilter/
        dx download "$DX_ASSETS_ID:/LOKI/LOKI-20160428-noSNPs.db" -o /usr/share/biofilter/loki.db

        python2.7 /usr/share/biofilter/biofilter.py -v -k /usr/share/biofilter/loki.db --gbv 37 -P input_file_snp_position.txt -a position gene upstream downstream
        awk 'BEGIN{OFS=FS="\t"}{gsub(":","\t",$2)}1' biofilter.position.gene-upstream-downstream | sed 's/chr//g' | sed 's/position/chr_37\tpos_37/g' > biofilter_snp_anno
        # Import biofilter result into
        biofilter_anno_col=$(head -n1 biofilter_snp_anno | \
        sed 's/^#/chr/g' | \
        sed 's/\//_/g' | \
        sed 's/%/_/g'  | \
        awk 'BEGIN{OFS=FS="\t"}{gsub("distance","up_distance",$7)}1' | \
        awk 'BEGIN{OFS=FS="\t"}{gsub("distance","down_distance",$9)}1' | \
        sed 's/\t/ varchar(255),/g' | \
        sed 's/$/ varchar(255)/g' )

        if $gene;
        then
            gene_query="left join biofilter_anno b on a.'$chr_col'=b.chr_37 and a.'$pos_col'=b.pos_37"
            gene_col=",IFNULL(group_concat(distinct gene), 'NA') as Gene"
            out_suffix="${out_suffix}_gene"
        fi


        if $up_gene; 
        then
            
            if [ "${gene_query}" != "" ]
            then
                up_gene_col=",IFNULL(group_concat(distinct b.upstream), 'NA') as 'Upstream_Gene', IFNULL(group_concat(distinct b.up_distance), 'NA') as 'Upstream_Distance'"
			    out_suffix="${out_suffix}_up-gene"
			else
				gene_query="left join biofilter_anno b on a.'$chr_col'=b.chr_37 and a."$pos_col"=b.pos_37"
				up_gene_col=",IFNULL(group_concat(distinct b.upstream), 'NA') as 'Upstream_Gene', IFNULL(group_concat(distinct b.up_distance), 'NA') as 'Upstream_Distance'"
				out_suffix="${out_suffix}_up-gene"
			fi
		
        fi

        if $down_gene; 
        then
            
            if [ "${gene_query}" != "" ]
            then
                down_gene_col=",IFNULL(group_concat(distinct b.downstream), 'NA') as 'Downstream_Gene', IFNULL(group_concat(distinct b.down_distance), 'NA') as 'Downstream_Distance'"
                out_suffix="${out_suffix}_down-gene"
            else
                gene_guery="left join biofilter_anno b on a."$chr_col"=b.chr_37 and a."$pos_col"=b.pos_37"
                down_gene_col=",IFNULL(group_concat(distinct b.downstream), 'NA') as 'Downstream_Gene', IFNULL(group_concat(distinct b.down_distance), 'NA') as 'Downstream_Distance'"
                out_suffix="${out_suffix}_down-gene"
            fi
        
        fi

        sqlite3 anno.db "create table biofilter_anno ($biofilter_anno_col)"

sqlite3 anno.db <<!
.separator \t
.import biofilter_snp_anno biofilter_anno
delete from biofilter_anno where rowid = 1;
!
    else
        echo "Skipping Gene and GWAS Annotations"
    fi
    #################################################################################


    ###########################################
    # Import input file into SQLite database 
    ###########################################

    # Get the column names from input file
    assoc_tbl_col=$(head -n1 input_file | sed 's/)//g' | sed 's/(//g' | sed 's/\t/ varchar(255),/g' | sed 's/$/ varchar(255)/g')

    # Create a table for input results file
    sqlite3 anno.db "create table assoc_result ($assoc_tbl_col)"

    # Insert input results to the table.
sqlite3 anno.db <<!
.separator \t
.import input_file assoc_result
delete from assoc_result where rowid = 1;
!
    #################################################################################

    ##################################
    # ICD-9 Description Annotation 
    ##################################

    # check if column names provided correctly else die
    if [ "$icd9_desc" == true ] && [ "$icd9_col" == "" ]
    then
        dx-jobutil-report-error "ERROR: ICD-9 code column name from the input file not provided!"
    elif [ "$icd9_desc" == true ] || [ "$(head -n 1 input_file | tr '\t' '\n' | grep "$icd9_col")" == "$icd9_col" ]
    then
        # Import ICD-9 code description table to the database
        dx download "$DX_ASSETS_ID:ICD9/icd9_codes_description.txt" -o icd9_code_desc.txt

sqlite3 anno.db <<!
create table icd9_code_desc (icd9_code varchar(8), desc text);
.header on
.separator "\t"
.import icd9_code_desc.txt icd9_code_desc
delete from icd9_code_desc where rowid = 1;
!

        #Build the query
        icd9_join="left join icd9_code_desc b on a."$icd9_col"=b.icd9_code"
        icd9_select=",IFNULL(b.desc, 'NA') as icd9_description"
        out_suffix="${out_suffix}_icd9-desc"
    else
        echo "Skipping ICD-9 description annotation.."
    fi
    #################################################################################


    ######################################
    # Odds Ratio and Confidence Interval 
    ######################################

    if $or_val;
    then
        or_val_query=",exp(var1_beta) as odds_ratio,round(exp(log(exp(Var1_beta)) - 1.96*Var1_SE),3) || ',' || round(exp(log(exp(Var1_beta)) + 1.96*Var1_SE),3) as OR_95CI"
        out_suffix="${out_suffix}_or"
    fi

    #################################################################################


    ############################
    # Case-Control Annotation 
    ############################

    if $case_control_num;
    then
        case_control_query=",Num_Cases as Cases, Num_NonMissing-Num_Cases as Controls"
        out_suffix="${out_suffix}_case-control"
    fi

    #################################################################################


    #########################################
    # Final Query for annotation requested 
    #########################################

echo "
sqlite3 anno.db <<!
.load /usr/local/lib/libsqlitefunctions.so
.headers on
.mode tabs
.output ${input_filename}${out_suffix}
select a.* ${icd9_select} ${case_control_query} ${or_val_query} ${gene_col} ${up_gene_col} ${down_gene_col} ${gwas_query} ${grasp_query} from assoc_result a ${icd9_join} ${gene_query} ${gwas_join} ${grasp_join} group by a.rowid;
!
"

sqlite3 anno.db <<!
.load /usr/local/lib/libsqlitefunctions.so
.headers on
.mode tabs
.output ${input_filename}${out_suffix}
select a.* ${icd9_select} ${case_control_query} ${or_val_query} ${gene_col} ${up_gene_col} ${down_gene_col} ${gwas_query} ${grasp_query} from assoc_result a ${icd9_join} ${gene_query} ${gwas_join} ${grasp_join} group by a.rowid;
!

    #################################################################################

    # add .tsv extension to the output file
    output_filename="${input_filename}${out_suffix}.tsv"
    mv "${input_filename}${out_suffix}" "${output_filename}"

    # Upload output back to DNAnexus
    output_tsv=$(dx upload ${output_filename} --brief)
    dx-jobutil-add-output output_tsv "$output_tsv" --class=file

}
