<!-- dx-header -->
# Annotate GWAS, PheWAS Associations (DNAnexus Platform App)
<!-- /dx-header -->

## What does this app do?

App primarily annotate the GWAS or PheWAS association results with Genes, GWAS Catalog, Odds Ratio, Phenotype long description (mainly ICD-9 codes).
Also, app can annotate following irrespective of association results:

- SNP with Gene: Provide chromosome and basepair location as input file and select gene annotation options.
- SNP with EBI GWAS catalog and GRASP catalog: Provide chromosome and basepair location as input file and select GWAS catalog annotation options.
- ICD-9 code with long description: Provide list of ICD-9 codes as input file and select ICD-9 annotation option.

## What data are required for this app to run?

App requires a tab separated file.

## Usage
### Options
**ICD-9 Description**

- **ICD-9 description:** Allows to add ICD-9 code long descriptions to the
- **ICD-9 Column name:** It is required, if "ICD-9 description" is true. It allows app to know which column in input file contains ICD-9 codes.

</br>

**Gene and GWAS Catalog**

- **EBI GWAS catalog:** Add traits found in EBI GWAS catalog for SNPs position in input file.
- **GRASP:** Add traits found in GRASP GWAS catalog for SNPs position in input file. GRASP has more thorough list of published GWAS results but it also include results with p-value < 0.01 and that makes it a bigger database than EBI.
- **GRASP p-value threshold:** As stated above about p-value threshold in GRASP, we suggest to use p-value threshold of 1E-05 to reduce the search space and also provide equivalent comparision with EBI catalog. You can change the p-value to adjust search space.
- **Gene:** Map gene to chromosome-basepair position in input file.
- **Upstream Gene:** Map upstream gene to chromosome-basepair position in input file and provide distance to SNP.
- **Downstream Gene:** Map downstream gene to chromosome-basepair position in input file and provide distance to SNP.
- **Chromosome column name:** It is required if any of the option in "Gene and GWAS Catalog" is true. It allows app to know which column in input file contain chromosome numbers.
- **Basepair position column name:** It is required if any of the option in "Gene and GWAS Catalog" is true. It allows app to know which column in input file contain basepair position.

</br>

**Statistics**

- **Odds Ratio:** Calculate odd ratio from the beta and standard error in the input file
- **Case-Control Numbers:**

### Example Screenshots
**Step 1**
<img src="https://github.com/GeisingerBTI/dnanexus/tree/master/apps/association_result_annotation/assets/step1.png">
![alt text](../../assets/step1.png)
**Step 2**
![alt text](https://github.com/GeisingerBTI/dnanexus/tree/master/apps/association_result_annotation/assets/step2.png)



## What does this app output?

App outputs all the columns from input file with extra columns requested for annotation.

## How does this app work?

App run on a custom bash script that process the input file and import it into a SQLite database. Then, SQLite queries are used to add the requested annotations to the input file. It also uses [Biofilter 2.4](http://ritchielab.psu.edu/files/RL_software/biofilter-manual-2.4.pdf) to annotate SNP to Genes.
