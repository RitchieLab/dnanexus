# PLATO Single Variant Analysis

## What does this app do?
This app performs genetic associations using PLATO (the PLatform for the Analysis, Translation and Organization of large-scale data). PLATO is a software tool equipped to handle whole-genome sequence data for hundreds of thousands of samples to explore complexity and genetic architecture underlying common traits using phenome-wide association studies (PheWAS), Genome-wide association studies (GWAS), genetic interactions, environment-wide association studies (EWAS); gene-environment interactions (GxE), and copy number and rare variant analyses.

Although PLATO is a platform to perform different kinds of analyses, this app specifically focuses on single outcome (GWAS) and multiple outcome (PheWAS) analyses.

## Where can I find PLATO software manual?
You can find the detailed PLATO software manual here: https://ritchielab.org/software/plato-download

## What data are required for this app to run?
There are two input options required to run the app:
- **Genotype Files**: It is required to input genotype files in binary PLINK format (`.bed`,`.bim`,`.fam`). Please refer to the [PLINK documentation](http://zzz.bwh.harvard.edu/plink/) for information on file format.
- **Phenotype File**: An input phenotype file must be a tab-delimited, where first row indicates column headers and each subsequent row represents values for an independent sample. The first two columns must be the Family ID (FID) and Individual ID (IID) for the sample as in PLINK .FAM file, and subsequent columns contain the numeric/categorical values for each phenotype.

    | FID | IID | Pheno 1 | Pheno 2 |
    | ------ | ------ | ------ | ------ |
    | 1 |   1   |   1   |   23.2    |
    | 2 |   2 |   0   |   10.2    |
    | 3 |   3 |   0   |   9.2    |
    | 4 |   4 |   1   |   5.3    |

Other optional inputs:
  - Continuous Covariates: A tab-separated file with columns: `FID`,`IID`,`Cov1`,`Cov2`,`....`
    - Here values of the covariate should be numerical continuous.
  - Categorical Covariates: A tab-separated file with columns: `FID`,`IID`,`Cat_Cov1`,`Cat_Cov2`,`....`.
    - Here entry of non-numeric classes are allowed and will automatically generate a dummy encoding.
  - Sample File: User-specified list of sample ID's to include in the analysis. The file format should be a space/tab-delimited text file with family IDs (`FID`) in the first column and individual IDs (`IID`) in the second column.
  - SNP file: User-specified list of SNPs to include in the analysis. The file format should be either one `RSID` per line or positions in rage format (`chr`,`start`,`stop`,`label`)

## Options

### REGRESSION OPTIONS
**Regression**:
PLATO app can perform various kind of regression methods:
  - **Regress-auto** : It is default regression method for the app. It will perform either a linear or logistic regression based on the phenotype provided to the app. It works well when you have quantitative and categorical outcomes in the input phenotype file.
  - **Linear**: This option performs a simple (ordinary least-squares) regression between SNP and outcome/s.
  - **Logistic**: This option performs logistic regression between SNP and outcome/s. The logistic regression is performed using a Newton-Rhapson iteratively reweighted least squares algorithm.
  - **Firth**:  This option perform Firth bias reduction to the models being tested (for details see Firth, 1993 and Heinze and Shempel, 2002). This option can be helpful when using many categorical covariates as control variables. Also, this can be helpful when testing a large number of dichotomous phenotypes with widely varying case/control ratios, as you may see in a PheWAS.

**Association Type**:
  - **PheWAS**: Using this option, app will perform association test iteratively on all traits/outcomes provided in the input Phenotype file. It is same as  `--phewas` command in PLATO (See PLATO Manual)
  - **GWAS**: Use this option to perform association test on a single outcome/trait like in GWAS analysis. When this option is used, user need to specify column name of outcome/trait in *Outcome* option of the app.

**Outcome**: Provide column name of trait/outcome specifed in phenotype file (Only required if Assocaition Type: GWAS)

**Covariates**: Provides a comma-separated list of traits column names to include as covariates in each model.

**Lowmem**: If this argument is given, PLATO will use temporary files to output some of the results in an attempt to save memory. This can be a useful option when testing millions of models, but it may increase runtime due to the slow disk speeds when compared to memory. The DEFAULT is set to *True*.

**Enter value used to denote missingness in phenotype or covariate input files**: This option is required. **"NA"** is most commonly used to denote missingness in the phenotype data.

**Correction Type**: Select  multiple test correction strategies to use in reporting results from the the list. Currently implemented strategies are limited to Bonferroni and FDR.


### PARALLEL JOB OPTIONS
**Phenotype per job**: Provide the number of phenotype to run in each sub job. This option will split the analysis into multiple cloud instances and each instance will run models for number of phenoypte provided.
For example, the input phenoytpe file contain 10 phenotypes and user provide **"2"** in *Phenoypte per job* option, then app will run five different subjobs with 2 phenotype each and once all the proccesses are complete, it will gather all the result into single output.

It is reccomended to use this option when running PheWAS as it will exponentially reduce the runtime.

### Advanced Options
**`WARNING: If used, leave above "Regression Options" EMPTY. Please refer to documentation in PLATO manual prior to using this option.`**

**PLATO command-line options**: Command-line regression and global options that will be supplied directly to the PLATO. For example, you can run analyses such as `logistic` and `linear` with options `--permutation` or `--interactions`.

Here are some examples of the Advanced options provided by PLATO:
  - **Permutation**:
  ```
  logistic --phewas --covariates sex,age,pc1,pc2,pc3,pc4 --permutations 10000 --output out.txt
  ```
  - **Interactions**:
  ```
  logistic --interaction --pairwise --output out.txt
  ```
  - **Genomic Inflation Factor**: This option enables the calculation and reporting of the genomic inflation factor (GIF).
  ```
  logistic --inflation --inflation-method GC --output out.txt
  ```

User don't need to specify `load-data` or `load-trait` commands here since app still automate that step and provide those input file to PLATO

## FILTER OPTIONS
**MAF Threshold**: It will exculde the marker if they do not meet the minor allele frequency thresholds provided in this option.

**Case Threshold**: It will exculude the assocaition results if they do not meet the threshold of number cases provided in this option. It is only applicable to case-control PheWAS analysis.

## Common
**Ouput Filename**: The filename of the output of the regression. If not provided, default file name `output.txt` will be used.

## What does this app output?

The output of the regression is given as a column-separated file, with one line per model tested, sorted by p-value:

  - **Outcome**: Outcome name (if not using the default phenotypic status).
  - **Var1_ID**: IDs for the variables(SNP) of interest.
  - **Var1_Pos**: Chromosome and basepair location of the variant.
  - **Var1_Allele**: Allele encoded for the model.
  - **MAF**: Allele frequency of the encoded allele. Note that his value may be different for the same marker in different models because of different patterns of missingness.
  - **Num_NonMissing**: The number of non-missing samples for the model.
  - **Analysis_Type**: Type of regression performed for each model.
  - **Num_Cases**: Number of samples with case status. It will be set to "0" for linear regression models
  - **Num_Iter**: Number of least-squares steps used before declaring nonconvergence of the method
  - **Converged:** If model was converged. 1=YES;0=NO. Estimates from the model that did not converge are questionable (only relevant to the logistic regression output)
  -  **Raw_LRT_pval**: Likelihood ratio test p-value (Full model - Reduced model)
  -  **Var1_Pvalue**: Model p-value
  -  **Var1_Beta**: Model direction of effect
  -  **Var1_SE**: Model Standard Error
  -  **Overall_Pval_(LRT)**: Same as Raw LRT p-value. But, sometime in PLATO when a model does not converge it is set to 1 and user can refer to raw LRT p-value but in case of non-convergence the estimates are questionable.
