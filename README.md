# Diablo ACMG
Tool for automated classification of genetic variants according to ACMG criteria (implemented with 
[open-cravat](https://github.com/KarchinLab/open-cravat)). This tool takes VCF file as input, annotates variants with 
several [modules](#required-modules), assigns ACMG criteria to them and return annotated file in TSV format. 

## Table of Contents

1. [Installation](#installation)
2. [Usage](#Usage)
3. [Required modules](#required-modules)

## Installation

I am currently working on deploying Diablo ACMG as conda package.

For now, the easiest way to install Diablo annotate is to git clone repository and install dependencies from .yml file.

```
git clone https://github.com/TheSergeyPixel/Diablo_ACMG.git
conda env create -f /path/to/cloned/repo/diablo_annotate.yml
```
Make sure, that you've installed all the requirements for open-cravat (check [Required modules](#required-modules) section)

## Usage 
Basic usage:

```
python Diablo_annotate.py -i /path/to/your/vcf.gz -o /output/file/name -d /path/to/db
```
You can also test your installation via running Diablo ACMG on test sample from core directory: <br/>
```
python Diablo_annotate.py -i test/IDTEPF.vcf.gz -o test/test.tsv
```
It's important to note, that **we don't support multi-sample VCFs** yet. 

Description of options:

| **Option** | **Description**                                                                                                                           |
|------------|-------------------------------------------------------------------------------------------------------------------------------------------|
| ```-i```   | Path to your input VCF file or TSV annotated file <br/>(latter case assumes, that you have already annotated your file with open-cravat). |
| ```-o```   | Path to your '.tsv' output file                                                                                                           |
| ```-s```   | You can provide a number of lines (variants),<br/>that you would like to see in your sorted output.                                       |
| ```-d```   | Path to the **folder** with files required for annotation<br/>(default is 'db', located in the repo).                                     |
| ```-t```   | Number of parallel processes. **NOTE**: Has high RAM requirements.                                                                        |
| ```-S```   | If 'True', Create a separate file with predicted splice variants (SpliceAI > 0.7)<br/> in the same directory as output.                   |



## Required Modules

| **Module name**  |
|------------------|
| ClinVar          |
| dbSNP            |
| Fathmm           |
| Fathmm_mkl       |
| GenoCanyon       |
| GERP             |
| GnomAD v3.1.1    |
| InterPro         |
| LRT              |
| MetaLR           |
| MetaSVM          |
| MutationAssessor |
| MutationTaster   |
| OMIM             |
| Polyphen2        |
| Provean          |
| SIFT             |
| Siphy            |
| SpliceAI         |
| HPO              |

These modules can be installed with the following command (make sure 
[open-cravat](https://github.com/KarchinLab/open-cravat) is installed):
```
oc module install clinvar dbsnp fathmm fathmm_mkl genocanyon gerp gnomad3 interpro lrt metalr metasvm mutation_assessor
mutationtaster omim polyphen2 provean sift siphy spliceai hpo
```
Make sure, that you have all the base open-cravat modules installed.
```
oc module install-base
```
Feel free to report any issues, we are still testing our software. Any suggestions would be highly appreciated :) <br/> 
<br/>
Due to complicated installation of open-cravat and large size of its modules, we will try to create a docker image, 
that is going to use open-cravat servers for annotation.
