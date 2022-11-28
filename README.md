# Diablo ACMG
Tool for automated classification of genetic variants according to ACMG criteria (implemented with 
[open-cravat](https://github.com/KarchinLab/open-cravat)). It takes VCF file as input, annotates variants with several
[modules](), assigns ACMG criteria to them and return annotated file in TSV format. 

## Table of Contents

1. [Usage](#Usage)
2. [Required databases](#Required Modules)

## Usage 



## Required Modules

| **Module name**  |
|------------------|
|      ClinVar     |
|       dbSNP      |
|      Fathmm      |
|    Fathmm_mkl    |
|    GenoCanyon    |
|       GERP       |
|   GnomADv3.1.1   |
|     InterPro     |
|        LRT       |
|      MetaLR      |
|      MetaSVM     |
| MutationAssessor |
|  MutationTaster  |
|       OMIM       |
|     Polyphen2    |
|      Provean     |
|       SIFT       |
|       Siphy      |
|     SpliceAI     |
|        HPO       |

These modules can be installed with the following command (make sure 
[open-cravat](https://github.com/KarchinLab/open-cravat) is installed):
```
oc module install clinvar dbsnp fathmm fathmm_mkl genocanyon gerp gnomad3 interpro lrt metalr metasvm mutation_assessor
mutationtaster omim polyphen2 provean sift siphy spliceai hpo
```

