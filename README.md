# Inferring gender from WGS data

[call_gender.py](call_gender.py) uses the number of reads per chromosome (based on `samtools idxstats`) and genotypes of SNPs (based on `bcftools query`) in the non-PAR of X chromosome to infer samples gender.

The method is described in the Supplementary Material of [Whole-genome sequencing of rare disease patients in a national healthcare system](https://www.biorxiv.org/content/10.1101/507244v1). Please cite this paper if you use this code in any form.

Samtools and bcftools are required to run the script.

The paths to BAM and VCF files together with sample ID and supposed gender should be provided in the [manifest file](input/manifest.txt). Example configuration can be found in [default.conf](default.conf).
