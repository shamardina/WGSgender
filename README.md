# Inferring gender from WGS data

[call_gender.py](call_gender.py) uses the number of reads per chromosome (based on `samtools idxstats`) and genotypes of SNPs (based on `bcftools query`) in the non-PAR of X chromosome to infer samples gender.

The method is described in the Supplementary Material of [Whole-genome sequencing of patients with rare diseases in a national health system](https://doi.org/10.1038/s41586-020-2434-2). Please cite this paper if you use this code in any form:

Turro, E., Astle, W.J., Megy, K. *et al.* Whole-genome sequencing of patients with rare diseases in a national health system. *Nature* **583**, 96-102 (2020). [https://doi.org/10.1038/s41586-020-2434-2](https://doi.org/10.1038/s41586-020-2434-2)

Version matching the paper: [v1.0.0-nature](https://github.com/shamardina/WGSgender/releases/tag/v1.0.0-nature).

Samtools and bcftools are required to run the script.

The paths to BAM and VCF files together with sample ID and supposed gender should be provided in the [manifest file](input/manifest.txt). Example configuration can be found in [default.conf](default.conf).
