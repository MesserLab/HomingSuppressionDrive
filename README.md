#### Data and code for Wang E, Metzloff M, Langmüller AM, Clark AG, Messer PW, Champer J., A homing suppression gene drive with multiplexed gRNAs maintains high drive conversion efficiency and avoids functional resistance alleles, bioRxiv, May 2021.

---

### HSD_AML_20210305_code_call-ML.R:
R-code to quantify the fitness costs of the gene drive from the observed frequency trajectories in our D. melanogaster cage populations. Uses functions implemented in HSD-AML-20210305-ML-v2.R

### HSD_AML_20210519_code_call-ML-r1.R:
R-code to estimate r1 resistance allele formation rate from the observed frequency trajectories in our D. melanogaster cage populations. Uses functions implemented in HSD-AML-20210519-ML-v5.R

### HSD-AML-20210305-ML-v2.R:
R-code of previously developed maximum likelihood framework (https://doi.org/10.1534/genetics.118.301893) modified to model a homing suppression drive.

### HSD-AML-20210519-ML-v5.R:
R-code of previously developed maximum likelihood framework (https://doi.org/10.1534/genetics.118.301893) modified to model a homing suppression drive & estimate r1 resistance allele formation rate.

### ../data/:
Data used by the maximum likelihood framework. HSD_AML_20210201_data_trans-matrix-FILLED.csv contains expected offspring count for each genotype combination (2 independent loci (drive/off-target)). HSD_AML_20210518_data_trans-matrix-FILLED contains expected offspring count for each genotype combination (1 drive locus, 4 alleles: wild type, drive, r1, r2). raw_cage[123].txt contain the raw genotype counts of our D. melanogaster cage populations.

### !NOTE!
Cage 1 in the manuscript = cage 2 in the R-code.
Cage 2 in the manuscript = cage 1 in the R-code.
Control cage in the manuscript = cage 3 in the R-code.
