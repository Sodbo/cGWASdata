# A network-based conditional genetic association analysis of the human metabolome

The repository stores all code and data necessary for reproducing all results of cGAS paper 
(doi:https://doi.org/10.1101/096982).

## Data description

The results of the paper were obtained from analysis of GWAS summary data 
(for more details see paper). **Data** folder stores GWAS summary data, 
SNP information file and matrices with pairwise correlation coefficients 
between metabolite's measurements. **Code** folder stores R functions for 
applying cGWAS to summary gwas data. 

### GWAS summary data

All gwas stored in text files with headings:

|SNP|beta|se |Z  |P  |
|---|----|---|---|---|
|rs11224105|0.17656|0.03382|5.21984|1.79e-07|


where **SNP** is the marker name; **beta** - effect size;  **se** - standard error of the beta; 
<br> **Z** is $Z_{value}=\frac{beta}{se}$; **P** is P-value of association.
<br>Sample size is N=1785.


Information about alleles, frequencies, imputation quality, physical position etc could be found in **SNP_information.txt.gz**.
<br> Heading of the **SNP_information.txt.gz** file:

|chr|SNP|pos|A1|A2|freq|R2_impute_info|hw|varg_1785|
|---|---|---|--|--|----|--------------|--|---------|
|1|rs11804171|713682|A|T|0.94829|0.911636|0.09440|0.09479|
|1|rs2977670|713754|G|C|0.05172|0.91039|0.09676|0.09470|

where **chr** is chromosome; **SNP** - SNP name; **pos** - position (r37 assembly); **A1** - effective allele; **A2** - reference allele; **freq** - effective allele frequency; **R2_imputie_info** - imputation quality; **hw** - Hardy-Weinberg equilibrium P-value; **varg_1785** - exact variance of the SNP (needed for calculation of cGWAS).

### Folders

- uGWAS:
Stores files with univariate GWAS results obtained by OmicABEL (doi:10.12688/f1000research.4867.1), filtered by P-value<1e-6. Genomic control (GC) correction wasn't applied.

- GGM-cGWAS:
Stores GGM-cGWAS results filtered by P-value<1e-6. GC correction was applied. GC Lambdas are stored in "GGM_cGWAS_gc_lambda.txt" file.

- BN-cGWAS:
Stores BN-cGWAS results filtered by P-value<1e-6. GC correction was applied. GC Lambdas are stored in "BN_cGWAS_gc_lambda.txt" file.

- uGWAS_snps_from_paper:
Stores uGWAS results for all SNPs mentioned in the paper obtained by "lm" function in R. GC correction wasn't applied.


### Matrices

- 20171207_corr_matrix.txt: Pearson correlation matrix for 151 metabolites.

- 20171207_partial_corr_matrix.txt: partial correlation matrix for 151 metabolites. We used "ppcor" R package for calculations.

- 20171207_partial_corr_pvalues_matrix.txt: partial correlation P-value matrix for 151 metabolites.

- 20171207_biochemical_distances.txt: biochemical distances used for BN-cGWAS. This matrix was produced in work of Krumsiek et.al, 2011 (Krumsiek, J., Suhre, K., Illig, T., Adamski, J., & Theis, F. J. (2011). Gaussian graphical modeling reconstructs pathway reactions from high-throughput metabolomics data. BMC Systems Biology, 5(1), 21. https://doi.org/10.1186/1752-0509-5-21).

### Scripts

Scripts are located in "scripts" folder.

- exact_cGWAS_functions.R: the main function for calculation of cGWAS using uGWAS results and correlation matrices. All descriptions are in file. 

- main_BN.R: script runs the cGAS analysis where covariates where selected based on biochemical distances between metabolites

- main_GGM.R: script runs the cGAS analysis where covariates where selected based on partial correlations between metabolites.

- Suppl.table.1A.R: script that creates Supplementary Table 1A, that compares results of uGAS and BN-cGAS

- Suppl.table.1B.R: script that creates Supplementary Table 1B, that compares results of BN-cGAS and GGM-cGAS

- Suppl.table.2.R: script that creates Supplementary Table 2, where all loci found by GGM-cGAS and uGAS are collected

- Figure_1.R: scripts reproduce Figure 1 presented in the paper

- Figure_2.R: scripts reproduce Figure 2 presented in the paper

- clumping_functions.R: fuinctions for performing clumping of GWAS results

- utilities.R: utilities function that are used by exact_cGWAS_functions


### Citation

If you use the exact_cGWAS_functions.R and the other procedures, respectively,
please cite:

Tsepilov, Y. A., Sharapov, S. Z., Zaytseva, O. O., Krumsiek, J., Prehn, C., Adamski, J., … Aulchenko, Y. S. (2018). A network-based conditional genetic association analysis of the human metabolome. Submitted
