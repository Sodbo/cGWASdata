
# Loading functions
source('../code/exact_cGWAS.R')

# Set the responce

response <- 'lysoPC.a.C20:4'

# Set the covariates

covariates <- c('lysoPC.a.C20:3')

# Loading correlation matrix for 151 metabolites

cor_mat <- read.table('../data/20171207_corr_matrix.txt')

# Loading information about SNPs, where column 'var_1785' 
# contains variance of the SNP

snp_info <- data.table::fread('../data/30_SNP_information.txt',
	data.table = FALSE)

exact_cGWAS(
	CovM = cor_matrix,
	all_varg = snp_info$var_1785,
	response = response,
	N = 1785,
	covariates = covariates,
	cn_b = "beta_SNP",
	cn_snp = "SNP",
	cn_se = "se_SNP",
	output_threshold=1,
	good_snps = NULL,
	all_CR = 1,
	correction = TRUE,
	path_uGWAS = '../data/uGWAS_snps_from_paper'
)