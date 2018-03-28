
# Loading functions
source('../code/exact_cGWAS.R')

# Set the responce

response <- 'lysoPC.a.C20.4'

# Set the covariates

covariates <- c('lysoPC.a.C20.3')

# Loading correlation matrix for 151 metabolites

cor_matrix <- read.table('../data/20171207_corr_matrix.txt')

# Loading information about SNPs, where column 'var_1785' 
# contains variance of the SNP

snp_info <- data.table::fread('../data/30_SNP_information.txt',
	data.table = FALSE)

all_varg <- snp_info$varg_1785

names(all_varg) <- snp_info$SNP

exact_cGWAS(
	CovM = cor_matrix,
	all_varg = all_varg,
	response = response,
	N = 1785,
	covariates = covariates,
	cn_b = "beta",
	cn_snp = "SNP",
	cn_se = "se",
	output_threshold = 1,
	good_snps = NULL,
	#correction = TRUE,
	all_CR = rep(1, length(all_varg)),
	path_uGWAS = '../data/uGWAS_snps_from_paper'

)