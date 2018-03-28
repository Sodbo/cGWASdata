# Sodbo Sharapov

# Loading functions
source('../code/exact_cGWAS.R')

# Loading correlation matrix for 151 metabolites

cor_matrix <- read.table('../data/20171207_corr_matrix.txt')

# Loading information about SNPs, where column 'var_1785' 
# contains variance of the SNP and CR contains call rate

snp_info <- data.table::fread('../data/30_SNP_information.txt',
	data.table = FALSE)

all_varg <- snp_info$varg_1785

names(all_varg) <- snp_info$SNP

all_CR <- snp_info$CR

names(all_CR) <- snp_info$SNP

# Load matrix with GGM partial correlations and their P-values

ggm_mat <- read.table('../data/20171207_partial_corr_matrix.txt', 
	head = TRUE, 
	stringsAsFactors = FALSE
	)

ggm_mat_p <- read.table('../data/20171207_partial_corr_pvalues_matrix.txt', 
	head = TRUE, 
	stringsAsFactors = FALSE
	)

# Intersect GGM matrix with correlaton matrix

ggm_mat <- ggm_mat[rownames(ggm_mat) %in% rownames(cor_matrix),]

ggm_mat <- ggm_mat[,colnames(ggm_mat) %in% colnames(cor_matrix)]

ggm_mat_p <- ggm_mat_p[rownames(ggm_mat_p) %in% rownames(cor_matrix),]

ggm_mat_p <- ggm_mat_p[,colnames(ggm_mat_p) %in% colnames(cor_matrix)]

# Set the threshold 

par_cor_p_thre <- 0.01/(151*150/2)

ggm_mat[ggm_mat_p >= par_cor_p_thre] <- 0

diag(ggm_mat) <- 0

# Create named list with bn covariates
# name of element is reponse variable
# variabels are covariates

list_resp_cov <- apply(ggm_mat,1,function(x) 
	ifelse(sum(x!=0)>0, return(names(x)[x!=0]), return(NA))
	)

# Load lambdas GC

lambda <- read.table('../data/GGM_cGWAS_gc_lambda.txt', 
	stringsAsFactors = FALSE,
	head = TRUE)

trait_names <- lambda$trait

lambda <- lambda$gc_lambda

names(lambda) <- trait_names

rm(trait_names)

# Check whether folder for results of cGAS on BN exists

if(!dir.exists('../results/GGM'))
	dir.create('../results/GGM/')

for(trait in names(list_resp_cov)){

	response <- trait

	covariates <- list_resp_cov[[response]]

	if(is.na(covariates[1]))
		covariates <- NULL
	
	res <- exact_cGWAS(
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
		all_CR = all_CR,
		path_uGWAS = '../data/uGWAS_snps_from_paper'

	)

	res$results$Pval_GC <- pchisq(res$results$chi2/lambda[trait], df = 1, lower.tail = FALSE)

	write.table(data.frame(SNP = res$snps, res$results), 
		quote = FALSE,
		row.names = FALSE,
		sep = '\t',
		file = paste0('../results/GGM/',trait,'.txt')
		)

}