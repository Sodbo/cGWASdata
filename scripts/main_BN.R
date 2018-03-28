
# Loading functions
source('../code/exact_cGWAS.R')

# Loading correlation matrix for 151 metabolites

cor_matrix <- read.table('../data/20171207_corr_matrix.txt')

# Loading information about SNPs, where column 'var_1785' 
# contains variance of the SNP

snp_info <- data.table::fread('../data/30_SNP_information.txt',
	data.table = FALSE)

all_varg <- snp_info$varg_1785

names(all_varg) <- snp_info$SNP

# Load matrix with biochemical distances

bn_mat <- read.table('../data/20171207_biochemical_distances.txt', 
	head = TRUE, 
	stringsAsFactors = FALSE
	)

bn_mat <- bn_mat[rownames(bn_mat) %in% rownames(cor_matrix),]

bn_mat <- bn_mat[,colnames(bn_mat) %in% colnames(cor_matrix)]

# Create named list with bn covariates
# name of element is reponse variable
# variabels are covariates

list_resp_cov <- apply(bn_mat,1,function(x) 
	ifelse(sum(x==1)>0, return(names(x)[x==1]), return(NA))
	)

for(trait in names(list_resp_cov)){

	response <- trait

	covariates <- list_resp_cov[[response]]

	if(is.na(covariates))
		covariates <- NULL

	print(trait)

	print(covariates)
	
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
		all_CR = rep(1, length(all_varg)),
		path_uGWAS = '../data/uGWAS_snps_from_paper'

	)

	write.table(data.frame(SNP = res$snps, res$results), 
		quote = FALSE,
		row.names = FALSE,
		sep = '\t',
		file = paste0('../results/',trait,'.tsv')
		)

}