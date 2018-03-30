# Sodbo Sharapov

# Load table 

tab1 <- read.table(
	file='data/SNPs_for_suppl_tab1.txt',
	head = TRUE,
	stringsAsFactors = FALSE
	)

# Load lambdas GC for uGWAS
lambda_ugas <- read.table('data/uGWAS_gc_lambda.txt', 
	head = TRUE, 
	stringsAsFactors = FALSE)


# Load lambdas GC for GGM GWAS

lambda_ggm <- read.table('data/GGM_cGWAS_gc_lambda.txt', 
	stringsAsFactors = FALSE,
	head = TRUE)

trait_names <- lambda_ggm$trait

lambda_ggm <- lambda_ggm$gc_lambda

names(lambda_ggm) <- trait_names

# Load lambdas GC for BN GWAS

lambda_bn <- read.table('data/BN_cGWAS_gc_lambda.txt', 
	stringsAsFactors = FALSE,
	head = TRUE)

trait_names <- lambda_bn$trait

lambda_bn <- lambda_bn$gc_lambda

names(lambda_bn) <- trait_names

rm(trait_names)

for(i in 1:nrow(tab1)){

	trait <- tab1$trait[i]

	snp <- tab1$SNP[i]

	metab_file_u <- paste0('data/uGWAS_snps_from_paper/',trait,'.txt')

	metab_file_bn <- paste0('results/BN/',trait,'.txt')

	metab_file_ggm <- paste0('results/GGM/',trait,'.txt')

	metab_u <- read.table(
		file = metab_file_u,
		head = TRUE,
		stringsAsFactors = FALSE)

	metab_bn <- read.table(
		file = metab_file_bn,
		head = TRUE,
		stringsAsFactors = FALSE)

	metab_ggm <- read.table(
		file = metab_file_ggm,
		head = TRUE,
		stringsAsFactors = FALSE)

	tab1$uGAS_b[i] <- metab_u$beta[metab_u$SNP == snp]

	tab1$uGAS_se[i] <- metab_u$se[metab_u$SNP == snp]

	tab1$uGAS_p[i] <- pchisq((metab_u$Z[metab_u$SNP == snp])^2/lambda_ugas$gc_lambda[lambda_ugas$trait == trait],
		df = 1,
		lower.tail = FALSE
		)


	tab1$bdGAS_b[i] <- metab_bn$b[metab_bn$SNP == snp]

	tab1$bdGAS_se[i] <- metab_bn$se[metab_bn$SNP == snp]

	tab1$bdGAS_p[i] <- metab_bn$Pval[metab_bn$SNP == snp]

	tab1$bdGAS_p[i] <- metab_bn$chi2[i] / lambda_bn

	tab1$bdGAS_p[i] <- pchisq(tab1$bdGAS_p[i], df =1 , lower.tail = FALSE)


	tab1$cGAS_b[i] <- metab_ggm$b[metab_ggm$SNP == snp]

	tab1$cGAS_se[i] <- metab_ggm$se[metab_ggm$SNP == snp]

	tab1$cGAS_p[i] <- metab_ggm$chi2[i] / lambda_ggm

	tab1$cGAS_p[i] <- pchisq(tab1$cGAS_p[i], df =1 , lower.tail = FALSE)


}

write.table(tab1, 
	quote = FALSE,
	row.names = FALSE,
	sep='\t',
	dec=','
)