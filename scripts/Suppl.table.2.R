# Sodbo Sharapov

# Load table 

tab2 <- read.table(
	file='data/SNPs_for_suppl_tab2.txt',
  stringsAsFactors = FALSE,
	head = TRUE
	)

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

for(i in 1:nrow(tab2)){

	trait <- tab2$Metabolite[i]

	snp <- tab2$SNP[i]

	metab_file_u <- paste0('data/uGWAS_snps_from_paper/',trait,'.txt')

	metab_file_ggm <- paste0('results/GGM/',trait,'.txt')

	metab_u <- read.table(
		file = metab_file_u,
		head = TRUE,
		stringsAsFactors = FALSE)

	metab_cgas <- read.table(
		file = metab_file_ggm,
		head = TRUE,
		stringsAsFactors = FALSE)

	tab2$u_beta[i] <- metab_u$beta[metab_u$SNP == snp]

	tab2$u_se[i] <- metab_u$se[metab_u$SNP == snp]

	tab2$u_Z[i] <- tab2$u_beta[i] / tab2$u_se[i]

	tab2$u_p[i] <- pchisq((tab2$u_Z[i])^2/lambda_ugas$gc_lambda[lambda_ugas$trait == trait],
		df = 1,
		lower.tail = FALSE
		)


	tab2$c_beta[i] <- metab_cgas$b[metab_cgas$SNP == snp]

	tab2$c_se[i] <- metab_cgas$se[metab_cgas$SNP == snp]

	tab2$c_z[i] <- tab2$c_beta[i] / tab2$c_se[i]

	tab2$c_p[i] <- (tab2$c_z[i])^2 / lambda_ggm[trait]
	
	tab2$c_p[i] <- pchisq(tab2$c_p[i], df = 1, lower.tail = FALSE)

}

write.table(tab2, 
	quote = FALSE,
	row.names = FALSE,
#	dec = ',',
	sep = '\t',
	file = 'results/Suppl_table_2.csv'
	)