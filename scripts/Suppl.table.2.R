# Sodbo Sharapov

# Load table 

tab2 <- xlsx::read.xlsx(
	file='~/Dropbox/Yakov/20180326/Supplementary Table 2.SSh.xlsx',
	sheetIndex=2, 
	startRow=3, 
	endRow=33
	)

lambda_ugas <- read.table('../data/uGWAS_gc_lambda.txt', 
	head = TRUE, 
	stringsAsFactors = FALSE)

for(i in 1:nrow(tab2)){

	trait <- tab2$Metabolite[i]

	snp <- tab2$SNP[i]

	metab_file_u <- paste0('../data/uGWAS_snps_from_paper/',trait,'.txt')

	metab_file_ggm <- paste0('../results/GGM/',trait,'.txt')

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

	tab2$c_p[i] <- metab_cgas$Pval_GC[metab_cgas$SNP == snp]

}

write.table(tab2[,c('u_beta','u_se','u_Z','u_p')], 
	quote=FALSE,
	row.names=FALSE,
	dec=',',
	sep='\t')


write.table(tab2[,c('c_beta','c_se','c_z','c_p')], 
	quote=FALSE,
	row.names=FALSE,
	dec=',',
	sep='\t')