# Sodbo Sharapov

# Load table 

tab1 <- xlsx::read.xlsx(
	file='~/Dropbox/Yakov/20180326/Supplementary Table 1.SSh.xlsx',
	sheetIndex = 2, 
	startRow = 3, 
	endRow = 25
	)

lambda_ugas <- read.table('../data/uGWAS_gc_lambda.txt', 
	head = TRUE, 
	stringsAsFactors = FALSE)

for(i in 1:nrow(tab1)){

	trait <- tab1$trait[i]

	snp <- tab1$SNP[i]

	metab_file_u <- paste0('../data/uGWAS_snps_from_paper/',trait,'.txt')

	metab_file_bn <- paste0('../results/BN/',trait,'.txt')

	metab_file_ggm <- paste0('../results/GGM/',trait,'.txt')

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

	tab1$bdGAS_p[i] <- metab_bn$Pval_GC[metab_bn$SNP == snp]


	tab1$cGAS_b[i] <- metab_ggm$b[metab_ggm$SNP == snp]

	tab1$cGAS_se[i] <- metab_ggm$se[metab_ggm$SNP == snp]

	tab1$cGAS_p[i] <- metab_ggm$Pval_GC[metab_ggm$SNP == snp]

}

write.table(tab1[,c('uGAS_b','uGAS_se','uGAS_p')], 
	quote = FALSE,
	row.names = FALSE,
	sep='\t',
	dec=','
	)

write.table(tab1[,c('bdGAS_b','bdGAS_se','bdGAS_p')], 
	quote = FALSE,
	row.names = FALSE,
	sep='\t',
	dec=','
	)

write.table(tab1[,c('cGAS_b','cGAS_se','cGAS_p')],
	quote = FALSE,
	row.names = FALSE,
	sep='\t',
	dec=','
	)