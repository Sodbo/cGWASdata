# Sodbo Sharapov

# Load table 
tab1 <- xlsx::read.xlsx(
	file='/home/sodbo/Dropbox/Yakov/20180326/Supplementary Table 1.SSh.xlsx',
	sheetIndex = 2, 
	startRow = 3, 
	endRow = 25
	)

for(i in 1:nrow(tab1)){

	metab_file_bn <- paste0('../results/BN/',tab1$trait[i],'.txt')

	metab_file_ggm <- paste0('../results/GGM/',tab1$trait[i],'.txt')

	metab_bn <- read.table(
		file = metab_file_bn,
		head = TRUE,
		stringsAsFactors = FALSE)

	metab_ggm <- read.table(
		file = metab_file_ggm,
		head = TRUE,
		stringsAsFactors = FALSE)

	tab1$bdGAS_b[i] <- metab_bn$b[metab_bn$SNP == tab1$SNP[i]]

	tab1$bdGAS_se[i] <- metab_bn$se[metab_bn$SNP == tab1$SNP[i]]

	tab1$bdGAS_p[i] <- metab_bn$Pval_GC[metab_bn$SNP == tab1$SNP[i]]

	tab1$cGAS_b[i] <- metab_ggm$b[metab_ggm$SNP == tab1$SNP[i]]

	tab1$cGAS_se[i] <- metab_ggm$se[metab_ggm$SNP == tab1$SNP[i]]

	tab1$cGAS_p[i] <- metab_ggm$Pval_GC[metab_ggm$SNP == tab1$SNP[i]]

}

write.table(tab1[,c('bdGAS_b','bdGAS_se','bdGAS_p')], quote=FALSE,row.names=FALSE)

write.table(tab1[,c('cGAS_b','cGAS_se','cGAS_p')], quote=FALSE,row.names=FALSE)