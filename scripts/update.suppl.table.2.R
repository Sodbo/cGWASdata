# Sodbo Sharapov

# Load table 
tab2 <- xlsx::read.xlsx(
	file='~/Dropbox/Yakov/20180326/Supplementary Table 2.SSh.xlsx',
	sheetIndex=2, 
	startRow=3, 
	endRow=33
	)

for(i in 1:nrow(tab2)){

	metab_file <- paste0('../results/GGM/',tab2$Metabolite[i],'.txt')

	metab_cgas <- read.table(
		file = metab_file,
		head = TRUE,
		stringsAsFactors = FALSE)

	tab2$c_beta[i] <- metab_cgas$b[metab_cgas$SNP == tab2$SNP[i]]

	tab2$c_se[i] <- metab_cgas$se[metab_cgas$SNP == tab2$SNP[i]]

	tab2$c_z[i] <- tab2$c_beta[i] / tab2$c_se[i]

	tab2$c_p[i] <- metab_cgas$Pval_GC[metab_cgas$SNP == tab2$SNP[i]]

}

print(tab2[,c('SNP','Metabolite','c_beta','c_se','c_z','c_p')])

print('\n')

write.table(tab2[,c('c_beta','c_se','c_z','c_p')], 
	quote=FALSE,
	row.names=FALSE)