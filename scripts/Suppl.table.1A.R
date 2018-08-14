# Scritp for reproducing Supplementary Table 1A

Pval_thr <- 5e-8/105

# Load table with biochemical distances between 153 metabolites
# and restrict it to a subset of 105 metabolites for which at least the
# one-reaction-step immediate biochemical neighbors are known

bn_dist <- read.table(
	'data/20171207_biochemical_distances.txt',
	header = TRUE,
	row.names = 1
	)

bn_dist[bn_dist != 1] <- 0

# Get list of 105 metabolites

traits <- names(which(apply(bn_dist,2,sum)>0))

traits_exists <- sub('.txt','',list.files('data/uGWAS'))

traits <- intersect(traits,traits_exists)

rm(traits_exists)

snp_info <- data.table::fread('zcat data/SNP_information.txt.zip', data.table = FALSE)

snp_info$maf <- pmin(1 - snp_info$freq, snp_info$freq)

minfreq=0.1
minP22=0
minR2=0.3
hw=1e-6
CR=0.95
#
gut_snps <- snp_info$SNP[which(snp_info[,"maf"]*as.numeric(snp_info[,"R2_impute_info"])>=minfreq &
                                      as.numeric(snp_info[,"R2_impute_info"])>=minR2 &
                                      snp_info[,"hw"]>=hw)]

source('scripts/clumping_functions.R')

locus_table <- function_for_making_full_table_without_gcv('data/uGWAS/', 
                                                          traits, 
                                                          snp_info, 
                                                          thr=5e-8, 
                                                          delta = 5e5)

locus_table_2 <- function_for_shlop_24_10_2013(
  locus_table,
  p_value="P-value",
  pos="Position",
  snp="SNP",
  delta=5e5,
  chr="Chromosome")

locus_table_2[locus_table_2$`P-value`< 5e-8/151,]

	"metab_file_bn <- paste0('results/BN/',trait,'.txt')
	metab_bn <- read.table(
		file = metab_file_bn,
		head = TRUE,
		stringsAsFactors = FALSE
		)

	metab_bn <- metab_bn[metab_bn$Pval < Pval_thr,]

	locus_tab_bn <- rbind(locus_tab_bn,metab_bn) "

	"tab1$uGAS_b[i] <- metab_u$beta[metab_u$SNP == snp]

	tab1$uGAS_se[i] <- metab_u$se[metab_u$SNP == snp]

	tab1$uGAS_p[i] <- pchisq((metab_u$Z[metab_u$SNP == snp])^2/lambda_ugas$gc_lambda[lambda_ugas$trait == trait],
		df = 1,
		lower.tail = FALSE
		)

	tab1$bdGAS_b[i] <- metab_bn$b[metab_bn$SNP == snp]

	tab1$bdGAS_se[i] <- metab_bn$se[metab_bn$SNP == snp]

	tab1$bdGAS_p[i] <- metab_bn$chi2[metab_bn$SNP == snp] / lambda_bn[trait]
	
	tab1$bdGAS_p[i] <- pchisq(tab1$bdGAS_p[i], df=1 , lower.tail = FALSE)

	tab1$bd_noise_comp <- log10((tab1$uGAS_se / tab1$bdGAS_se)^2)
	
	tab1$cGAS_b[i] <- metab_ggm$b[metab_ggm$SNP == snp]

	tab1$cGAS_se[i] <- metab_ggm$se[metab_ggm$SNP == snp]

	tab1$cGAS_p[i] <- metab_ggm$chi2[metab_ggm$SNP == snp] / lambda_ggm[trait]

	tab1$cGAS_p[i] <- pchisq(tab1$cGAS_p[i], df=1 , lower.tail = FALSE)"

"write.table(tab1, 
	quote = FALSE,
	row.names = FALSE,
	sep = '\t',
#	dec = ',',
	file = 'results/Suppl_table_1A.csv'
)"


"# Load lambdas GC for uGWAS
lambda_ugas <- read.table(
	'data/uGWAS_gc_lambda.txt', 
	head = TRUE, 
	stringsAsFactors = FALSE,
	row.names = 1
	)


# Load lambdas GC for GGM GWAS

lambda_ggm <- read.table(
	'data/GGM_cGWAS_gc_lambda.txt', 
	stringsAsFactors = FALSE,
	head = TRUE,
	row.names = 1
	)

lambda_ggm <- lambda_ggm[traits,]

names(lambda_ggm) <- traits

# Load lambdas GC for BN GWAS

lambda_bn <- read.table(
	'data/BN_cGWAS_gc_lambda.txt', 
	stringsAsFactors = FALSE,
	head = TRUE,
	row.names = 1
	)

lambda_bn <- lambda_bn[traits,]

names(lambda_bn) <- traits"