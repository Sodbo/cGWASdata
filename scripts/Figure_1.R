# Sodbo Sharapov, Yakov Tsepilov (c)

# Create directory for storing a figure

if(!dir.exists('../results'))
  dir.create('../results')

# pdf('../results/figure_1.pdf')

# Load table with list of SNPs and metabolites for which we want to create figure 1

assoc_2_plot <- read.table("../data/BN_snps_traits_covariates.txt",
                      header = TRUE,
                      stringsAsFactors = FALSE,
                      sep="\t",
                      fill = T)

# Load SNP infor with var(g)

snp_info <- read.table('../data/30_SNP_information.txt',
                       head = TRUE,
                       stringsAsFactors = FALSE)

var_snp <- snp_info$varg_1785

names(var_snp) <- snp_info$SNP

rm(snp_info)

xv=NULL

yv1=yv2=yv3=yv4=yv5=NULL

for (i in 1:nrow(assoc_2_plot)){
  
  # Assign trait and SNP into variabels
  snp <- assoc_2_plot$SNP[i]
  
  trait <- assoc_2_plot$Metabolite[i]
  
  # Load file with uGAS statistics
  ugas_file <- paste0('../data/uGWAS_snps_from_paper/',trait,'.txt')
  
  ugas_stats <- read.table(ugas_file,
                           head = TRUE,
                           stringsAsFactors = FALSE)
  
  # Get SE and Chi2 from uGAS stats
  
  beta_snp_u <- ugas_stats$b[ugas_stats$SNP == snp]
  
  se_snp_u <- ugas_stats$se[ugas_stats$SNP == snp]
  
  chi2_snp_u <- (ugas_stats$Z[ugas_stats$SNP == snp])^2
  
  rm(ugas_file, ugas_stats)
  
  # Load file with biochemical network cGAS statistics
  bngas_file <- paste0('../results/BN/',trait,'.txt')
  
  bngas_stats <- read.table(bngas_file,
                           head = TRUE,
                           stringsAsFactors = FALSE)
  
  # Get SE and Chi2 from biochemical network cGAS statistics
  
  se_snp_c <- bngas_stats$se[bngas_stats$SNP == snp]
  
  chi2_snp_c <- bngas_stats$chi2[bngas_stats$SNP == snp]
  
  
  # Get betas for SNP and covariates
  
  beta_snp_c <- bngas_stats$b[bngas_stats$SNP == snp]
  
  if(is.numeric(bngas_stats$b_cov[bngas_stats$SNP == snp])){
    
    beta_cov_c <- bngas_stats$b_cov[bngas_stats$SNP == snp]
  
    } else{
      
      beta_cov_c <- unlist(strsplit(bngas_stats$b_cov[bngas_stats$SNP == snp],";"))
    
      beta_cov_c <- as.numeric(beta_cov_c)
      
  }
  
  # Assign list of covariates
  
  covariates <- unlist(strsplit(bngas_stats$Covariates[i],";"))
  
  rm(bngas_file, bngas_stats)
  
  beta_cov_u <- NULL
  
  for(Covariate in covariates) {
    
    cov_file <- paste0('../data/uGWAS_snps_from_paper/',Covariate,'.txt')
    
    cov_stats <- read.table(cov_file,
                            head = TRUE,
                            stringsAsFactors = FALSE)
    
    beta_cov_u <- c(beta_cov_u, cov_stats$beta[cov_stats$SNP == snp])
    
  }
  
  cor_trait_snp  <- beta_snp_u * sqrt(var_snp[snp])
  
  cor_cov_snp <- beta_cov_u * sqrt(var_snp[snp])
  
  rm(cov_file, Covariate, cov_stats)
  
  # Load table with correlation matrix for phenotypes
  
  cor_matrix <- read.table('../data/20171207_corr_matrix.txt')
  
  cor_trait_cov <- as.numeric(cor_matrix[trait, covariates])
  
  rm(cor_matrix)
  
  f1 <- log10( (se_snp_u / se_snp_c )^2 )
  f2 <- log10( ( 1 - ( beta_cov_c %*% cor_cov_snp ) / cor_trait_snp )^2 )
  f2_ <- log10( ( 1 - ( cor_trait_cov %*% cor_cov_snp ) / cor_trait_snp )^2 )
  
  xv <- c(xv, log10 (chi2_snp_c / chi2_snp_u) )
  
  yv1 <- c(yv1, ( f1 + f2 ) )
  yv2 <- c(yv2, ( f1 + f2_) )
  yv3 <- c(yv3, ( f1) )
  yv4 <- c(yv4, ( f2) )
  yv5 <- c(yv5, ( f2_) )
  
  
}

out=cbind(xv,yv1,yv2,yv3,yv4,yv5,gene=assoc_2_plot$Gene,col=colors()[200])
out[out[,"gene"]=="ACADM","col"]=colors()[575]
out[out[,"gene"]=="SLC22A4","col"]=colors()[119]
out[out[,"gene"]=="PLEKHH1","col"]=colors()[132]

plot(xv,yv1,col=out[,"col"],xlab='log10(Chi2c/Chi2u)',ylab="Components",main="Decomposition plot",pch="*",
     xlim=c(min(xv),max(xv)),ylim=c(min(yv1,yv2,yv3,yv4,yv5,na.rm=T),max(yv1,yv2,yv3,yv4,yv5,na.rm=T)),cex=4.5,bg=colors()[233])

points(xv,yv3,col=out[,"col"],pch=17,cex=2)
points(xv,yv4,col=out[,"col"],pch=15,cex=2)
abline(v=0,col="grey92")
abline(h=0,col="grey92")

rr=seq(from=min(xv),to=max(xv),length.out =100)

a=lm(yv1~xv)$coef[1];b=lm(yv1~xv)$coef[2]
points(rr,a+b*rr,col="darkgrey",type="l")

rr=seq(min(xv),max(xv),by=0.2)
a=lm(yv3~xv)$coef[1];b=lm(yv3~xv)$coef[2]
points(rr,a+b*rr,col="black",lty=3,type="l") #,type="c"

a=lm(yv4~xv)$coef[1];b=lm(yv4~xv)$coef[2]
points(rr,a+b*rr,col="black",lty=3,type="l") # ,type="c"

#abline(a=lm(yv1~xv)[1],b=lm(yv1~xv)[2],col="grey")
#abline(a=lm(yv3~xv)[1],b=lm(yv3~xv)[2],col="grey")
#abline(a=lm(yv4~xv)[1],b=lm(yv4~xv)[2],col="grey")

#points(xv,yv5,col="blue")
#legend(x=1,y=-0.3,legend=c("ACADM","SLC22A4","PLEKHH1","all others"),
#       lty=c(1,1),lwd=c(2.5,2.5),col=c("green","brown","turquoise","pink"),
#       cex=1)

ll=c("ACADM","SLC22A4","PLEKHH1")
ll=sapply(ll,FUN=function(x){eval(parse(text=paste("expression(italic(",x,"))",sep="")))})
ll=c(ll,"other loci")
legend(x=1,y=-0.3,legend=ll,
       pch=16,col=c(out[1,"col"],out[3,"col"],out[9,"col"],out[14,"col"]),
       cex=1)




#text(x=xv,y=yv1,labels=out[,"gene"],cex=0.5)
#cvz=c(2,4,5,9)
cvz=c(2,13,3,14)

#points(x=out[cvz,"xv"],y=rep(2.6,length(cvz)),col="darkgreen",pch="*",cex=4)

ll=sapply(out[cvz,"gene"],FUN=function(x){eval(parse(text=paste("expression(italic(",x,"))",sep="")))})
text(x=xv[cvz],y=apply(cbind(yv3[cvz],yv1[cvz]),1,max),labels=ll,cex=1,offset = 0.5,pos=3)  

i=1
for (i in 1:length(cvz)){
  
  ymin=min(yv1[cvz[i]],yv3[cvz[i]],yv4[cvz[i]],0)
  ymax=max(yv1[cvz[i]],yv3[cvz[i]],yv4[cvz[i]])
  points(x=rep(xv[cvz[i]],2),y=c(ymin,ymax),type="l",col="darkgreen",lwd=2)
}


cvz=c(10,12)
#points(x=out[cvz,"xv"],y=rep(2.6,length(cvz)),col="red",pch="*",cex=4)


ll=sapply(out[cvz,"gene"],FUN=function(x){eval(parse(text=paste("expression(italic(",x,"))",sep="")))})
text(x=xv[cvz],y=apply(cbind(yv3[cvz],yv1[cvz]),1,max),labels=ll,cex=1,offset = 0.5,pos=3)  
     
i=1
for (i in 1:length(cvz)){
  
  ymin=min(yv1[cvz[i]],yv3[cvz[i]],yv4[cvz[i]],0)
  ymax=max(yv1[cvz[i]],yv3[cvz[i]],yv4[cvz[i]],0)
  points(x=rep(xv[cvz[i]],2),y=c(ymin,ymax),type="l",col="red",lwd=2)
}

# dev.off()

