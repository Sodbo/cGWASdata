# Sodbo Sharapov, Yakov Tsepilov (c)

if(!dir.exists('../results'))
  dir.create('../results')

# pdf('../results/figure_1.pdf')

load("../data/data_explained_var_new.RData")

xb_snps <- read.table("../data/BN_snps_traits_covariates.txt",header = T,stringsAsFactors = F,sep="\t",fill = T)

xv=NULL

yv1=yv2=yv3=yv4=yv5=NULL

for (i in 1:nrow(xb_snps)){
  
  # Assign trait and SNP into variabels
  snp <- xb_snps$SNP[i]
  
  trait <- xb_snps$Metabolite[i]
  
  # Assign list of covariates
  cvrts <- unlist(strsplit(xb_snps$Covariates[i],";"))
  
  dta[is.na(dta[,snp]),snp]=mean(dta[,snp],na.rm=T)
  sdta=cbind(snp=dta[,snp],phedata[,c(cvrts,trait)])

  flac=paste(trait,"~snp+",paste(cvrts,collapse="+"))

  func=lm(flac,data=as.data.frame(sdta))
  
  # Load file with uGAS statistics
  ugas_file <- paste0('../data/uGWAS_snps_from_paper/',trait,'.txt')
  
  ugas_stats <- read.table(ugas_file,
                           head = TRUE,
                           stringsAsFactors = FALSE)
  
  # Get SE and Chi2 from uGAS stats
  se_u <- ugas_stats$se[ugas_stats$SNP == snp]
  
  chi2_u <- (ugas_stats$Z[ugas_stats$SNP == snp])^2
  
  rm(ugas_file, ugas_stats)
  
  # Load file with biochemical network cGAS statistics
  bngas_file <- paste0('../results/BN/',trait,'.txt')
  
  bngas_stats <- read.table(bngas_file,
                           head = TRUE,
                           stringsAsFactors = FALSE)
  
  # Get SE and Chi2 from biochemical network cGAS statistics
  se_c <- bngas_stats$se[bngas_stats$SNP == snp]
  
  chi2_c <- bngas_stats$chi2[bngas_stats$SNP == snp]
  
  rm(bngas_file, bngas_stats)

  bcg=summary(func)$coef["snp","Estimate"]
  bcc=summary(func)$coef[cvrts,"Estimate"]
  
  ryg=cor(sdta[,trait],sdta[,"snp"])=bug*sqrt(varg)
  rcg=cor(sdta[,cvrts],sdta[,"snp"])
  ryc=cor(sdta[,trait],sdta[,cvrts])
  
  f1 <- log10( (se_u / se_c ) ^ 2 )
  f2=log10((1-(bcc%*%rcg)/ryg)^2)
  f2_=log10((1-(ryc%*%rcg)/ryg)^2)
  
  xv <- c(xv,log10(chi2_c / chi2_u))
  
  yv1=c(yv1,(f1+f2))
  yv2=c(yv2,(f1+f2_))
  yv3=c(yv3,(f1))
  yv4=c(yv4,(f2))
  yv5=c(yv5,(f2_))
  
  
}


out=cbind(xv,yv1,yv2,yv3,yv4,yv5,gene=xb_snps$Gene,col=colors()[200])
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

