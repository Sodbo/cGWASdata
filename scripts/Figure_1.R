#setwd("/Users/tsepilov/Yandex.Disk.localized/Work/Projects/multivarate project_2013_2015/data/")
#load("/Users/tsepilov/Yandex.Disk.localized/Work/Projects/multivarate project_2013_2015/data/data_explained_var.RData")

setwd("/mnt/Disk_D/YD_ux/Work/Projects/multivarate project_2013_2015/data/")
#load("//mnt/Disk_D/YD_ux/Work/Projects/multivarate project_2013_2015/data/data_explained_var.RData")
load("/mnt/Disk_D/YD_ux/Work/Projects/multivarate project_2013_2015/data/data_explained_var_new.RData")

xb=read.table("20161111_table_cGWAS_uGWAS_bdGWAS.txt",header = T,stringsAsFactors = F,sep="\t",fill = T)
xb_snps=read.table("20171205_table1.txt",header = T,stringsAsFactors = F,sep="\t",fill = T)

xb_1=paste(xb[,1],xb[,2],sep="_")
xb_2=paste(xb_snps[,3],xb_snps[,4],sep="_")

xb=xb[xb_1%in%xb_2,]
xb_1=paste(xb[,1],xb[,2],sep="_")
rownames(xb)=xb_1
xb=xb[xb_2,]

nrow(xb)
length(table(xb[,"Locus"]))


i=1

xv=NULL
yv1=yv2=yv3=yv4=yv5=NULL

for (i in (1:dim(xb)[1])){
  snp=xb[i,"SNP"]
  trait=xb[i,"trait"]
  cvrts=unlist(strsplit(xb[i,"bd_cvrts"],";"))
  dta[is.na(dta[,snp]),snp]=mean(dta[,snp],na.rm=T)
  sdta=cbind(snp=dta[,snp],phedata[,c(cvrts,trait)])
  
  N=dim(dta)[1]
  
  
  flau=paste(trait,"~snp")
  if (length(cvrts)>0){
    flac=paste(trait,"~snp+",paste(cvrts,collapse="+"))
  } else flac=flau
  
  funu=lm(flau,data=as.data.frame(sdta))
  func=lm(flac,data=as.data.frame(sdta))
  
  x2u=summary(funu)$coef["snp","t value"]^2
  x2c=summary(func)$coef["snp","t value"]^2
  
  seu=summary(funu)$coef["snp","Std. Error"]
  sec=summary(func)$coef["snp","Std. Error"]
  
  bcg=summary(func)$coef["snp","Estimate"]
  bcc=summary(func)$coef[cvrts,"Estimate"]
  
  ryg=cor(sdta[,trait],sdta[,"snp"])
  rcg=cor(sdta[,cvrts],sdta[,"snp"])
  ryc=cor(sdta[,trait],sdta[,cvrts])
  
  f1=log10((seu/sec)^2)
  f2=log10((1-(bcc%*%rcg)/ryg)^2)
  f2_=log10((1-(ryc%*%rcg)/ryg)^2)
  
  xv=c(xv,log10(x2c/x2u))
  yv1=c(yv1,(f1+f2))
  yv2=c(yv2,(f1+f2_))
  yv3=c(yv3,(f1))
  yv4=c(yv4,(f2))
  yv5=c(yv5,(f2_))
  
  
}


out=cbind(xv,yv1,yv2,yv3,yv4,yv5,gene=xb[,"Gene"],col=colors()[200])
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
t#ext(x=xv[cvz],y=apply(cbind(yv3[cvz],yv1[cvz]),1,max),labels=ll,cex=1,offset = 0.5,pos=3)  

i=1
for (i in 1:length(cvz)){
  
  ymin=min(yv1[cvz[i]],yv3[cvz[i]],yv4[cvz[i]],0)
  ymax=max(yv1[cvz[i]],yv3[cvz[i]],yv4[cvz[i]])
  points(x=rep(xv[cvz[i]],2),y=c(ymin,ymax),type="l",col="darkgreen",lwd=2)
}


cvz=c(10,12)
#points(x=out[cvz,"xv"],y=rep(2.6,length(cvz)),col="red",pch="*",cex=4)


ll=sapply(out[cvz,"gene"],FUN=function(x){eval(parse(text=paste("expression(italic(",x,"))",sep="")))})
#text(x=xv[cvz],y=apply(cbind(yv3[cvz],yv1[cvz]),1,max),labels=ll,cex=1,offset = 0.5,pos=3)  
     
i=1
for (i in 1:length(cvz)){
  
  ymin=min(yv1[cvz[i]],yv3[cvz[i]],yv4[cvz[i]],0)
  ymax=max(yv1[cvz[i]],yv3[cvz[i]],yv4[cvz[i]],0)
  points(x=rep(xv[cvz[i]],2),y=c(ymin,ymax),type="l",col="red",lwd=2)
}