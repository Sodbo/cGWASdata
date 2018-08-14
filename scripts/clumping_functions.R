function_for_making_full_table_without_gcv <- function(path_sumstats, snp_col, beta_col, se_col, Pval_col, traits, snp_info, thr = 5e-8, delta = 2.5e5){
		pheindex = 1
		locus_table = as.data.frame(NULL)
		j = 0
		for (trait in traits){
		  #print(trait)
			file_name = paste0(path_sumstats,trait,".txt")
			Zx <- data.table::fread(file_name, data.table=FALSE)
			Zx <- Zx[,c(snp_col, beta_col, se_col, Pval_col)]
			colnames(Zx)=c("SNP","beta_SNP","se_SNP","P-value")
			
			Zx=Zx[Zx[,"P-value"]<=thr,]

			kostil=FALSE
			if (is.vector(Zx)){
					kostil=TRUE
				} else {
					if (nrow(Zx)>0){
						kostil=TRUE
					}
				}

			if (kostil) {
				snps=Zx[,"SNP"]
				snps=unique(snps)
				index=match(snps,Zx[,"SNP"])
				Zx=Zx[index,]
				
				index=match(snps,snp_info$SNP)
				Chromosome=snp_info[index,"chr"]
				Position=snp_info[index,"pos"]
				Zx=cbind(Zx,Chromosome,Position)
				Zx=cbind(Zx,trait)
							
				snps=Zx[,"SNP"]
				rownames(Zx)=snps
				snps=intersect(snps,gut_snps)
				
				if (length(snps)>1){	
					Zx=Zx[snps,]
					Zx=Zx[order(Zx$Chromosome,Zx$Position),]
					i=2
					while (i<=(dim(Zx)[1])){
						if ((abs(Zx$Position[i]-Zx$Position[i-1])<=delta)&(Zx$Chromosome[i]==Zx$Chromosome[i-1])){
							if (Zx[i,"P-value"]<=Zx[i-1,"P-value"]){
								Zx=Zx[-(i-1),]
							} else{
								Zx=Zx[-i,]
							}
						} else{
							i=i+1
						}
					}
					locus_table=rbind(locus_table,Zx)
				}
				if (length(snps)==1){
					Zx=Zx[snps,]
					locus_table=rbind(locus_table,Zx)
				}
			}
			
			#print(locus_table)
		}
		locus_table[,"P-value"]=as.numeric(locus_table[,"P-value"])
		locus_table[,"Position"]=as.numeric(locus_table[,"Position"])
		locus_table=locus_table[order(locus_table$Chromosome,locus_table$Position),]
		return(locus_table)
	}
	
	function_for_shlop_24_10_2013=function(locus_table,p_value="imp_P_value",pos="imp_pos",snp="imp_SNP",delta=5e5,chr="chr"){
		locus_table[,p_value]=as.numeric(locus_table[,p_value])
		locus_table[,pos]=as.numeric(locus_table[,pos])
		Zx <-locus_table
		Zx=Zx[order(Zx[,chr],Zx[,pos]),]
		n_traits=1
		Zx=cbind(Zx,n_traits)
		i=2
		while (i<=(dim(Zx)[1])){
			if ((abs(Zx[i,pos]-Zx[i-1,pos])<=delta)&(Zx[i,chr]==Zx[i-1,chr])){
				if (Zx[i,p_value]<=Zx[i-1,p_value]){
					Zx=Zx[-(i-1),]
					Zx[i-1,"n_traits"]=Zx[i-1,"n_traits"]+1
				} else{
					Zx=Zx[-i,]
					Zx[i-1,"n_traits"]=Zx[i-1,"n_traits"]+1
				}
			} else{
				i=i+1
			}
		}
		locus_table<-Zx
		rownames(locus_table)=as.character(locus_table[,snp])
		return(locus_table)
	}