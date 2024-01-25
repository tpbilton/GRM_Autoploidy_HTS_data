##### File for performing the analysis of manuscript 
##### "Construction of relatedness matrices in autopolyploid populations using low depth high-throughput sequencing data"

##### Notes: #####
## Supplementary files associated with the paper need to be downloaded from the webpage of the publisher.
## For this script, the following files are required:
## - Supplementary_File3.csv
## - Supplementary_file4.tab.ra
## - Supplementary_File5.csv
## - Supplementary_file6.tab.ra

########################
#### SNP-array data ####
calcGRM_vanRaden <- function(geno, ploid=4, phat, thres=0.01){
  nInd = nrow(geno)
  snpsubset <- which(phat < 1-thres & phat > thres)
  nsnpsub <- length(snpsubset)
  phat <- phat[snpsubset]
  geno = geno[,snpsubset]
  genon0 <- geno - ploid*rep.int(phat, rep(nInd, nsnpsub))
  genon0[is.na(geno)] <- 0
  P0 <- matrix(phat,nrow=nInd,ncol=nsnpsub,byrow=T)
  P1 <- 1-P0
  P0[is.na(geno)] <- 0
  P1[is.na(geno)] <- 0
  div0 <- ploid*tcrossprod(P0,P1)
  GRM <- tcrossprod(genon0)/div0
  return(GRM)
}

geno.filt = as.matrix(read.csv("Supplementary_File3.csv", row.names = 1))
phat_snp = colMeans(geno.filt, na.rm=T)/4
grm_vr = calcGRM_vanRaden(geno.filt, phat=phat_snp)

######################
##### GBS data: ######
library(GUSrelate)

## Load the data in R
ra <- readRA("Supplementary_file4.tab.ra")

# Note: Supplementary_file4.ra.tab.gz is in RA format. 
#       This format has format [reference]/[alternate] instead of a genotype.
#       E.g., 5/3 means that there are 5 reads for the reference allele and 3 reads for the alternate allele                                                     

## write out sample file for the GUSrelate package
write.csv(data.frame(ID=ra$.__enclos_env__$private$indID, Ploidy=4), file = "SampleFile.csv", row.names=FALSE, quote=FALSE)

## Create GRM object
grm_pot = makeGRM(ra, samfile="SampleFile.csv", filter=list(MAF=NULL, MISS=NULL))

## Index high and low depth SNPs
indx_high = grm_pot$.__enclos_env__$private$snpdepth  > 50
indx_low = grm_pot$.__enclos_env__$private$snpdepth  < 50

#### GUSrelate
grm_pot$computeGRM("VR_high", method="VanRaden", ep=0, snpsubset=which(indx_high), filter=list(MAF=NULL, PVALUE=NULL, MISS=NULL))
grm_pot$computeGRM("VR_low", method="VanRaden", ep=0, snpsubset=which(indx_low), filter=list(MAF=NULL, PVALUE=NULL, MISS=NULL))

grm_vr_gbs_high = grm_pot$.__enclos_env__$private$GRM$VR_high$GRM
grm_vr_gbs_low = grm_pot$.__enclos_env__$private$GRM$VR_low$GRM

#### AGHmatrix:
library(AGHmatrix)
agh_computeGRM = function(ref, alt){
  ratio = ref/(ref + alt)
  #ratio[which((ref + alt) < 20)] = NA
  ratio[which(is.na(ratio))] = -9
  return(Gmatrix(ratio, ratio=T))
}
ref = grm_pot$.__enclos_env__$private$ref
alt = grm_pot$.__enclos_env__$private$alt
agh_high = agh_computeGRM(ref[,which(indx_high)], alt[,which(indx_high)])
agh_low = agh_computeGRM(ref[,which(indx_low )], alt[,which(indx_low)])

#### Cericola et al:
GRM_Cericola = function(ref, alt, ploid=4, phat){
  ratio = ref/(ref + alt)
  depth = ref + alt
  #depth[depth == 0] = NA
  snpdepth = rowMeans(depth, na.rm = T)
  #ratio_mean = colMeans(ratio, na.rm=T)
  rmean_mat = matrix(phat, nrow=nrow(ratio), ncol=ncol(ratio), byrow=T)
  miss_indx = which(is.na(ratio), arr.ind=T) 
  ratio[miss_indx] = phat[miss_indx[,2]]
  #P0 = matrix(phat, nrow=nrow(ratio), ncol=length(phat), byrow=T)
  #P0[depth == 0] = 0
  Z = tcrossprod(ratio-rmean_mat,ratio-rmean_mat)/(1/ploid*sum(phat*(1-phat)))
  #Z = ploid*tcrossprod(ratio-rmean_mat,ratio-rmean_mat)/(tcrossprod(P0,1-P0))
  diag(Z) = diag(Z)*(1 - (ploid-1)/(snpdepth + ploid - 1))
  return(Z)
} 
cer_high = GRM_Cericola(ref[,which(indx_high)], alt[,which(indx_high)], phat=grm_pot$.__enclos_env__$private$pfreq[indx_high])
cer_low = GRM_Cericola(ref[,which(indx_low)], alt[,which(indx_low)], phat = grm_pot$.__enclos_env__$private$pfreq[indx_low])


#######################################################

## Add sample IDs to row and columns of the GRMs for GBS data. 
colnames(grm_vr_gbs_high) = colnames(grm_vr_gbs_low) = colnames(agh_low) = colnames(agh_high) = colnames(cer_low) = colnames(cer_high) =  grm_pot$.__enclos_env__$private$GRM$VR_low$indID
rownames(grm_vr_gbs_high) = rownames(grm_vr_gbs_low) = rownames(agh_low) = rownames(agh_high) = rownames(cer_low) = rownames(cer_high) =  grm_pot$.__enclos_env__$private$GRM$VR_low$indID

## Read in the list matching IDs
ID_list = read.csv("Supplementary_File5.csv")
ID_list_ord = ID_list[match(substr(colnames(grm_vr_gbs_high),1,9), ID_list[,2]),]

## match genotype GRM with GBS GRMs
neword = match(ID_list_ord[,1], colnames(grm_vr))
grm_vr_ord = grm_vr[neword, neword]

##############################################################
## Bland-altmann plots comparing the relatedness estimates ###

methods = c("CHIP","GUS (High)","AGH (High)","Cer (High)","GUS (Low)","AGH (Low)","Cer (Low)")
library(MethComp)

## Relatedness comparisons
png("Figure5.png", width=800*300/72, height=800*300/72, res=300)
noff=sum(upper.tri(grm_vr_ord))
offdiag = function(x) x[upper.tri(x)]
df_rr <- rbind(data.frame(meth="CHIP", item=1:noff, y=offdiag(grm_vr_ord)),
               data.frame(meth="GUS (High)", item=1:noff, y=offdiag(grm_vr_gbs_high)),
               data.frame(meth="AGH (High)", item=1:noff, y=offdiag(agh_high)),
               data.frame(meth="Cer (High)", item=1:noff, y=offdiag(cer_high)),
               data.frame(meth="GUS (Low)", item=1:noff, y=offdiag(grm_vr_gbs_low)),
               data.frame(meth="AGH (Low)", item=1:noff, y=offdiag(agh_low)),
               data.frame(meth="Cer (Low)", item=1:noff, y=offdiag(cer_low)))
df_rr$meth = factor(df_rr$meth, levels=c("CHIP","GUS (High)","AGH (High)","Cer (High)","GUS (Low)","AGH (Low)","Cer (Low)"))
plot(Meth(df_rr), diff.range = 1.25)
mina = 0.01
maxa = 0.13
dista = (1-mina*2)/length(methods)
xset = 1/4*dista
xset2 = 1/5*dista
yset = 9/10*dista
yset2 = 1/10*dista
for(i in 1:(length(methods)-1)){
  for(j in (i+1):length(methods)){
    temp_cor = format(round(cor(df_rr$y[which(df_rr$meth == methods[i])], df_rr$y[which(df_rr$meth == methods[j])]),2), nsmall = 2)
    temp_ccc = format(round(epiR::epi.ccc(df_rr$y[which(df_rr$meth == methods[i])], df_rr$y[which(df_rr$meth == methods[j])])$rho.c[1],2),nsmall=2)
    text(paste0("r = ",temp_cor), x=mina + (i-1)*dista + xset, y=mina + (length(methods)-j)*dista + yset, xpd=NA)
    text(paste0("ccc = ",temp_ccc), x=maxa + (i-1)*dista - xset2, y=mina + (length(methods)-j)*dista + yset2, xpd=NA)
  }
}
dev.off()

## Self Relatedness comparions
png("Figure6.png", width=800*300/72, height=800*300/72, res=300)
ndiag=length(diag(grm_vr_ord))
df_sr <- rbind(data.frame(meth="CHIP", item=1:ndiag, y=diag(grm_vr_ord)),
               data.frame(meth="GUS (High)", item=1:ndiag, y=diag(grm_vr_gbs_high)),
               data.frame(meth="AGH (High)", item=1:ndiag, y=diag(agh_high)),
               data.frame(meth="Cer (High)", item=1:ndiag, y=diag(cer_high)),
               data.frame(meth="GUS (Low)", item=1:ndiag, y=diag(grm_vr_gbs_low)),
               data.frame(meth="AGH (Low)", item=1:ndiag, y=diag(agh_low)),
               data.frame(meth="Cer (Low)", item=1:ndiag, y=diag(cer_low)))
df_sr$meth = factor(df_sr$meth, levels=c("CHIP","GUS (High)","AGH (High)","Cer (High)","GUS (Low)","AGH (Low)","Cer (Low)"))
plot(Meth(df_sr), diff.range = 1.25)
mina = 0.01
maxa = 0.13
dista = (1-mina*2)/length(methods)
xset = 1/4*dista
xset2 = 1/5*dista
yset = 9/10*dista
yset2 = 1/10*dista
for(i in 1:(length(methods)-1)){
  for(j in (i+1):length(methods)){
    temp_cor = format(round(cor(df_sr$y[which(df_sr$meth == methods[i])], df_sr$y[which(df_sr$meth == methods[j])]),2), nsmall = 2)
    temp_ccc = format(round(epiR::epi.ccc(df_sr$y[which(df_sr$meth == methods[i])], df_sr$y[which(df_sr$meth == methods[j])])$rho.c[1],2),nsmall=2)
    text(paste0("r = ",temp_cor), x=mina + (i-1)*dista + xset, y=mina + (length(methods)-j)*dista + yset, xpd=NA)
    text(paste0("ccc = ",temp_ccc), x=maxa + (i-1)*dista - xset2, y=mina + (length(methods)-j)*dista + yset2, xpd=NA)
  }
}
dev.off()


#################################
### Supplementary Figure:

ra_noHWE <- readRA("Supplementary_file6.tab.ra")

## Create GRM object
grm_pot_noHWE = makeGRM(ra_noHWE, samfile="SampleFile.csv", filter=list(MAF=NULL, MISS=NULL))

## Index high and low depth SNPs
indx_high_noHWE = grm_pot_noHWE$.__enclos_env__$private$snpdepth  > 50
indx_low_noHWE = grm_pot_noHWE$.__enclos_env__$private$snpdepth  < 50

#### GUSrelate
grm_pot_noHWE$computeGRM("VR_high_noHWE", method="VanRaden", ep=0, snpsubset=which(indx_high_noHWE), filter=list(MAF=NULL, PVALUE=NULL, MISS=NULL))
grm_pot_noHWE$computeGRM("VR_low_noHWE", method="VanRaden", ep=0, snpsubset=which(indx_low_noHWE), filter=list(MAF=NULL, PVALUE=NULL, MISS=NULL))

grm_vr_gbs_high_noHWE = grm_pot_noHWE$.__enclos_env__$private$GRM$VR_high_noHWE$GRM
grm_vr_gbs_low_noHWE = grm_pot_noHWE$.__enclos_env__$private$GRM$VR_low_noHWE$GRM

#######################################################

## Add sample IDs to row and columns of the GRMs for GBS data. 
colnames(grm_vr_gbs_high_noHWE) = colnames(grm_vr_gbs_low_noHWE) =  grm_pot_noHWE$.__enclos_env__$private$GRM$VR_high_noHWE$indID
rownames(grm_vr_gbs_high_noHWE) = rownames(grm_vr_gbs_low_noHWE) =  grm_pot_noHWE$.__enclos_env__$private$GRM$VR_high_noHWE$indID

############################################################

methods = c("CHIP","GUS (High)","GUS (Low)","GUS_noHWE\n(High)","GUS_noHWE\n(Low)")

## Relatedness comparisons
png("FigureS5.png", width=800*300/72, height=800*300/72, res=300)
noff=sum(upper.tri(grm_vr_ord))
offdiag = function(x) x[upper.tri(x)]
df_rr <- rbind(data.frame(meth="CHIP", item=1:noff, y=offdiag(grm_vr_ord)),
               data.frame(meth="GUS (High)", item=1:noff, y=offdiag(grm_vr_gbs_high)),
               data.frame(meth="GUS (Low)", item=1:noff, y=offdiag(grm_vr_gbs_low)),
               data.frame(meth="GUS_noHWE\n(High)", item=1:noff, y=offdiag(grm_vr_gbs_high_noHWE)),
               data.frame(meth="GUS_noHWE\n(Low)", item=1:noff, y=offdiag(grm_vr_gbs_low_noHWE)))
df_rr$meth = factor(df_rr$meth, levels=c("CHIP","GUS (High)","GUS (Low)","GUS_noHWE\n(High)","GUS_noHWE\n(Low)"))
plot(Meth(df_rr), diff.range = 1.25)
mina = 0.01
maxa = 0.13
dista = (1-mina*2)/length(methods)
xset = 1/4*dista
xset2 = 1/5*dista
yset = 9/10*dista
yset2 = 1/10*dista
for(i in 1:(length(methods)-1)){
  for(j in (i+1):length(methods)){
    temp_cor = format(round(cor(df_rr$y[which(df_rr$meth == methods[i])], df_rr$y[which(df_rr$meth == methods[j])]),2), nsmall = 2)
    temp_ccc = format(round(epiR::epi.ccc(df_rr$y[which(df_rr$meth == methods[i])], df_rr$y[which(df_rr$meth == methods[j])])$rho.c[1],2),nsmall=2)
    text(paste0("r = ",temp_cor), x=mina + (i-1)*dista + xset, y=mina + (length(methods)-j)*dista + yset, xpd=NA)
    text(paste0("ccc = ",temp_ccc), x=maxa + (i-1)*dista - xset2, y=mina + (length(methods)-j)*dista + yset2, xpd=NA)
  }
}
dev.off()

## Self Relatedness comparisons
png("FigureS6.png", width=800*300/72, height=800*300/72, res=300)
ndiag=length(diag(grm_vr_ord))
df_sr <- rbind(data.frame(meth="CHIP", item=1:ndiag, y=diag(grm_vr_ord)),
               data.frame(meth="GUS (High)", item=1:ndiag, y=diag(grm_vr_gbs_high)),
               data.frame(meth="GUS (Low)", item=1:ndiag, y=diag(grm_vr_gbs_low)),
               data.frame(meth="GUS_noHWE\n(High)", item=1:ndiag, y=diag(grm_vr_gbs_high_noHWE)),
               data.frame(meth="GUS_noHWE\n(Low)", item=1:ndiag, y=diag(grm_vr_gbs_low_noHWE)))
df_sr$meth = factor(df_sr$meth, levels=c("CHIP","GUS (High)","GUS (Low)","GUS_noHWE\n(High)","GUS_noHWE\n(Low)"))
plot(Meth(df_sr), diff.range = 1.25)
mina = 0.01
maxa = 0.13
dista = (1-mina*2)/length(methods)
xset = 1/4*dista
xset2 = 1/5*dista
yset = 9/10*dista
yset2 = 1/10*dista
for(i in 1:(length(methods)-1)){
  for(j in (i+1):length(methods)){
    temp_cor = format(round(cor(df_sr$y[which(df_sr$meth == methods[i])], df_sr$y[which(df_sr$meth == methods[j])]),2), nsmall = 2)
    temp_ccc = format(round(epiR::epi.ccc(df_sr$y[which(df_sr$meth == methods[i])], df_sr$y[which(df_sr$meth == methods[j])])$rho.c[1],2),nsmall=2)
    text(paste0("r = ",temp_cor), x=mina + (i-1)*dista + xset, y=mina + (length(methods)-j)*dista + yset, xpd=NA)
    text(paste0("ccc = ",temp_ccc), x=maxa + (i-1)*dista - xset2, y=mina + (length(methods)-j)*dista + yset2, xpd=NA)
  }
}
dev.off()
