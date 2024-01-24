
## Simulation parameters
run = as.integer(commandArgs(trailingOnly=T))[1]
ploid <- as.integer(commandArgs(trailingOnly=T))[2]
meanDepth <- as.integer(commandArgs(trailingOnly=T))[3] 
#alpha =  1/14    # Double reduction parameter
epsilon = 0.001 #as.integer(commandArgs(trailingOnly=T))[4] # sequencing error parameter
nSires <- 20  # Number of sires
nDams <- 5      # Number of dams per sire
nProg <- 2      # Number of progeny for each pair of parents
nSelf <- 1      # Number of progeny from mating of siblings
nInd <- nSires+(nDams*nSires)*(1+nProg)+ nSelf*choose(nProg,2)*nDams*nSires # Number of individuals in total
seqEffort = 10000000*4
nSnps <- round(seqEffort/(meanDepth*nInd))  # number of SNPs
set = paste0("ploid",ploid,"_depth",meanDepth)
seed = round(run*39+meanDepth*17+ploid*111+nSnps*113)

## filenames

locfolder = "./"
pedfile = paste0(locfolder,"Potato_sim",run,"_",set,".ped")
pfreqfile = paste0(locfolder,"pfreq_sim",run,"_",set,".txt")
genfile = paste0(locfolder,"sim",run,"_",set,".gen")
chromfile = paste0(locfolder,"sim",run,"_",set,".chrom")
parfile = paste0(locfolder,"sim",run,"_",set,".par")
mapfile = paste0(locfolder,"sim",run,"_",set,".map")
outfile = paste0(locfolder,"out_sim",run,"_",set)
###############################
#### Generate the pedigree file:

## Genomic parameters
nChroms = 12
chr_len <- c(98,77,92,113,114,76,72,70,102,75,77,78)
centro <- c(35,5,10,35,50,10,25,20,40,22,45,30)
chr_names <- paste0("Chr",1:12)

## simulate the founder allele frequencies
set.seed(seed)
p <- runif(n = nSnps, min=0, max = 1)
write(p, file=pfreqfile, ncolumns = 1)

#################################################
### Specify pedigree structure and write to file:
nFound = (nSires+nDams*nSires)*2
nPar = nSires + (nSires*nDams)
simped = matrix(nrow=nInd+nFound, ncol=3)
simped[,1] = c(paste0("Founder",1:nFound),paste0("Sire",1:nSires),paste0("Dam",1:(nDams*nSires)),paste0("Offspring",1:(2*nDams*nSires)), paste0("Inbred",1:(nDams*nSires)))

# pedigree of the parents from the founders
simped[1:nPar + nFound,2:3] = simped[1:nFound,1] 

# parents of progeny
progpos <- 1:(nProg*nDams*nSires) + nPar + nFound
simped[progpos,2] = simped[as.vector(sapply(1:nSires,rep,times=nDams*nProg)),1]
simped[progpos,3] = simped[as.vector(sapply(1:(nDams*nSires) + nSires,rep,times=nProg)),1]

# parents of the inbred individuals
selfpos = 1:(nDams*nSires) + nProg*nDams*nSires + nPar + nFound
simped[selfpos,2:3] = matrix(simped[progpos,1], ncol=2, byrow=T) 

# write out the pedigree file
colnames(simped) <- c("Name","Parent1","Parent2")
write.table(simped, file=pedfile, row.names=F, sep="\t", quote=F)

#################################################
## Generate files for PedigreeSim to undertake the simulation
sim_found = t(sapply(p, function(x) rbinom(n=nFound*ploid*2, size=2, prob=x)))
samGeno = cbind(paste0("Mar",1:nSnps), sim_found)
colnames(samGeno) = c("marker", as.vector(sapply(1:nFound, function(x) paste0(simped[x,1],"_",1:(ploid*2)))))
write.table(samGeno, file=genfile, quote=FALSE, row.names=F, sep="\t")
  
## linkage map file
nsnp_per_chrom =  table(sample(1:nChroms,prob=chr_len/sum(chr_len), size=nSnps, replace=T)) # round(chr_len/sum(chr_len)*nSnps)
pos = unlist(sapply(1:nChroms, function(x) round(seq(0,chr_len[x], length.out = nsnp_per_chrom[x]+2)[-c(1,nsnp_per_chrom[x]+2)],6)))
linkfile <- cbind(marker=paste0("Mar",1:nSnps), chromosome=rep(chr_names,nsnp_per_chrom),position=pos)
write.table(linkfile, file=mapfile, quote=FALSE, row.names=F, sep="\t")

## Chromosome file
cat(c("chromosome",	"length",	"centromere",	"prefPairing",	"quadrivalents"),"\n",sep="\t", file=chromfile)
for(chr in 1:nChroms){
  cat(c(chr_names[chr],chr_len[chr], round(centro[chr],2), "0.0", "0.5"),"\n", sep="\t", file=chromfile, append=T)
}

## parameter file
parinfo <- paste(paste0("PLOIDY = ",ploid*2),
             paste0("MAPFUNCTION = KOSAMBI"),
             "MISSING = NA",
             paste0("CHROMFILE = ",chromfile),
             paste0("PEDFILE = ", pedfile),
             paste0("MAPFILE = ", mapfile),
             paste0("FOUNDERFILE = ", genfile),
             paste0("OUTPUT = ", outfile,"\n"),
             "ALLOWNOCHIASMATA = 0",
             "NATURALPAIRING = 1",
             sep="\n")
cat(parinfo, file=parfile)

## Run PedigreeSim
system(paste0("java -jar PedigreeSim.jar ",parfile))

## read in the allele dosages ...
geno <- t(read.table(paste0(outfile,"_alleledose.dat"), header=T))
colnames(geno) <- geno[1,]
geno <- geno[-1,]
indx <- which(!(substr(rownames(geno),1,7) %in% c("Founder")))
geno <- geno[indx,]
geno2 <- matrix(as.integer(geno), nrow=nrow(geno), ncol=ncol(geno))
## generate the GBS data
depth <- matrix(rnbinom(length(geno2),mu=meanDepth, size=2*ploid),ncol=nSnps, nrow=nrow(geno2))
Acounts <- matrix(rbinom(nrow(geno2)*ncol(geno2),depth,geno2/(2*ploid)),ncol=ncol(geno2))
genoGBS <- (2*ploid)*Acounts/depth
## simulate sequencing errors
Bcounts <- depth - Acounts
aCountsFinal <- matrix(rbinom(nrow(geno2)*ncol(geno2),Acounts,prob=1-epsilon),ncol=ncol(geno2)) + matrix(rbinom(nrow(geno2)*ncol(geno2),Bcounts,prob=epsilon),ncol=ncol(geno2))
genoGBS_ep <- (2*ploid)*aCountsFinal/depth

#### Construct the GRMs:

## True IBD inheritance:
founderAllele = read.table(paste0(outfile,"_founderalleles.dat"), header=T)
IBD_all_ind = t(founderAllele[,(nFound*ploid*2+2):ncol(founderAllele)])

IBDgmat <- function(IBD, ploid){
  nInd <- nrow(IBD)/ploid
  gmat <- matrix(NA, nrow=nInd, ncol=nInd)
  for(i in 1:nInd){
    mat1 <- IBD[rep(1:ploid + (i-1)*ploid,ploid),]
    for(j in i:nInd){
      mat2 <- IBD[rep(1:ploid + (j-1)*ploid,rep(ploid,ploid)),]
      theta <- sum(mat1==mat2)/(ncol(IBD)*ploid^2)
      gmat[i,j] <- theta*ploid
    }
  }
  gmat[lower.tri(gmat)] <- t(gmat)[lower.tri(gmat)]
  return(gmat)
}
IBDtrue = IBDgmat(IBD_all_ind, ploid*2)

## GUSrelate
library(GUSrelate)
GUS_gmat_err = computeGRM(aCountsFinal, depth-aCountsFinal, ploid=ploid*2, phat=colMeans(genoGBS_ep, na.rm=T)/(ploid*2), ep=epsilon)

indx_mat = matrix(rep(1:nInd,rep(nInd,nInd)), nrow=nInd, ncol=nInd)
output = cbind(Row=indx_mat[lower.tri(indx_mat, diag=T)], Column=t(indx_mat)[lower.tri(indx_mat, diag=T)],
               IBD=IBDtrue[lower.tri(IBDtrue, diag=T)],
               GUS_err=GUS_gmat_err[lower.tri(GUS_gmat_err, diag=T)])
write.csv(output, file=paste0(locfolder, "files/", "run",run,"_ploid",ploid*2,"_depth",meanDepth,"_ep",epsilon,".csv"), row.names=F)

#pfreqfile = paste0(locfolder,"pfreq_sim",run,".txt")
#genfile = paste0(locfolder,"sim",run,".gen")
#chromfile = paste0(locfolder,"sim",run,".chrom")
#parfile = paste0(locfolder,"sim",run,".par")
#mapfile = paste0(locfolder,"sim",run,".map")
#outfile = paste0(locfolder,"out_sim",run)

file.remove(c(pedfile,pfreqfile,genfile,chromfile,parfile,mapfile))
file.remove(c(paste0(outfile,".hsa"), paste0(outfile,".hsb"), paste0(outfile,"_allAlleledose.dat"), paste0(outfile,"_meioticconfigs.dat"),
              paste0(outfile,"_alleledose.dat"),paste0(outfile,"_founderalleles.dat"),paste0(outfile,"_genotypes.dat")))

