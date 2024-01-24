########## Analysis of Simulation results
## Author: Timothy Bilton

library(foreach)
library(doParallel)

ploidLevels = c(1,2,3,4)
Depths = c(5,25,50)
epsilon = c(0,0.001,0.01)
nRuns = 500
nMet = 5

computeEst = function(est1, est2){
  return(c(mean(est1 - est2), mean((est1-est2)^2), mean((est1-est2)^2) - mean(est1 - est2)^2, cor(est1,est2)))
}
estvalues = c("Bias","MSE","Var","Correlation")
comb <- function(x, ...) {  
      mapply(rbind,x,...,SIMPLIFY=FALSE)
}

### Recreate pedigres structure
nSires <- 20    # Number of sires
nDams <- 5      # Number of dams per sire
nProg <- 2      # Number of progeny for each pair of parents
nSelf <- 1      # Number of progeny from mating of siblings
nInd <- nSires+(nDams*nSires)*(1+nProg)+ nSelf*choose(nProg,2)*nDams*nSires # Number of individuals in total
### Specify pedigree structure and write to file:
nFound = (nSires+nDams*nSires)*2
nPar = nSires + (nSires*nDams)
simped = matrix(nrow=nInd+nFound, ncol=3)
simped[,1] = c(paste0("Founder",1:nFound),paste0("Sire",1:nSires),paste0("Dam",1:(nDams*nSires)),paste0("Offspring",1:(2*nDams*nSires)), paste0("Inbred",1:(nDams*nSires)))
# pedigree of the parents from the founders
simped[1:nPar + nFound,2:3] = simped[1:nFound,1] 
# parents of progeny
progpos <- 1:(nProg*nDams*nSires) + nPar + nFound
simped[progpos,2] = simped[as.vector(sapply(1:nSires,rep,times=nDams*nProg))+nFound,1]
simped[progpos,3] = simped[as.vector(sapply(1:(nDams*nSires) + nSires,rep,times=nProg))+nFound,1]
# parents of the inbred individuals
selfpos = 1:(nDams*nSires) + nProg*nDams*nSires + nPar + nFound
simped[selfpos,2:3] = matrix(simped[progpos,1], ncol=2, byrow=T) 
# write out the pedigree file
colnames(simped) <- c("Name","Parent1","Parent2")
simped = simped[-c(1:240),]

## Indexes for the specific relationships
i_dad = paste(match(simped[,2],simped[,1]), 1:nrow(simped), sep=":")
i_mum = paste(match(simped[,3],simped[,1]), 1:nrow(simped), sep=":")

## Index the different relationships:
res = read.csv(paste0("files/run1_ploid1_Depth5_ep0.csv"))

res_ind = paste(res[,1], res[,2], sep=":")

## Self-relatedness All
sr_all = which(res$Row == res$Column)
rr_all = which(res$Row != res$Column)
## Self-relatedness Gen1
sr_g1 = which(res$Row == res$Column & res$Row < 121)
## Self-relatedness Gen2
sr_g2 = which(res$Row == res$Column & res$Row < 321 & res$Row > 120)
## Self-relatedness Gen3
sr_g3 = which(res$Row == res$Column & res$Row > 320)
## relatedness Parent-offspring (G1-G2)
rr_g1 = which(res_ind %in% i_dad[121:320] | res_ind %in% i_mum[121:320])
## relatedness Parent-offspring (G2-G3)
rr_g2 = which(res_ind %in% i_dad[321:420] | res_ind %in% i_mum[321:420])
## relatedness full-sibs (G2)
rr_g3 = which(res_ind %in% paste(seq(121,320,2),seq(122,320,2), sep=":"))
## relatedness grandparents-offspring (G1-G3)
rr_g4 = which(res_ind %in% c(paste(match(simped[,2],simped[,1])[match(simped[,2],simped[,1])[321:420]], 321:420, sep=":"),
                             paste(match(simped[,3],simped[,1])[match(simped[,3],simped[,1])[321:420]], 321:420, sep=":")))
## relatedness half-sibs (G2: different dam, same sire)
rr_g5 = which(res_ind %in% unlist(sapply(121:320, function(x) {
  indPar = ceiling((x-120)/2)
  if((indPar%%5) != 0){
    indOPar  = ((indPar%%5):5 + 5*((indPar-1)%/%5))[-1]
    indHS = 120 + seq(min(indOPar)-0.5, max(indOPar), 0.5)*2
    paste(x,indHS, sep=":")
  } else return(character(0))
})))
## relatedness half-cousins (G3: different granddam, same grandsire)
rr_g6 = which(res_ind %in% as.vector(sapply(seq(321,420,5), function(x) {tmp = t(combn(x:(x+4),2)); paste(tmp[,1],tmp[,2], sep=":")})))
## relatedness Parent-half nibling (G2-G3): Different family 
rr_g7 = which(res_ind %in% sapply(121:320, function(x) {
  indPar = ceiling((x-120)/2)
  indOff = setdiff(seq(320,419,5)[ceiling(indPar/5)] + 1:5,indPar+320)
  paste(x, indOff, sep=":")
}))
## relatedness Unrelated
rr_g8 = setdiff(1:nrow(res), c(sr_g1,sr_g2,sr_g3,rr_g1,rr_g2,rr_g3,rr_g4,rr_g5,rr_g6,rr_g7))


relate_groups = c("sr_all","rr_all","sr_g1","sr_g2","sr_g3","rr_g1","rr_g2","rr_g3","rr_g4","rr_g5","rr_g6","rr_g7","rr_g8")

rm(res)

doParallel::registerDoParallel(20)
res_all = list()
for(ploid in ploidLevels){
  for(depth in Depths){
    for(ep in epsilon){
      tmp = foreach(run = 1:nRuns, .combine="comb") %dopar% {
        set = paste0("run",run,"_ploid",ploid,"_Depth",depth,"_ep",ep)
        res = read.csv(paste0("files/",set,".csv")) 
        out = sapply(relate_groups, function(x) {
          apply(res[get(x),c("Geno","GUS","GUS_err","AGH","Cer")], 2, computeEst, est2=res$IBD[get(x)])
        }, simplify=FALSE)
        return(out)  
      }
    res_all = c(res_all,list(lapply(tmp, function(x) cbind(x, Quantity=estvalues,  ploid=ploid, depth=depth, ep=ep))))
    }
  }
}
## Combine results across the different simulations
res = sapply(1:length(res_all[[1]]), function(x) do.call("rbind",lapply(res_all, function(y) y[[x]])), simplify=F)
names(res) = relate_groups

## results to a file
save(res, file="Sim_poly_res.Rdata")


# ### Plot the results:
library(ggplot2)
res_df = lapply(res, function(x) {
  tmp = as.data.frame(tidyr::pivot_longer(as.data.frame(x, ),  c("Geno","GUS","GUS_err","AGH","Cer"), names_to="Method"))
  tmp$value = as.numeric(tmp$value)
  tmp$ploid = factor(tmp$ploid, levels=ploidLevels)
  tmp$depth = factor(tmp$depth, levels=Depths)
  tmp$ep = factor(tmp$ep, levels=epsilon)
  tmp$Method = factor(tmp$Method, levels=c("Geno","GUS","GUS_err","AGH","Cer"))
  tmp$Quantity = factor(tmp$Quantity, levels=estvalues)
  return(tmp)
})

res_df_mean = lapply(res_df, function(x) cbind(value=as.vector(tapply(x$value,list(x$ploid, x$depth, x$ep, x$Method, x$Quantity), mean, na.rm=T)), expand.grid(ploid=levels(x$ploid), depth=levels(x$depth), ep=levels(x$ep), Method=levels(x$Method), Quantity=levels(x$Quantity))))

res_df_mean = lapply(res_df_mean, function(x){
  index = which(x$Quantity == "MSE")
  x$value[index] = sqrt(x$value[index])
  x$Quantity = as.character(x$Quantity)
  x$Quantity[index] = "RMSE"
  x$Quantity = factor(x$Quantity, levels=unique(sort(x$Quantity)), labels=c("Bias","Correlation","Root Mean Square Error (RMSE)", "Variance"))
  x$Ploidy = factor(as.numeric(as.character(x$ploid))*2)
  x = x[,which(colnames(x) != "ploid")]
  return(x)
})

### Combined MSE:
res_df_mean_comb = do.call("rbind", sapply(3:13, function(x) data.frame(res_df_mean[[x]], Relationship=names(res_df_mean)[x]), simplify = FALSE))
res_df_mean_comb$Relationship = factor(res_df_mean_comb$Relationship, levels=sort(unique(res_df_mean_comb$Relationship)), 
                                       labels=c("PO (G1-G2)","PO (G2-G3)","FS (G2)", "GO (G1-G3)","HS (G2)","HC (G2)","PHN (G2-G3)", "UR","SR (G1)","SR (G2)","SR (G3)"))
## Reorder levels:
res_df_mean_comb$Relationship = factor(res_df_mean_comb$Relationship,
                                       levels = c(rev(c("PO (G2-G3)","PO (G1-G2)","FS (G2)", "GO (G1-G3)","HS (G2)","HC (G2)","PHN (G2-G3)", "UR")),c("SR (G1)","SR (G2)","SR (G3)")))

eptext = paste0("\u03B5 = ", levels(res_df_mean_comb$ep))
names(eptext) = levels(res_df_mean[[1]]$ep)
depthtext = paste0("Depth = ", levels(res_df_mean_comb$depth))
names(depthtext) = levels(res_df_mean[[1]]$depth)

#### MSE
p = ggplot(subset(res_df_mean_comb, res_df_mean_comb$Quantity == "Root Mean Square Error (RMSE)" & res_df_mean_comb$Relationship %in% c("SR (G1)","SR (G2)","SR (G3)")), 
           aes(y=value, x=Relationship, col=Method, shape=Ploidy, group=Method)) + 
  geom_hline(yintercept = 0, linetype="dashed") + #geom_vline(mapping = aes(xintercept = 1.5), linetype="dotted") +
  geom_point(position = position_dodge(width=0.75), size=2) + facet_grid(depth~ep, labeller = labeller(ep = eptext, depth = depthtext)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust=0.5), legend.position = "top") + ylab("Root Mean Square Error (RMSE)")
p = p + geom_vline(xintercept = c(1.5,2.5), linetype="dotted")
ggsave(plot=p, filename=paste0("Figure3.png"), width=8, height=9, dpi=300)

p = ggplot(subset(res_df_mean_comb, res_df_mean_comb$Quantity == "Root Mean Square Error (RMSE)" & !(res_df_mean_comb$Relationship %in% c("SR (G1)","SR (G2)","SR (G3)"))), 
           aes(y=value, x=Relationship, col=Method, shape=Ploidy, group=Method)) + 
  geom_hline(yintercept = 0, linetype="dashed") + #geom_vline(mapping = aes(xintercept = 1.5), linetype="dotted") +
  geom_point(position = position_dodge(width=0.75), size=2) + facet_grid(depth~ep, labeller = labeller(ep = eptext, depth = depthtext)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust=0.5), legend.position = "top") + ylab("Root Mean Square Error (RMSE)")
p = p + geom_vline(xintercept = seq(1.5,7.5,1), linetype="dotted")
ggsave(plot=p, filename=paste0("Figure2.png"), width=10.6, height=7, dpi=300)

#### Bias
p = ggplot(subset(res_df_mean_comb, res_df_mean_comb$Quantity == "Bias" & res_df_mean_comb$Relationship %in% c("SR (G1)","SR (G2)","SR (G3)")), 
           aes(y=value, x=Relationship, col=Method, shape=Ploidy, group=Method)) + 
  geom_hline(yintercept = 0, linetype="dashed") + #geom_vline(mapping = aes(xintercept = 1.5), linetype="dotted") +
  geom_point(position = position_dodge(width=0.75), size=2) + facet_grid(depth~ep, labeller = labeller(ep = eptext, depth = depthtext)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust=0.5), legend.position = "top") + ylab("Bias")
p = p + geom_vline(xintercept = seq(1.5,2.5,1), linetype="dotted")
ggsave(plot=p, filename=paste0("FigureS3.png"), width=8, height=9, dpi=300)

p = ggplot(subset(res_df_mean_comb, res_df_mean_comb$Quantity == "Bias" & !(res_df_mean_comb$Relationship %in% c("SR (G1)","SR (G2)","SR (G3)"))), 
           aes(y=value, x=Relationship, col=Method, shape=Ploidy, group=Method)) + 
  geom_hline(yintercept = 0, linetype="dashed") + #geom_vline(mapping = aes(xintercept = 1.5), linetype="dotted") +
  geom_point(position = position_dodge(width=0.75), size=2) + facet_grid(depth~ep, labeller = labeller(ep = eptext, depth = depthtext)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust=0.5), legend.position = "top") + ylab("Bias")
p = p + geom_vline(xintercept = seq(1.5,7.5,1), linetype="dotted")
ggsave(plot=p, filename=paste0("FigureS1.png"), width=10.6, height=7, dpi=300)

#### Variance
p = ggplot(subset(res_df_mean_comb, res_df_mean_comb$Quantity == "Variance" & res_df_mean_comb$Relationship %in% c("SR (G1)","SR (G2)","SR (G3)")), 
           aes(y=value, x=Relationship, col=Method, shape=Ploidy, group=Method)) + 
  geom_hline(yintercept = 0, linetype="dashed") + #geom_vline(mapping = aes(xintercept = 1.5), linetype="dotted") +
  geom_point(position = position_dodge(width=0.75), size=2) + facet_grid(depth~ep, labeller = labeller(ep = eptext, depth = depthtext)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust=0.5), legend.position = "top") + ylab("Variance")
p = p + geom_vline(xintercept = seq(1.5,2.5,1), linetype="dotted")
ggsave(plot=p, filename=paste0("FigureS4.png"), width=8, height=9, dpi=300)

p = ggplot(subset(res_df_mean_comb, res_df_mean_comb$Quantity == "Variance" & !(res_df_mean_comb$Relationship %in% c("SR (G1)","SR (G2)","SR (G3)"))), 
           aes(y=value, x=Relationship, col=Method, shape=Ploidy, group=Method)) + 
  geom_hline(yintercept = 0, linetype="dashed") + #geom_vline(mapping = aes(xintercept = 1.5), linetype="dotted") +
  geom_point(position = position_dodge(width=0.75), size=2) + facet_grid(depth~ep, labeller = labeller(ep = eptext, depth = depthtext)) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "top") + ylab("Variance")
p = p + geom_vline(xintercept = seq(1.5,7.5,1), linetype="dotted")
ggsave(plot=p, filename=paste0("FigureS2.png"), width=10.6, height=7, dpi=300)


#################################
## table of relationships: (Table 1)
#################################

res = read.csv(paste0("files/run1_ploid1_Depth5_ep0.csv"))

library(AGHmatrix)
nFound = (nSires+nDams*nSires)*2
nPar = nSires + (nSires*nDams)
simped = matrix(nrow=nInd+nFound, ncol=3)
simped[,1] = c(paste0("Founder",1:nFound),paste0("Sire",1:nSires),paste0("Dam",1:(nDams*nSires)),paste0("Offspring",1:(2*nDams*nSires)), paste0("Inbred",1:(nDams*nSires)))

# pedigree of the parents from the founders
simped[1:nPar + nFound,2:3] = simped[1:nFound,1] 

# parents of progeny
progpos <- 1:(nProg*nDams*nSires) + nPar + nFound
simped[progpos,2] = simped[as.vector(sapply(1:nSires,rep,times=nDams*nProg))+nFound,1]
simped[progpos,3] = simped[as.vector(sapply(1:(nDams*nSires) + nSires,rep,times=nProg)) + nFound,1]

# parents of the inbred individuals
selfpos = 1:(nDams*nSires) + nProg*nDams*nSires + nPar + nFound
simped[selfpos,2:3] = matrix(simped[progpos,1], ncol=2, byrow=T) 
simped[is.na(simped)] = 0
simped = as.data.frame(simped)

for(i in 1:ncol(simped))
  simped[,i] <- factor(simped[,i])

Amat <- Amatrix(simped, ploidy = 4, w=0, slater = T)
Amat_sub = Amat[241:ncol(Amat),241:ncol(Amat)]
table(diag(Amat_sub))
table(Amat_sub[upper.tri(Amat_sub)])

table(Amat_sub[as.matrix(res[sr_g1,1:2])])
table(Amat_sub[as.matrix(res[sr_g2,1:2])])
table(Amat_sub[as.matrix(res[sr_g3,1:2])])

table(Amat_sub[as.matrix(res[rr_g1,1:2])])
table(Amat_sub[as.matrix(res[rr_g2,1:2])])
table(Amat_sub[as.matrix(res[rr_g3,1:2])])
table(Amat_sub[as.matrix(res[rr_g4,1:2])])
table(Amat_sub[as.matrix(res[rr_g5,1:2])])
table(Amat_sub[as.matrix(res[rr_g6,1:2])])
table(Amat_sub[as.matrix(res[rr_g7,1:2])])
table(Amat_sub[as.matrix(res[rr_g8,1:2])])

