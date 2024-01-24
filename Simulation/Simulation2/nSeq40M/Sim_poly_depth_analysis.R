########## Analysis of Simulation results
########## GUSrelate paper
## Author: Timothy Bilton

library(foreach)
library(doParallel)
library(tidyr)
library(dplyr)

ploidLevels = c(1,2,3,4)*2
Depths = c(1:30,35,40)
epsilon = 0.001
nSnps = round(10000000*4/(Depths*420))
nRuns = 500
#nMet = 5

computeEst = function(est1, est2){
  #return(c(mean(est1 - est2), mean((est1-est2)^2), mean((est1-est2)^2)+mean(est1 - est2)^2, cor(est1,est2)))
  return(c(mean((est1-est2)^2), cor(est1,est2)))
}
estvalues = c("MSE","Correlation") #c("Bias","Var","MSE","Correlation")
comb <- function(x, ...) {  
      mapply(rbind,x,...,SIMPLIFY=FALSE)
}

doParallel::registerDoParallel(20)
res_all = list()
for(ploid in ploidLevels){
  for(depth in Depths){
    tmp = foreach(run = 1:nRuns, .combine="c") %dopar% {
      set = paste0("run",run,"_ploid",ploid,"_depth",depth,"_ep",epsilon)
      res = read.csv(paste0("files/",set,".csv")) 
      sr = which(res$Row == res$Column)
      rr = which(res$Row != res$Column)
      out = computeEst(res$IBD, res$GUS_err)
      out = sapply(c("sr","rr"), function(x) {
        computeEst(res$IBD[get(x)], res$GUS_err[get(x)])
      }, simplify=FALSE)
      out = cbind(Value=unlist(out), Relationship=rep(c("Self-Relatedness","Relatedness"),c(2,2)), Run=run) 
      return(list(out))
    }
    res_all = c(res_all,lapply(tmp, function(x) cbind(x, Quantity=estvalues,  ploid=ploid, depth=depth, ep=epsilon, nSnps=nSnps[which(depth == Depths)])))
  }
}

## Combine results across the different simulations
res = do.call("rbind",res_all)

## results to a file
save(res, file="Sim_poly_depth_res.Rdata")




