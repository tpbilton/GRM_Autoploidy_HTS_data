
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())

## parameter values
ploidLevels = c(1,2,3,4)*2
Depths = c(1:30,35,40)
epsilon = 0.001
nRuns = 500

## load results from the simulations
load("nSeq40M/Sim_poly_depth_res.Rdata")
res = data.frame(res, SeqDepth = "40M")
res_all = res

load("nSeq10M/Sim_poly_depth_res.Rdata")
res = data.frame(res, SeqDepth = "10M")
res_all = rbind(res_all,res)

## Reshape data:
res_all$Value = as.numeric(res_all$Value)
res_all$depth = as.numeric(res_all$depth)
res_all$nSnps = as.numeric(res_all$nSnps)
res_df = res_all %>% group_by(SeqDepth, Relationship, Quantity, ploid, depth, nSnps) %>% summarise(Mean = mean(Value))

## Figure 4
p = ggplot(subset(res_df, res_df$Quantity == "MSE"), aes(y=sqrt(Mean), x=depth, col=SeqDepth, linetype=SeqDepth, shape=ploid)) +
  geom_point() + geom_line() + facet_wrap(~Relationship, scales="free_y", ncol=2) +
  ylab("Root Mean Square Error (RMSE)") + xlab("Average Read Depth") +
  labs(shape="Ploidy", col="Total Seq\nEffort", linetype="Total Seq\nEffort")
ggsave(plot = p, filename="Figure4.png", dpi = 300, width=10, height=5) 

## table of minimum depth:
#tab_mindepth = pivot_wider(res_df %>% group_by(SeqDepth, ploid,Relationship) %>% summarise(Min=depth[which.min(Mean)]), names_from = c("Relationship","SeqDepth"), values_from = "Min", names_sort = T)

