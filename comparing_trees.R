library("ape")
library("ggplot2")
library("ggtree")


treefg <- "PARAMS_COMMUNITY1_ACCEPTED/Outputs/phylogenetic_Run1.newick"
tree <- read.tree(file=treefg)
ggtree(tree,mrsd="2018-12-12")+ theme_tree2()





results <- read.csv(file = 'phylogenetic_community1.txt',sep=',',stringsAsFactors = TRUE,)


stats <- read.csv(file='PARAMS_COMMUNITY1_ACCEPTED/Outputs/statistics_runs.csv')

par(mfrow=c(2,2))

for (i in seq(2,length(stats))){
  min <-min( min(stats[,i]),results[,i])
  max <- max (max(stats[,i]),results[,i])
  plot(stats[,i], ylab=names(stats)[i],  ylim=(c (min,max)))+
  abline(h=results[,i], aes(linestyle="dashed", color="red"))
}

#sp <- ggplot(data=stats, aes(data[,1],data[,i])) +geom_point()
#sp <- sp + geom_hline(data=stats, yintercept=results[,i],color="red", linetype="dashed")


params_sex <- read.csv(file = 'PARAMS_COMMUNITY1_ACCEPTED/param_processed_patch0_partnerships.csv',sep=' ',stringsAsFactors = TRUE,)

assortativity<-head(params_sex$assortativity,500)
c_multiplier <-head(params_sex$c_multiplier,500)
brakup_scale<- head(params_sex$breakup_scale_multiplier_overall,500)




#for (branch in tosave_branch){
#  plot(assortativity,tosave_branch)  
#}
