library('ape')
library('phangorn')
library('zoom')
library('phytools')
library('adephylo')
library('phylobase')
library('phyloTop')
library('castor')

args <- commandArgs(trailingOnly = TRUE)


filename='timetree_com1.nexus'
filename <- args[[1]]
outname <- args[[2]]
run_number <- args[[3]]
#run_number <- 0
#filename="Outputs/phylogenetic_Run1.newick"
#print(run_number)
#outname <- 'phylogenetic_community1.txt'
#tree_unworked <- read.nexus(filename)
tree_unworked <- read.tree(filename);
# plotTree.singletons(tree)
rooted<-multi2di(tree_unworked, tol=1e-3)
tree<-collapse.singles(tree_unworked,root.edge=TRUE)
distances <- distRoot(tree,tips="all", method="patristic")

#####
# Topology related statistics
####
descendants <- balance(tree)
colless <- sum(abs(descendants[,2]-descendants[,1]))
sackin <- sum(distRoot(tree,tips="all", method="patristic"))

depths <- getDepths(tree)
max_depth <- max(depths$tipDepths)

widths <- integer(max_depth)
rlen_tips <- rle(sort(depths$tipDepths))
rlen_nodes <- rle(sort(depths$nodeDepths))

for (i in 1:max_depth)
{
  if (i %in% rlen_tips$values)
  {
    matchi <- match(i, rlen_tips$values)
    tipswidth <- rlen_tips$lengths[matchi]
    widths[i] <-  tipswidth
  }
  
  if (i %in% rlen_nodes$values)
  {
    matchi <- match(i, rlen_nodes$values)
    nodeswidth <- rlen_nodes$lengths[matchi]
    widths[i] <- widths[i] + nodeswidth
  }  
}



max_width <-max(widths)

WD_ratio <- max_width/max_depth

DeltaW <- max(abs(diff(widths)))

max_ladder <- max(ladderSizes(tree)$ladderNodes)

staircaseness_1 <- length(which((descendants[,1]-descendants[,2])!=0))/length(descendants[,1])

staircaseness_2 <- mean(pmax(descendants[,1],descendants[,2])/pmin(descendants[,1],descendants[,2]))


######
# Summary statistics based on branch lengths, see 
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005416
####
max_H <- max(distances)
min_H <- min(distances)
a_BL_mean <- mean(distances)
a_BL_median <- median(distances)
a_BL_var <- var(distances)


trees <- split_tree_at_height(tree,height=round(max_depth/3), by_edge_count = TRUE)
i_BL_mean_i <- 0
i_BL_var_i <- 0
i_BL_mean_e <- 0
i_BL_var_e <- 0

Ntot_i <- 0
Ntot_e <- 0
for (i in seq(1,trees$Nsubtrees))
{
  subtree = trees$subtrees[[i]]$tree
  
  internal_branches = subtree$edge.length[subtree$edge[,2]>Ntip(subtree)]
  external_branches = subtree$edge.length[subtree$edge[,2]<Ntip(subtree)]
  print(mean(internal_branches))
  if (length(internal_branches)>0){
    i_BL_mean_i = i_BL_mean_i + sum(internal_branches) 
    i_BL_var_i = i_BL_var_i + sum(internal_branches**2)
    Ntot_i = Ntot_i +length(internal_branches)
  }
  if (length(external_branches)>0){
    i_BL_mean_e = i_BL_mean_e + sum(external_branches)
    i_BL_var_e = i_BL_var_e + sum(external_branches**2)
    Ntot_e = Ntot_e + length(external_branches)
  }
}
i_BL_mean_e = i_BL_mean_e/Ntot_e
i_BL_mean_i = i_BL_mean_i/Ntot_i

i_BL_var_e = i_BL_var_e/Ntot_e
i_BL_var_i = i_BL_var_i/Ntot_i

i_BL_var_i <- i_BL_var_i - i_BL_mean_i**2
i_BL_var_e <- i_BL_var_e - i_BL_mean_e**2


#tree = generate_random_tree(list(birth_rate_intercept=1), max_tips=1000)$tree
results = count_lineages_through_time(tree, Ntimes=100,include_slopes = T)

# plot classical LTT curve
#plot(results$times, results$lineages, type="l", xlab="time", ylab="# clades")
# }

max_L <- max(results$lineages)

tmax_arg <- which.max(results$lineages)

t_max_L <- results$times[tmax_arg]

slope_1 <- (max_L - results$lineages[1])/(t_max_L-results$times[1])

slope_2 <- (results$lineages[length(results$lineages)] - max_L)/(results$times[length(results$lineages)] - t_max_L)



####
# LTT stats
####

# calculate classical LTT curve
results = count_lineages_through_time(tree, Ntimes=100,include_slopes=TRUE)

max_L <- max(results$lineages)

tmax_arg <- which.max(results$lineages)

tmax <- results$times[tmax_arg]

t_max_L <- results$times[tmax_arg]

slope_1 <- (max_L - results$lineages[1])/(t_max_L-results$times[1])

slope_2 <- (results$lineages[length(results$lineages)] - max_L)/(results$times[length(results$lineages)] - t_max_L)

slope_ratio <- slope_1/slope_2




# Still missing: i_BL_median_i ie_BL_median_e, IL_nodes, mean_s_time, mean_b_time

tosave_branch <- c(run_number,max_H,min_H,a_BL_mean,a_BL_median, a_BL_var, i_BL_mean_i, i_BL_var_i,
            i_BL_mean_e, i_BL_var_e)

tosave_topology <- c(colless,sackin, WD_ratio, DeltaW, max_ladder, staircaseness_1,
                     staircaseness_2)

tosave_ltt <- c(max_L, t_max_L,slope_1,slope_2,slope_ratio)


savethis <- c(tosave_branch, tosave_topology,tosave_ltt)

savethis <- gsub(" ", "", savethis)
names <- c('#run','max_H','min_H','a_BL_mean','a_BL_median','a_BL_var','i_BL_mean_i','i_BL_var_i','i_BL_mean_e','i_BL_var_e',
           'colless','sackin','WD_ratio','DeltaW','max_ladder','staircaseness_1','staircaseness_2',
           'max_L','t_max_L','slope_1','slope_2','slope_ratio')
df <- data.frame(names,savethis)
#write.table(t(df),file=outname, col.names = FALSE, row.names=FALSE, quote = FALSE, append=TRUE)
write(savethis, outname,append=TRUE,ncolumns=length(names),sep=",")


