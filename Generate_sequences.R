library('ape')
library('phangorn')
filename='/home/fra/Desktop/Oxford/popart/Coatran/CoaTran/phylogenetic.newick'
tree <- ape::read.tree(filename);
sequences = simSeq(tree,l=50)
d = unlist(strsplit(tree$tip.label,'|',fixed=TRUE))
tree$tip.label = rep('x', length(d))#d[seq(2,length(d),3)]
#plot(tree)