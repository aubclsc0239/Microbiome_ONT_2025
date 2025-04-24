# this command prompt was written to build a phylogenetic tree using Fasttree
# could be done on the terminal or submitted in the queue depending on data size
# the command requires an aligned fasta file

module load FastTree
FastTree -nt fasta_align > tree.nwk
