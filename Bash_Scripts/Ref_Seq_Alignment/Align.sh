## this line of prompt is written for alignment using mafft
## this output an aligned fasta file down tree making
## ensure mafft is loaded

module load mafft
mafft fasta_file.fasta > fasta.align
