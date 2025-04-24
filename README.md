# 16S rRNA Microbiome EPI2ME to Phyloseq

## Project description
Since long reads technology has lower read quality than short reads, and several commercial pipelines or
packages are not designed to handle long reads appropriately, I plan to integrate several analytical
pipelines to process 16S read from ONT GridIon data. EPI2ME is an analytical pipeline designed for processing ONT reads.
However, while this pipeline performs demultiplexing and primer removal steps, merging of reads, followed by
taxonomic assignment and preliminary alpha diversity metrics, the software does not perform other intergative analysis
such as phylogenetic analysis for beta diversity. This porject intend to process EPI2ME output further and integrate into Phyloseq for intended downline analysis such as differential abundance, functional pathway analysis, and flexibility in graph generation. 

My proposed steps are highlighted
below:



## Links to analysis

- [Challenge 4 Analysis](coding_challenge_4/Challenge4.md)

## file tree

```r
fs::dir_tree()
```

```bash
├── Challenge1.R
├── Challenge2.R
├── Challenge3.R
├── Challenge4.Rmd
├── Challenge4.html
├── Challenge_1.R
├── MycotoxinData.csv
├── PLPA6820_2025.Rproj
├── README.md
└── coding_challenge_4
    ├── Challenge4.docx
    └── Rmarkdown.md
```
