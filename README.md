# Brazilian_Monkey
Fecal bacterial communities of wild black capuchin monkeys (*Sapajus nigritus*) from the Atlantic Forest biome in Southern Brazil are divergent from those of other non-human primates
## Directory organisation
This is our project architecture:
* `data`: input files
* `sources`: code files
* `README.md`: project description

### Sources
* `Monkey_phyloseq`: R code to generate composition plots (phylum, family and genera), alpha-diversity (boxplots, ANOVA and linear tests) and beta-diversity (PCA, clustering and permanova test);
* `OTUs_stat_Deseq2`:  R code to identify if there are  statistic differences between the two regions assessed in the study (Santa Cruz do Sul and São Sebastião do Cai) for phylum and family OTUs; 
* `picrust.sh`: this code is used to predict the functional potential of the bacterial community on the basis of 16S rRNA gene sequencing and the complete tutorial can be found [here](https://github.com/picrust/picrust2). To run this code, you need to install a conda environment providing picrust2-2.3.0_b.
* `functions_stat_Deseq2`: R code to identify if there are statistic differences between the two regions assessed in the study (Santa Cruz do Sul and São Sebastião do Cai) for the predicted functions.
