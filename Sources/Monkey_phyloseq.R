library(ggplot2) 
library(phyloseq) 
library(reshape2) 
library(ape) 
library(scales)
library(vegan) 
library(nlme)
library(dplyr)
library(reshape2)
library(gridExtra) 
library(parallel) 
library(permute) 
library(lattice)

scripts <- c("graphical_methods.R",
             "tree_methods.R",
             "plot_merged_trees.R",
             "specificity_methods.R",
             "ternary_plot.R",
             "richness.R",
             "edgePCA.R",
             "copy_number_correction.R",
             "import_frogs.R",
             "prevalence.R",
             "compute_niche.R")
urls <- paste0("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/R/", scripts)

for (url in urls) {
  source(url)
}

# working directory
wd <- list()

# The path for the data in my working directory
wd$data   <- "C:/Users/myUser/data/"

#add tables
otu_table <- read.csv2(paste(wd$data, "otu_table.csv", sep=""),row.names=1, dec=".")

tax_table <- read.csv2(paste(wd$data ,"tax_table.csv", sep=""),row.names=1, dec=".")

sample_data <- read.csv2(paste(wd$data, "sample_data.csv", sep=""),row.names=1, dec=".")

#convert to matrix
taxM <- as.matrix(tax_table)
otuM <- as.matrix(otu_table)

#convert to phyloseq OTU and TAX
OTU = otu_table(otuM, taxa_are_rows = TRUE)
TAX = tax_table(taxM)
monkey = phyloseq(OTU, TAX) 

#add to phyloseq sample-data
sample_data(monkey) <-sample_data

#add to phyloseq tree
tree <- read_tree(wd$data + "Tree.nwk")
phy_tree(monkey) <- tree

#rarefy_even_depth downsamples/normalizes all samples to the same depth and prunes OTUs that disappear from all samples as a result. 
monkeyRare <- rarefy_even_depth(monkey, rngseed = 1121983) 
# check with sample_sums 
sample_sums(monkeyRare)[1:5]

monkey <- monkeyRare 

#Customization: plot_composition function - Phylum
p <- plot_composition(monkey, "Kingdom", "Bacteria", "Phylum", numberOfTaxa = 5, fill = "Phylum") 
p <- p + scale_fill_manual(values=c("darkslategray4", "darkseagreen", "springgreen4", "darkslategray3", "darkseagreen2", " gray23"))
p <- p + scale_colour_manual(values=c ("darkslategray4", "darkseagreen", "springgreen4", "darkslategray3", "darkseagreen2", " gray23"))
p <- p + facet_wrap(vars(Region), scales = "free_x", nrow = 1) 
plot(p)

#How to select the 20 most abundant Family ? 
p <- plot_composition(monkey, "Kingdom", "Bacteria", "Family", numberOfTaxa = 20, fill = "Family") 
p <- p + scale_fill_manual(values=c("palegreen4", "palegreen3", "palegreen", "goldenrod1", "darkorange", "darkorange3", "orange3",  "red3", "firebrick1","tomato", "lightpink", "hotpink2", "magenta", "magenta4", "mediumorchid1", " mistyrose2", " mistyrose3", " navajowhite2", " navajowhite3", " navajowhite4", " gray23"))
p <- p + scale_colour_manual(values=c("palegreen4", "palegreen3", "palegreen", "goldenrod1", "darkorange", "darkorange3", "orange3", "red3", "firebrick1","tomato", "lightpink", "hotpink2", "magenta", "magenta4", "mediumorchid1", " mistyrose2", " mistyrose3", " navajowhite2", " navajowhite3", " navajowhite4", " gray23"))
p <- p + facet_wrap(vars(Region), scales = "free_x", nrow = 1) 
plot(p)

#How to select the 20 most abundant genera? 
p <- plot_composition(monkey, "Kingdom", "Bacteria", "Genus", numberOfTaxa = 20, fill = "Genus", x = "SampleID") 
p <- p + scale_fill_manual(values=c("palegreen4", "palegreen3", "palegreen", "goldenrod1", "darkorange", "darkorange3", "orange3",  "red3", "firebrick1","tomato", "lightpink", "hotpink2", "magenta", "magenta4", "mediumorchid1", " mistyrose2", " mistyrose3", " navajowhite2", " navajowhite3", " navajowhite4", " gray23"))
p <- p + scale_colour_manual(values=c("palegreen4", "palegreen3", "palegreen", "goldenrod1", "darkorange", "darkorange3", "orange3", "red3", "firebrick1","tomato", "lightpink", "hotpink2", "magenta", "magenta4", "mediumorchid1", " mistyrose2", " mistyrose3", " navajowhite2", " navajowhite3", " navajowhite4", " gray23"))
p <- p + facet_wrap(vars(Region), scales = "free_x", nrow = 1) 
plot(p)

#Alpha-diversity
p <- plot_richness(monkey, color = "Region", x = "Region", measures = c("Observed", "Chao1", "Simpson", "InvSimpson")) 
#plot as boxplot 
p <- p + geom_boxplot(aes(fill = Region), alpha=0.2) + theme(axis.text.x = element_blank())
colors <- c("blue", "red" )
p <- p + scale_color_manual(values = colors)
p <- p + scale_fill_manual(values = colors)
plot(p) 

#statistica
alpha.diversity <- estimate_richness(monkey, measures = c("Observed", "Chao1", "Simpson", "InvSimpson"))

#ANOVA Observed
data <- cbind(sample_data(monkey), alpha.diversity) 
monkey.anova.Observed <- aov(Observed ~ Region, data)  
summary(monkey.anova.Observed) 

#ANOVA Chao1
data <- cbind(sample_data(monkey), alpha.diversity) 
monkey.anova.Chao1 <- aov(Chao1 ~ Region, data)  
summary(monkey.anova.Chao1) 

#ANOVA Simpson
data <- cbind(sample_data(monkey), alpha.diversity) 
monkey.anova.Simpson <- aov(Simpson ~ Region, data)  
summary(monkey.anova.Simpson) 

#ANOVA InvSimpson
data <- cbind(sample_data(monkey), alpha.diversity) 
monkey.anova.InvSimpson <- aov(InvSimpson ~ Region, data)  
summary(monkey.anova.InvSimpson) 

# Performing a linear model regression to assess the effect of the factor Region and the effect of sampling on alpha diversity measures
#Observed
ps_alpha_div <- estimate_richness(monkey, split = TRUE, measure = "Observed")
ps_alpha_div$SampleID <- rownames(ps_alpha_div) %>% as.factor()
ps_samp <- sample_data(monkey) %>% unclass() %>% data.frame() %>% left_join(ps_alpha_div, by = "SampleID") %>% melt(measure.vars = "Observed", variable.name = "diversity_measure", value.name = "alpha_diversity")
diversity_means <- ps_samp %>%  group_by(SampleID) %>% summarise(mean_div = mean(alpha_diversity)) %>% arrange(mean_div)
ps_samp$SampleID <- factor(ps_samp$SampleID, diversity_means$SampleID)
alpha_div_model_Observed <- lme(fixed = alpha_diversity ~ Region, data = ps_samp, random = ~ 1 | SampleID)
summary(alpha_div_model_Observed)

#Chao1
ps_alpha_div <- estimate_richness(monkey, split = TRUE, measure = "Chao1")
ps_alpha_div$SampleID <- rownames(ps_alpha_div) %>% as.factor()
ps_samp <- sample_data(monkey) %>% unclass() %>% data.frame() %>% left_join(ps_alpha_div, by = "SampleID") %>% melt(measure.vars = "Chao1", variable.name = "diversity_measure", value.name = "alpha_diversity")
diversity_means <- ps_samp %>%  group_by(SampleID) %>% summarise(mean_div = mean(alpha_diversity)) %>% arrange(mean_div)
ps_samp$SampleID <- factor(ps_samp$SampleID, diversity_means$SampleID)
alpha_div_model_Chao1 <- lme(fixed = alpha_diversity ~ Region, data = ps_samp, random = ~ 1 | SampleID)
summary(alpha_div_model_Chao1)

#Simpson
ps_alpha_div <- estimate_richness(monkey, split = TRUE, measure = "Simpson")
ps_alpha_div$SampleID <- rownames(ps_alpha_div) %>% as.factor()
ps_samp <- sample_data(monkey) %>% unclass() %>% data.frame() %>% left_join(ps_alpha_div, by = "SampleID") %>% melt(measure.vars = "Simpson", variable.name = "diversity_measure", value.name = "alpha_diversity")
diversity_means <- ps_samp %>%  group_by(SampleID) %>% summarise(mean_div = mean(alpha_diversity)) %>% arrange(mean_div)
ps_samp$SampleID <- factor(ps_samp$SampleID, diversity_means$SampleID)
alpha_div_model_Simpson <- lme(fixed = alpha_diversity ~ Region, data = ps_samp, random = ~ 1 | SampleID)
summary(alpha_div_model_Simpson)

#InvSimpson
ps_alpha_div <- estimate_richness(monkey, split = TRUE, measure = "InvSimpson")
ps_alpha_div$SampleID <- rownames(ps_alpha_div) %>% as.factor()
ps_samp <- sample_data(monkey) %>% unclass() %>% data.frame() %>% left_join(ps_alpha_div, by = "SampleID") %>% melt(measure.vars = "InvSimpson", variable.name = "diversity_measure", value.name = "alpha_diversity")
diversity_means <- ps_samp %>%  group_by(SampleID) %>% summarise(mean_div = mean(alpha_diversity)) %>% arrange(mean_div)
ps_samp$SampleID <- factor(ps_samp$SampleID, diversity_means$SampleID)
alpha_div_model_InvSimpson <- lme(fixed = alpha_diversity ~ Region, data = ps_samp, random = ~ 1 | SampleID)
summary(alpha_div_model_InvSimpson)


#PCoA
ord <- ordinate(monkey, method = "MDS", distance = "jaccard") 
p <- plot_ordination(monkey, ord, color = "Region") 
p <- p + theme_bw() + ggtitle("PCoA - Jaccard") + geom_text(mapping = aes(label = SampleID), size = 4, vjust=1.5) +  geom_point(size = 3) 
colors <- c("blue" ,"red")
p <- p + scale_color_manual(values = colors)
plot(p) 

#Clustering
dist.jac <- distance(monkey, method = "jaccard")
Region <- get_variable(monkey, "Region")
tipColor <- col_factor(colors, levels = levels(Region), Region)(Region)
clust.jac <- as.phylo(hclust(dist.jac, method = "ward.D2"))  
plot(clust.jac, tip.color = tipColor, direction = "downwards",  main = "Ward.D2-Jaccard")

#Beta-diversity - Permanova
set.seed(1)
metadata <- as(sample_data(monkey), "data.frame")
permanova.test <- adonis(distance(monkey, method="jaccard") ~ Region,  data = metadata)
summary(permanova.test)
permanova.test
