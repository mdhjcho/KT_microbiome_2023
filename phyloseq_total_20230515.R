library("phyloseq")
library("vegan")
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(tidyverse)
library(edgeR); packageVersion("edgeR")    
library(limma); packageVersion("limma") 

BiocManager::install("DEFormats")
library(DEFormats); packageVersion("DEFormats") 

library(DESeq2); packageVersion("DESeq2")

BiocManager::install("ANCOMBC")
library(ANCOMBC)
library(Maaslin2)
install.packages("VennDiagram")
library(VennDiagram)

library("ggplot2"); packageVersion("ggplot2")
require(RColorBrewer)
require(randomcoloR)
require(MASS)

sessionInfo()

# Define a default theme for ggplot graphics.
theme_set(theme_bw())


# Install qiime2R
install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")
library(qiime2R)

library(devtools)
install_github("umerijaz/microbiomeSeq") 
library(microbiomeSeq)
library(plyr)

sessionInfo()

getwd()


metadata<-read.table("/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/sample-metadata_20210818.txt",sep="\t", header = T)

physeq_total<-qza_to_phyloseq(features = "table-2424.qza",
                              tree = "rooted-tree.qza", "taxonomy-silva-2424.qza",
                              metadata = "sample-metadata_20210818.txt")

ntaxa(physeq_total) #16171
nsamples(physeq_total) #220

sample_variables(physeq_total)
head(sample_data(physeq_total))


metadata1<-data.frame(sample_data(physeq_total))

# Create table, number of features for each phyla
table(tax_table(physeq_total)[,"Phylum"], exclude=NULL)

# Examining the number of reads for each sample
sample_sums(physeq_total) # to calculate the total number of reads for each sample
sort(sample_sums(physeq_total)) # min read 639, R164_0, R58_12, R80_12, R114_12, R56_12 <1000 read
hist(sample_sums(physeq_total), main="Histogram: Read Counts", xlab="Total Reads", 
     border="blue", col="green", las=1, breaks=12)

metadata1$total_reads <- sample_sums(physeq_total)

# Examining the OTU table
ntaxa(physeq_total)
head(taxa_names(physeq_total))
head(taxa_sums(physeq_total))
asv_tab<-data.frame(otu_table(physeq_total)) # ASV table by sample_id

# examining the taxonomy
rank_names(physeq_total)
head(tax_table(physeq_total))
table(tax_table(physeq_total)[,2])
tax_tab<-data.frame(tax_table(physeq_total)) # Taxa table

# examining the reference sequences
head(refseq(physeq_total))

# Agglomerating and subsetting taxa
pt_phylum<-tax_glom(physeq_total,"Phylum")
pt_phylum
taxa_names(pt_phylum)
taxa_names(pt_phylum) <- tax_table(pt_phylum)[, 2] # tax name ??Á®?À±? 
taxa_names(pt_phylum)
otu_table(pt_phylum)[1:5, c(1:3, 5, 7)] 


# Subset taxa..
(pt_bacteroides <- subset_taxa(physeq_total, Genus == "Bacteroides"))
tax_table(pt_bacteroides)


## Filter ###

# Create table, number of features for each phyla
table(tax_table(physeq_total)[, "Phylum"], exclude = NULL)

# remove ambiguous phylum annotation
physeq_f<- subset_taxa(physeq_total, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
physeq_f
table(tax_table(physeq_f)[, "Phylum"], exclude = NULL)

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(physeq_f),
               MARGIN = ifelse(taxa_are_rows(physeq_f), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(physeq_f),
                    tax_table(physeq_f))

prevdf 


# Compute the total and average prevalences of the features in each phylum.
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Define phyla to filter
filterPhyla = c("Euryarchaeota", "Marinimicrobia_(SAR406_clade)","Myxococcota")

# Filter entries with unidentified Phylum.
physeq1 = subset_taxa(physeq_f, !Phylum %in% filterPhyla)
physeq1

# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(physeq1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(physeq_f),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(physeq_f)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
physeq2 = prune_taxa(keepTaxa, physeq_f)
physeq2 # results of 5% prevalance

# Define prevalence threshold as 10% of total samples
prevalenceThreshold1 = 0.1 * nsamples(physeq_f)
prevalenceThreshold1

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa1 = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold1)]
physeq3 = prune_taxa(keepTaxa1, physeq_f)
physeq3 # result of 10% Prevalance

# Define prevalence threshold as 1% of total samples
prevalenceThreshold2 = 0.01 * nsamples(physeq_f)
prevalenceThreshold2

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa2 = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold2)]
physeq4 = prune_taxa(keepTaxa2, physeq_f)
physeq4 # result of 10% Prevalance


# How many genera would be present after filtering?
length(get_taxa_unique(physeq1, taxonomic.rank = "Genus")) # 328,total data 
length(get_taxa_unique(physeq2, taxonomic.rank = "Genus")) # 61, prevalence 5%
length(get_taxa_unique(physeq3, taxonomic.rank = "Genus")) # 33, prevalence 10%
length(get_taxa_unique(physeq4, taxonomic.rank = "Genus")) # 124, prevalence 1%

## Real data analysis for total data
#1. remove read <1000
(physeq1_1k = prune_samples(sample_sums(physeq1) > 1000, physeq1))
physeq1_1k 

#2.filer, sample >10
physeq1_filter = filter_taxa(physeq1_1k, function(x) sum(x)>10, TRUE)
physeq1_filter 
(physeq1_filter = prune_taxa(taxa_sums(physeq1_filter) > 0, physeq1_filter))
any(taxa_sums(physeq1_filter) == 0)

#3. rarefying
?rarefy_even_depth
(ps_rare <- rarefy_even_depth(physeq1_filter, sample.size = 1500, rngseed = 711, replace = FALSE))
# R1_3, R111_12, R119_12, R53_12, R65_12 removed d/t low sample sized


#4. TSS
ps_relabund<-transform_sample_counts(ps_rare, function(x) x / sum(x))
ps_relabund

par(mfrow = c(1, 2))
title = "Sum of reads for each sample, ps_rare"
plot(sort(sample_sums(ps_rare), TRUE), type = "h", main = title, ylab = "reads", 
     ylim = c(0, 2000))
title = "Sum of reads for each sample, ps_relabund"
plot(sort(sample_sums(ps_relabund), TRUE), type = "h", main = title, ylab = "reads", 
     ylim = c(0, 2000))

# Fix month levels in sample_data
sample_data(ps_rare)$time <- factor(
  sample_data(ps_rare)$time, 
  levels = c("0m", "3m", "12m"))

# alpha-diversity
rich_T<-estimate_richness(ps_rare)
rich_T

desired_order<-c("0m","3m","12m")
(p <- plot_richness(ps_rare, x = "time", color = "time", measures = c("Observed", "Shannon")))
p$data$time<-as.character(p$data$time)
p$data$time<-factor(p$data$time, levels=desired_order)
p+geom_boxplot()
p


kruskal.test(rich_T$Observed~sample_data(ps_rare)$time)
pairwise.wilcox.test(rich_T$Observed,sample_data(ps_rare)$time,p.adjust.method="fdr")
? kruskal.test
? pairwise.wilcox.test

kruskal.test(rich_T$Shannon~sample_data(ps_rare)$time)
pairwise.wilcox.test(rich_T$Shannon,sample_data(ps_rare)$time,p.adjust.method="fdr")

# Violin plot of alph diversity
sample_colors <- c("0m"="lightgoldenrod", "3m"="darkorange","12m"="lightblue")
sample_types <- c("0m","3m","12m")
sample_labels <- c("0m"="0 month","3m"="3 month","12m"="12 month")

desired_order<-c("0m","3m","12m")
p$data$time<-as.character(p$data$time)
p$data$time<-factor(p$data$time, levels=desired_order)


 # Vionlin plot
p<-plot_richness(ps_rare,
                 x = "time",   #compare diversity of datatype
                 measures = c("Observed", "Shannon")) +  #choose diversity measures
  geom_violin(aes(fill = time), show.legend = FALSE)+ #make violin plot, set fill aes to Time
  geom_boxplot(width=0.1) + #add boxplot, set width
  theme_classic()+       #change theme to classic
  xlab(NULL)+    #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),  #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 15))+  #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+ #adjust headings
  scale_fill_manual(values = sample_colors)+   #set fill colors
  scale_x_discrete(          #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +        #add title
  theme(plot.title=element_text(size = 15, face = "bold", hjust = 0.5)) #change title size, face and position
p



# beta diversity: ???PCoA???on bray, plot type 2
#ordination
all_pcoa <- ordinate(
  physeq = ps_rare, 
  method = "PCoA", 
  distance = "bray")

p<-plot_ordination(
  physeq = ps_rare,           #phyloseq object
  ordination = all_pcoa)+     #ordination
  geom_point(aes(fill = time, shape = time), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22, 23))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

# beta diversity: ???NMDS???on bray, plot type 2

set.seed(1) # Let???s try an NMDS instead. For NMDS plots it???s important to set a seed since the starting positions of samples in the alogrithm is random.

#ordination
all_nmds <- ordinate(
  physeq = ps_rare, 
  method = "NMDS", 
  distance = "bray")

p<-plot_ordination(
  physeq = ps_rare,           #phyloseq object
  ordination = all_nmds)+     #ordination
  geom_point(aes(fill = time, shape = time), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22, 23))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

## Calculate bray curtis distance matrix: permanova
set.seed(1)
?phyloseq::distance
# Calculate bray curtis distance matrix
ps_bray<- phyloseq::distance(ps_rare, method = "bray")

# Adonis test
library(phyloseq)
library(vegan)
adonis(ps_bray ~ sample_data(ps_rare)$time) # 0.001

# Homogeneity of dispersion test
beta <- betadisper(ps_bray, sample_data(ps_rare)$time)
permutest(beta) # p 0.003 **


## 0 month analysis

# Subsetting samples and transforming counts
ps_0m<-subset_samples(physeq1, time =="0m")
ps_0m


#1. remove read <1000
(ps_0m_1k = prune_samples(sample_sums(ps_0m) > 1000, ps_0m))

#2.filer, sample >10
ps_0m_filter = filter_taxa(ps_0m_1k, function(x) sum(x)>10, TRUE)
ps_0m_filter
(ps_0m_filter = prune_taxa(taxa_sums(ps_0m_filter) > 0, ps_0m_filter))
any(taxa_sums(ps_0m_filter) == 0)

summarize_phyloseq(ps_0m_filter)


# Microbiome package analysis========================

#1. data operation

# Absolute abundances
otu.absolute <- abundances(ps_0m_filter)

# Relative abundances
otu.relative <- abundances(ps_0m_filter, "compositional")

# Taxonomy table
tax_tab_0m_p <- phyloseq::tax_table(ps_0m_filter)

# check 
tax_tab_0m_p[1:5,1:5]

# Number of taxa
n <- ntaxa(ps_0m_filter)

# Most abundant taxa
topx <- top_taxa(ps_0m_filter, n = 10)


#.2 Cleaning Taxonomy table

install.packages("devtools")
devtools::install_github("microsud/microbiomeutilities")
library(microbiomeutilities)

# Next, cleaning the "k__" and similar values.
# extending gsub to entire table 
tax_table(ps_0m_filter)[, colnames(tax_table(ps_0m_filter))] <- gsub(tax_table(ps_0m_filter)[, colnames(tax_table(ps_0m_filter))],     pattern = "[a-z]__", replacement = "")

colnames(tax_table(ps_0m_filter))
phyloseq::tax_table(ps_0m_filter_ex)[1:3,1:3] 

# Sequence -> ASV1 
ps_0m_filter_ex<-ps_0m_filter # replace sequence -> ASV1, ASV2...
taxa_names(ps_0m_filter_ex) <- paste0("ASV", seq(ntaxa(ps_0m_filter_ex)))

# There can also be spaces between names e.g [Prevotella] copri which can be changed to Prevotella copri (Failure, Not done)

taxa_names(ps_0m_filter_ex) <- gsub("\\[|\\]", "", taxa_names(ps_0m_filter_ex))

# Aggregate at genus level.
ps.gen <- phyloseq::tax_glom(ps_0m_filter_ex, "Genus", NArm = TRUE)
taxa_names(ps.gen)[1:5]

# Substitute these IDs with names of genus.
taxa_names(ps.gen) <- tax_table(ps.gen)[,"Genus"]

tax_table(ps.gen)[,"Genus"][1:5]

unique(tax_table(ps.gen)[,"Genus"] )
ps.gen
ps_0m_filter_ex


# Joining otu/asv table and taxonomy in one data frame
asv_tab_0m <- as.data.frame(abundances(ps_0m_filter_ex)) # get asvs/otus
asv_tab_0m$asv_id <- rownames(asv_tab_0m) # add a new column for ids
#tax_tab <- as.data.frame(tax_table(x)) # get taxonomy note: can be slow
tax_tab_0m <- as(ps_0m_filter_ex@tax_table,"matrix") # get taxonomy note as matrix
tax_tab_0m <- as.data.frame(tax_tab_0m) # convert to data frame
tax_tab_0m$asv_id <- rownames(tax_tab_0m) # add a new column for ids
asv_tax_tab_0m <- tax_tab_0m %>% 
  left_join(asv_tab_0m, by="asv_id") # join to get taxonomy and asv table

head(asv_tax_tab_0m)[,1:8]  # ASV table of 0m, annotation by taxonomy  


#3. Microbiome composition

install.packages("hrbrthemes")
library(hrbrthemes)
install.packages("gcookbook")
library(gcookbook)
library(tidyverse)
library(dplyr)
BiocManager::install("jeevanuDB")
library(jeevanuDB)
library(scales)
library(RColorBrewer)
install.packages("colorRamps")
library(colorRamps)


# Make sure we use functions from correct package
transform <- microbiome::transform


# Merge rare taxa to speed up examples
otu.relative_ex <- transform(ps_0m_filter_ex, "compositional")
pseq <- aggregate_rare(otu.relative_ex, level = "Genus", detection = 1/100, prevalence = 30/100)
pseq

sample_data(pseq)$Rejection_1YA <- factor(
  sample_data(pseq)$Rejection_1YA, 
  levels = c("No Rejection","Rejection"))

tax_pseq <- as(pseq@tax_table,"matrix") # get taxonomy note as matrix
tax_pseq  <- as.data.frame(tax_pseq ) #

colourCount <- length(unique(tax_pseq$Genus)) # number of levels
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

# Figure 1c, Averaged by group of rejection ===================

f1c <- plot_composition(pseq,
                      average_by = "Rejection_1YA") + 
  guides(fill = guide_legend(ncol = 1)) + 
  labs(
       y = "Relative abundance") 
print(f1c + scale_fill_manual("Genus", values = getPalette(colourCount)) + theme_bw())


ggsave("figure1c.pdf")
dev.off()

getwd()

#4. Figure 1d, venn diagram ======================
install.packages("eulerr")
library(eulerr)
sessionInfo()

# simple way to count number of samples in each group
table(meta(ps_0m_filter_ex)$Rejection_1YA, useNA = "always") # NR 64, R 33

# convert to relative abundances
otu.relative_ex<- transform(ps_0m_filter_ex, "compositional")

# Make a list of DiseaseStates.
Rejection_1YA <- unique(as.character(meta(otu.relative_ex)$Rejection_1YA))
print(Rejection_1YA)


# use the pseq.rel object created at the begening of this tutorial. 
taxa_names(otu.relative_ex)[1:5]


# format names
pseq.rel.f <- format_to_besthit(otu.relative_ex)
pseq.f <- format_to_besthit(ps_0m_filter_ex)
pseq.f

# check names
taxa_names(pseq.rel.f)[1:5]

# Write a for loop to go through each of the Rejection_1YA  one by one and combine identified core taxa into a list.

list_core <- c() # an empty object to store information

for (n in Rejection_1YA){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel.f, Rejection_1YA == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.3)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

print(list_core)
class(list_core)
save(list_core, file="list_core_vendiagram.RData")
capture.output(list_core, file = "list_core_vendiagram.csv")

# Specify colors and plot venn
# supplying colors in the order they appear in list_core
mycols <- c(NoRejection="#FF99FF", Rejection="#99CCFF") 

f1d<-plot(venn(list_core),
          fills = mycols)

f1d

ggsave("figure1d.pdf")
dev.off()


#6. prevalence
install.packages("knitr")
library(knitr)
library(microbiome)

pseq.rel.f
head(prevalence(pseq.rel.f, detection = 0.001, sort=TRUE))
head(prevalence(pseq.rel.f, detection = 0.3, sort=TRUE))

### core micobiota(detection=0.1%, prevalence 50%)

pseq.core <- core(pseq.rel.f, detection = 0.001, prevalence = 0.5)
# how many taxa are there?
taxa(pseq.core) %>% length()

# vizualize core microbiota

# core line plots
detections <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(0.05, 1, 0.05)
plot_core(pseq.rel.f, prevalences = prevalences, detections = detections, plot.type = "lineplot")


# core heatmaps
prevalences <- seq(0.05, 1, 0.05)
detections <- round(10^seq(log10(5e-3), log10(0.2), length = 10),3)
# define gray color palette
gray <- gray(seq(0, 1, length = 5))

p_core <- plot_core(
  pseq.rel.f, 
  plot.type = "heatmap", 
  colours = rev(brewer.pal(5, "Spectral")),
  prevalences = prevalences, 
  detections = detections,
  min.prevalence = 0.1)+
  theme(axis.text.x = element_text(size = 9))


print(p_core)
p_core
getwd()
ggsave("figure_coreplot.pdf")
dev.off()


# Core with absolute counts and horizontal view:
# and minimum population prevalence (given as percentage)
detections <- 10^seq(log10(1), log10(max(abundances(ps_0m_filter_ex))/100), length = 10)
library(RColorBrewer)
p_core2 <- plot_core(ps_0m_filter_ex, 
               plot.type = "heatmap",
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .1, 
               horizontal = TRUE)
print(p_core2)



#7.diffrenetial analysis-deseq2 (Not done)
library(DESeq2)
# start by converting phyloseq object to deseq2 format
ds2 <- phyloseq_to_deseq2(ps_0m_filter_ex, ~Rejection_1YA)
# Run DeSeq2 analysis (all taxa at once!)
dds <- DESeq(ds2) # error

# Compare nationalities based on DESeq2 results. Note that covariates are ignored.

res1 <- results(dds, contrast = c("Rejection_1YA", "G_Age", "sex"))
kable(head(as.data.frame(res1)))



#8. ACOMBC
pseq.f
summarize_phyloseq(pseq.f)

# perform the analysis 
ancombc_out = ancombc(
  phyloseq = pseq, 
  formula = "pheno", 
  p_adj_method = "fdr", 
  zero_cut = 1, # no prev filtering necessary anymore 
  lib_cut = 0, 
  group = "pheno", 
  struc_zero = TRUE, 
  neg_lb = TRUE, 
  tol = 1e-5, 
  max_iter = 100, 
  conserve = TRUE, 
  alpha = 0.05, 
  global = TRUE
)
# store the results in res 
res <- out$res



# Microbial pakagomeMarker ==c=======



# wit
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("microbiomeMarker")h library(microbiomeMarker)

base R
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa") 

## Agglomerate Taxa


ps_0m_filter1=  tax_glom(ps_0m_filter,  "Genus",  NArm  =  TRUE)

h1  =  0.4
ps_0m_filter2 =  tip_glom(ps_0m_filter,  h  =  h1)

multiPlotTitleTextSize  =  15
p2tree  =  plot_tree(ps_0m_filter,  method  =  "treeonly",
                     ladderize  =  "left",
                     title  =  "Before  Agglomeration")  + theme(plot.title  =  element_text(size  =  multiPlotTitleTextSize))
p3tree  =  plot_tree(ps_0m_filter1,  method  =  "treeonly",
                     ladderize  =  "left",  title  =  "By  Genus")  + theme(plot.title  =  element_text(size  =  multiPlotTitleTextSize))
p4tree  =  plot_tree(ps_0m_filter2,  method  =  "treeonly",
                     ladderize  =  "left",  title  =  "By  Height")  + theme(plot.title  =  element_text(size  =  multiPlotTitleTextSize))

grid.arrange(nrow  =  1,  p2tree,  p3tree,  p4tree)



#3. rarefying(for beta-diversity)
?rarefy_even_depth
(ps_0m_rare <- rarefy_even_depth(ps_0m_filter, sample.size = 1500, rngseed = 711, replace = FALSE))

#4. TSS
ps_relabund_0m_rare<-transform_sample_counts(ps_0m_rare, function(x) x / sum(x))
ps_relabund_0m_rare # rarefying data

ps_relabund_0m<-transform_sample_counts(ps_0m_filter, function(x) x / sum(x))
ps_relabund_0m # Non rarefying data

otu_table(ps_relabund_0m)[1:10, 1:10]

#3-1. Transform to relative abundance. Save as new object.
ps1.rel  =  transform_sample_counts(ps_0m_filter1,  function(x){x  /  sum(x)}) # relative abundance of agglomerate data

ps1.rel 

?psmelt
ps1ra_melt  <-  psmelt(ps1.rel)  ## data of Agglomerate Taxa

getwd()
write.table(ps1ra_melt,"ps1ra_melt.tsv",row.names = FALSE,sep = "\t")

# plot_bar

#1. 
ggplot(ps1ra_melt,  aes(x=Rejection_1YA,  y=Abundance,  fill=Phylum,  order=as.factor(Phylum)))+ geom_bar(stat  =  "identity",  color  =  "black", position = "stack")+
  ggtitle("Sample  -  Relative  Abuncance,  sorted  by  Phylum  abundance")  + theme(#legend.position = "none",
    axis.text.x  =  element_blank())  

#2.

plot_bar(ps_relabund_0m, fill="Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat="identity", position="stack") +
  labs(x = "", y = "Relative Abundance\n") +
  theme(panel.background = element_blank())

#3.
top20 <- names(sort(taxa_sums(ps_0m_filter), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps_0m_filter, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Rejection_1YA", fill="Genus") 


#5. How many genera would be present after filtering?
length(get_taxa_unique(ps_0m_filter, taxonomic.rank = "Genus")) # 223
length(get_taxa_unique(ps_0m_rare, taxonomic.rank = "Genus")) # 258


#6. Prevalence
library("microbiome")
prevalence(ps_relabund_0m, detection = 0, sort=TRUE) %>% kable()

saveRDS(ps_relabund_0m, "inputs_0m/ps_relabund_0m.RDS") # No rarefying relative abundance table
saveRDS(ps_0m_filter, "inputs_0m/ps_0m_filter.RDS") # No rarefying otu table
saveRDS(ps_relabund_0m_rare, "inputs_0m/ps_relabund_0m_rare.RDS") # rarefying relative abundance table

# Lastly, we can use the plot_richness function to calculate ??-diversity (Shannon, Chao) for each sample and then save the values for subsequent plotting and statistical analyses

pAlpha = plot_richness(ps_0m_filter, measures = c("Observed", "Shannon"))
alphadt = data.table(pAlpha$data)
write.table(alphadt, "inputs_0m/KT_alpha_div.txt", row.names=F, sep="\t") # No rarefying diversity

pAlpha_rare = plot_richness(ps_0m_rare, measures = c("Observed", "Shannon")) 
alphadt_rare = data.table(pAlpha_rare$data)
write.table(alphadt_rare, "inputs_0m/KT_alpha_div_rare.txt", row.names=F, sep="\t") # rarefying alpha diversity




### ???????? ??????..

# Fix month levels in sample_data
sample_variables(ps_0m_rare)
str(ps_0m_rare)

# renames variable
sample_data(ps_0m_rare)$Rejection_1YA <- factor(sample_data(ps_0m_rare)$Rejection_1YA, levels = list("Rejection" = "Rejection", "No Rejection" = "No Rejection", "No Rejection"="No rejection"))

sample_data(ps_0m_rare)$Rejection_1YA <- factor(
  sample_data(ps_0m_rare)$Rejection_1YA, 
  levels = c("No Rejection","Rejection"))

sample_data(ps_0m_rare)$Rejection_1YA_B <- factor(
  sample_data(ps_0m_rare)$Rejection_1YA_B, 
  levels = c("No Rejection","Rejection"))

sample_data(ps_0m_rare)$Rejection_1YA_all.B <- factor(
  sample_data(ps_0m_rare)$Rejection_1YA_all.B, 
  levels = c("No Rejection","Rejection"))

# alpha-diversity
rich_0m_rare<-estimate_richness(ps_0m_rare)
rich_0m_rare

rich_0m<-estimate_richness(ps_0m_filter)
rich_0m

 ## Rejection

(p <- plot_richness(ps_0m_rare, x = "Rejection_1YA", color = "Rejection_1YA", measures = c("Observed", "Shannon")))
p+geom_boxplot()
p

wilcox.test(rich_0m$Observed~sample_data(ps_0m_rare)$Rejection_1YA, p.adjust.method="fdr")
wilcox.test(rich_0m$Shannon~sample_data(ps_0m_rare)$Rejection_1YA,p.adjust.method="fdr")


# Violin plot of alph diversity
sample_colors <- c("No Rejection"="#2d51b5","Rejection"="#c70643")
sample_types <- c("No Rejection","Rejection")
sample_labels <- c("No Rejection"="No Rejection","Rejection"="Rejection")


p<-plot_richness(ps_0m_rare,
                 x = "Rejection_1YA",   #compare diversity of datatype
                 measures = c("Observed", "Shannon")) +  #choose diversity measures
  geom_violin(aes(fill = "Rejection_1YA"), show.legend = FALSE)+ #make violin plot, set fill aes to Time
  geom_boxplot(width=0.1) + #add boxplot, set width
  theme_classic()+       #change theme to classic
  xlab(NULL)+    #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),  #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 15))+  #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+ #adjust headings
  scale_fill_manual(values= sample_colors) +   #set fill colors
  scale_x_discrete(          #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +        #add title
  theme(plot.title=element_text(size = 15, face = "bold", hjust = 0.5)) #change title size, face and position
p




## G_Age

sample_data(ps_0m_rare)$G_Age <- factor(
  sample_data(ps_0m_rare)$G_Age, 
  levels = c("Twenty", "Thirty", "Forty","Fifty","Sixty"))

(p <- plot_richness(ps_0m_rare, x = "G_Age", color = "G_Age", measures = c("Observed", "Shannon")))
p+geom_boxplot()
p

kruskal.test(rich_0m$Observed~sample_data(ps_0m_rare)$G_Age)
pairwise.wilcox.test(rich_0m$Observed,sample_data(ps_0m_rare)$G_Age,p.adjust.method="fdr")

kruskal.test(rich_0m$Shannon~sample_data(ps_0m_rare)$G_Age)
pairwise.wilcox.test(rich_0m$Shannon,sample_data(ps_0m_rare)$G_Age,p.adjust.method="fdr")

sample_colors <- c("No rejection"="#2874C5", "Rejection"="#EABF00")
sample_types <- c("Twenty", "Thirty", "Forty","Fifty","Sixty")
sample_labels <- c("Twenty"="Twenty", "Thirty"="Thirty", "Forty"="Forty","Fifty"="Fifty","Sixty"="Sixty")

p<-plot_richness(ps_0m_rare,
                 x = "G_Age",   #compare diversity of datatype
                 measures = c("Observed", "Shannon")) +  #choose diversity measures
  geom_violin(aes(fill = "G_Age"), show.legend = FALSE)+ #make violin plot, set fill aes to Time
  geom_boxplot(width=0.1) + #add boxplot, set width
  theme_classic()+       #change theme to classic
  xlab(NULL)+    #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),  #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 15))+  #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+ #adjust headings
  scale_fill_manual(values= sample_colors) +   #set fill colors
  scale_x_discrete(          #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +        #add title
  theme(plot.title=element_text(size = 15, face = "bold", hjust = 0.5)) #change title size, face and position
p


# =========



## G_BMI

sample_data(ps_0m_rare)$G_BMI <- factor(
  sample_data(ps_0m_rare)$G_BMI, 
  levels = c("Underweight", "Normal", "Overweight","Obese"))

(p <- plot_richness(ps_0m_rare, x = "G_BMI", color = "G_BMI", measures = c("Observed", "Shannon")))
p+geom_boxplot()
p

kruskal.test(rich_0m$Observed~sample_data(ps_0m_rare)$G_BMI)
pairwise.wilcox.test(rich_0m$Observed,sample_data(ps_0m_rare)$G_BMI,p.adjust.method="fdr")

kruskal.test(rich_0m$Shannon~sample_data(ps_0m_rare)$G_BMI)
pairwise.wilcox.test(rich_0m$Shannon,sample_data(ps_0m_rare)$G_BMI,p.adjust.method="fdr")

sample_colors <- c("No rejection"="#2874C5", "Rejection"="#EABF00")
sample_types <- c("Underweight", "Normal", "Overweight","Obese")
sample_labels <- c("Underweight"="Underweight", "Thirty"="Thirty", "Forty"="Forty","Fifty"="Fifty","Sixty"="Sixty")

p<-plot_richness(ps_0m_rare,
                 x = "G_BMI",   #compare diversity of datatype
                 measures = c("Observed", "Shannon")) +  #choose diversity measures
  geom_violin(aes(fill = "G_BMI"), show.legend = FALSE)+ #make violin plot, set fill aes to Time
  geom_boxplot(width=0.1) + #add boxplot, set width
  theme_classic()+       #change theme to classic
  xlab(NULL)+    #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),  #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 15))+  #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+ #adjust headings
  scale_fill_manual(values= sample_colors) +   #set fill colors
  scale_x_discrete(          #change x-axis labels
    breaks = sample_types,
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +        #add title
  theme(plot.title=element_text(size = 15, face = "bold", hjust = 0.5)) #change title size, face and position
p

## ISA_before

sample_data(ps_0m_rare)$ISA_before <- factor(
  sample_data(ps_0m_rare)$ISA_before, 
  levels = c("No ISA", "ISA"))

(p <- plot_richness(ps_0m_rare, x = "ISA_before", color = "ISA_before", measures = c("Observed", "Shannon")))
p+geom_boxplot()
p


wilcox.test(rich_0m$Observed~sample_data(ps_0m_rare)$ISA_before, p.adjust.method="fdr")
wilcox.test(rich_0m$Shannon~sample_data(ps_0m_rare)$ISA_before,p.adjust.method="fdr")


sample_colors <- c("No rejection"="#2874C5", "Rejection"="#EABF00")
sample_types <- c("No ISA", "ISA")
sample_labels <- c("No ISA"="No ISA", "ISA"="ISA")

p<-plot_richness(ps_0m_rare,
                 x = "ISA_before",   #compare diversity of datatype
                 measures = c("Observed", "Shannon")) +  #choose diversity measures
  geom_violin(aes(fill = "ISA_before"), show.legend = FALSE)+ #make violin plot, set fill aes to Time
  geom_boxplot(width=0.1) + #add boxplot, set width
  theme_classic()+       #change theme to classic
  xlab(NULL)+    #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),  #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 15))+  #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+ #adjust headings
  scale_fill_manual(values= sample_colors) +   #set fill colors
  scale_x_discrete(          #change x-axis labels
    breaks = sample_types,
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +        #add title
  theme(plot.title=element_text(size = 15, face = "bold", hjust = 0.5)) #change title size, face and position
p

## Abts_before

sample_data(ps_0m_rare)$Abts_before <- factor(
  sample_data(ps_0m_rare)$Abts_before, 
  levels = c("No Abts", "Abts"))

(p <- plot_richness(ps_0m_rare, x = "Abts_before", color = "Abts_before", measures = c("Observed", "Shannon")))
p+geom_boxplot()
p


wilcox.test(rich_0m$Observed~sample_data(ps_0m_rare)$Abts_before, p.adjust.method="fdr")
wilcox.test(rich_0m$Shannon~sample_data(ps_0m_rare)$Abts_before,p.adjust.method="fdr")


sample_colors <- c("No rejection"="#2874C5", "Rejection"="#EABF00")
sample_types <- c("No Abts", "Abts")
sample_labels <- c("No Abts"="No Abts", "Abts"="Abts")

p<-plot_richness(ps_0m_rare,
                 x = "Abts_before",   #compare diversity of datatype
                 measures = c("Observed", "Shannon")) +  #choose diversity measures
  geom_violin(aes(fill = "Abts_before"), show.legend = FALSE)+ #make violin plot, set fill aes to Time
  geom_boxplot(width=0.1) + #add boxplot, set width
  theme_classic()+       #change theme to classic
  xlab(NULL)+    #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),  #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 15))+  #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+ #adjust headings
  scale_fill_manual(values= sample_colors) +   #set fill colors
  scale_x_discrete(          #change x-axis labels
    breaks = sample_types,
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +        #add title
  theme(plot.title=element_text(size = 15, face = "bold", hjust = 0.5)) #change title size, face and position
p

## Antiviral

sample_data(ps_0m_rare)$Antiviral <- factor(
  sample_data(ps_0m_rare)$Antiviral, 
  levels = c("No Antiviral", "Antiviral"))

(p <- plot_richness(ps_0m_rare, x = "Antiviral", color = "Antiviral", measures = c("Observed", "Shannon")))
p+geom_boxplot()
p


wilcox.test(rich_0m$Observed~sample_data(ps_0m_rare)$Antiviral, p.adjust.method="fdr")
wilcox.test(rich_0m$Shannon~sample_data(ps_0m_rare)$Antiviral,p.adjust.method="fdr")


sample_colors <- c("No rejection"="#2874C5", "Rejection"="#EABF00")
sample_types <- c("No Antiviral", "Antiviral")
sample_labels <- c("No Antiviral"="No Antiviral", "Antiviral"="Antiviral")

p<-plot_richness(ps_0m_rare,
                 x = "Antiviral",   #compare diversity of datatype
                 measures = c("Observed", "Shannon")) +  #choose diversity measures
  geom_violin(aes(fill = "Antiviral"), show.legend = FALSE)+ #make violin plot, set fill aes to Time
  geom_boxplot(width=0.1) + #add boxplot, set width
  theme_classic()+       #change theme to classic
  xlab(NULL)+    #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),  #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 15))+  #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+ #adjust headings
  scale_fill_manual(values= sample_colors) +   #set fill colors
  scale_x_discrete(          #change x-axis labels
    breaks = sample_types,
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +        #add title
  theme(plot.title=element_text(size = 15, face = "bold", hjust = 0.5)) #change title size, face and position
p


## Induction

sample_data(ps_0m_rare)$Induction <- factor(
  sample_data(ps_0m_rare)$Induction, 
  levels = c("No", "Basiliximab", "ATG","Both"))

(p <- plot_richness(ps_0m_rare, x = "Induction", color = "Induction", measures = c("Observed", "Shannon")))
p+geom_boxplot()
p

kruskal.test(rich_0m$Observed~sample_data(ps_0m_rare)$Induction)
pairwise.wilcox.test(rich_0m$Observed,sample_data(ps_0m_rare)$Induction,p.adjust.method="fdr")

kruskal.test(rich_0m$Shannon~sample_data(ps_0m_rare)$Induction)
pairwise.wilcox.test(rich_0m$Shannon,sample_data(ps_0m_rare)$Induction,p.adjust.method="fdr")



sample_colors <- c("No rejection"="#2874C5", "Rejection"="#EABF00")
sample_types <- c("No", "Basiliximab", "ATG","Both")
sample_labels <- c("No"="No", "Basiliximab"="Basiliximab", "ATG"="ATG","Both"="Both")

p<-plot_richness(ps_0m_rare,
                 x = "Induction",   #compare diversity of datatype
                 measures = c("Observed", "Shannon")) +  #choose diversity measures
  geom_violin(aes(fill = "Induction"), show.legend = FALSE)+ #make violin plot, set fill aes to Time
  geom_boxplot(width=0.1) + #add boxplot, set width
  theme_classic()+       #change theme to classic
  xlab(NULL)+    #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),  #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 15))+  #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+ #adjust headings
  scale_fill_manual(values= sample_colors) +   #set fill colors
  scale_x_discrete(          #change x-axis labels
    breaks = sample_types,
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +        #add title
  theme(plot.title=element_text(size = 15, face = "bold", hjust = 0.5)) #change title size, face and position
p


## Sex


(p <- plot_richness(ps_0m_rare, x = "sex", color = "sex", measures = c("Observed", "Shannon")))
p+geom_boxplot()
p

wilcox.test(rich_0m$Observed~sample_data(ps_0m_rare)$sex, p.adjust.method="fdr")
wilcox.test(rich_0m$Shannon~sample_data(ps_0m_rare)$sex,p.adjust.method="fdr")

sample_colors <- c("Female"="#2874C5", "Male"="#EABF00")
sample_types <- c("Female","Male")
sample_labels <- c("Female"="Female","Male"="Male")

p<-plot_richness(ps_0m_rare,
                 x = "Sex",   #compare diversity of datatype
                 measures = c("Observed", "Shannon")) +  #choose diversity measures
  geom_violin(aes(fill = "Sex"), show.legend = FALSE)+ #make violin plot, set fill aes to Time
  geom_boxplot(width=0.1) + #add boxplot, set width
  theme_classic()+       #change theme to classic
  xlab(NULL)+    #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),  #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 15))+  #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+ #adjust headings
  scale_fill_manual(values= sample_colors) +   #set fill colors
  scale_x_discrete(          #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +        #add title
  theme(plot.title=element_text(size = 15, face = "bold", hjust = 0.5)) #change title size, face and position
p



## Donor


(p <- plot_richness(ps_0m_rare, x = "Donor", color = "Donor", measures = c("Observed", "Shannon")))
p+geom_boxplot()
p

wilcox.test(rich_0m$Observed~sample_data(ps_0m_rare)$Donor, p.adjust.method="fdr")
wilcox.test(rich_0m$Shannon~sample_data(ps_0m_rare)$Donor,p.adjust.method="fdr")

sample_colors <- c("Female"="#2874C5", "Male"="#EABF00")
sample_types <- c("DDKT","LDKT")
sample_labels <- c("DDKT"="DDKT","LDKT"="LDKT")

p<-plot_richness(ps_0m_rare,
                 x = "Donor",   #compare diversity of datatype
                 measures = c("Observed", "Shannon")) +  #choose diversity measures
  geom_violin(aes(fill = "Donor"), show.legend = FALSE)+ #make violin plot, set fill aes to Time
  geom_boxplot(width=0.1) + #add boxplot, set width
  theme_classic()+       #change theme to classic
  xlab(NULL)+    #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),  #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 15))+  #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+ #adjust headings
  scale_fill_manual(values= sample_colors) +   #set fill colors
  scale_x_discrete(          #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +        #add title
  theme(plot.title=element_text(size = 15, face = "bold", hjust = 0.5)) #change title size, face and position
p

## DM=======


(p <- plot_richness(ps_0m_rare, x = "DM", color = "DM", measures = c("Observed", "Shannon")))
p+geom_boxplot()
p

wilcox.test(rich_0m$Observed~sample_data(ps_0m_rare)$DM, p.adjust.method="fdr")
wilcox.test(rich_0m$Shannon~sample_data(ps_0m_rare)$DM,p.adjust.method="fdr")

sample_colors <- c("No DM"="#2874C5", "DM"="#EABF00")
sample_types <- c("No DM","DM")
sample_labels <- c("No DM"="No DM","DM"="DM")

p<-plot_richness(ps_0m_rare,
                 x = "DM",   #compare diversity of datatype
                 measures = c("Observed", "Shannon")) +  #choose diversity measures
  geom_violin(aes(fill = "DM"), show.legend = FALSE)+ #make violin plot, set fill aes to Time
  geom_boxplot(width=0.1) + #add boxplot, set width
  theme_classic()+       #change theme to classic
  xlab(NULL)+    #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),  #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 15))+  #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+ #adjust headings
  scale_fill_manual(values= sample_colors) +   #set fill colors
  scale_x_discrete(          #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +        #add title
  theme(plot.title=element_text(size = 15, face = "bold", hjust = 0.5)) #change title size, face and position
p

# Desentization=======
sample_data(ps_0m_rare)$desensitization<- factor(
  sample_data(ps_0m_rare)$desensitization, 
  levels = c("No Desensitization", "Desensitization"))


(p <- plot_richness(ps_0m_rare, x = "desensitization", color = "desensitization", measures = c("Observed", "Shannon")))
p+geom_boxplot()
p

wilcox.test(rich_0m$Observed~sample_data(ps_0m_rare)$desensitization, p.adjust.method="fdr")
wilcox.test(rich_0m$Shannon~sample_data(ps_0m_rare)$desensitization,p.adjust.method="fdr")

sample_colors <- c("0"="#2874C5", "1"="#EABF00")
sample_types <- c("No Desensitization", "Desensitization")
sample_labels <- c("No Desensitization"="No Desensitization","Desensitization"="Desensitization")

p<-plot_richness(ps_0m_rare,
                 x = "desensitization",   #compare diversity of datatype
                 measures = c("Observed", "Shannon")) +  #choose diversity measures
  geom_violin(aes(fill = "desensitization"), show.legend = FALSE)+ #make violin plot, set fill aes to Time
  geom_boxplot(width=0.1) + #add boxplot, set width
  theme_classic()+       #change theme to classic
  xlab(NULL)+    #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),  #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 15))+  #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+ #adjust headings
  scale_fill_manual(values= sample_colors) +   #set fill colors
  scale_x_discrete(          #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +        #add title
  theme(plot.title=element_text(size = 15, face = "bold", hjust = 0.5)) #change title size, face and position
p


# Modality=======

sample_data(ps_0m_rare)$Modality<- factor(
  sample_data(ps_0m_rare)$Modality, 
  levels = c("Preemptive", "HD","PD","HDPD"))


(p <- plot_richness(ps_0m_rare, x = "Modality", color = "Modality", measures = c("Observed", "Shannon")))
p+geom_boxplot()
p

kruskal.test(rich_0m$Observed~sample_data(ps_0m_rare)$Modality)
pairwise.wilcox.test(rich_0m$Observed,sample_data(ps_0m_rare)$Modality,p.adjust.method="fdr")

kruskal.test(rich_0m$Shannon~sample_data(ps_0m_rare)$Modality)
pairwise.wilcox.test(rich_0m$Shannon,sample_data(ps_0m_rare)$Modality,p.adjust.method="fdr")


sample_colors <- c("0"="#2874C5", "1"="#EABF00")
sample_types <- c("Preemptive", "HD","PD","HDPD")
sample_labels <- c("Preemptive"="Preemptive","HD"="HD","PD"="PD","HDPD"="HD+PD")

p<-plot_richness(ps_0m_rare,
                 x = "Modality",   #compare diversity of datatype
                 measures = c("Observed", "Shannon")) +  #choose diversity measures
  geom_violin(aes(fill = "Modality"), show.legend = FALSE)+ #make violin plot, set fill aes to Time
  geom_boxplot(width=0.1) + #add boxplot, set width
  theme_classic()+       #change theme to classic
  xlab(NULL)+    #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),  #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 15))+  #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+ #adjust headings
  scale_fill_manual(values= sample_colors) +   #set fill colors
  scale_x_discrete(          #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +        #add title
  theme(plot.title=element_text(size = 15, face = "bold", hjust = 0.5)) #change title size, face and position
p

# Abts_ativiral====

sample_data(ps_0m_rare)$Abts_ativiral<- factor(
  sample_data(ps_0m_rare)$Abts_ativiral, 
  levels = c("0", "1"))


(p <- plot_richness(ps_0m_rare, x = "Abts_ativiral", color = "Abts_ativiral", measures = c("Observed", "Shannon")))
p+geom_boxplot()
p

wilcox.test(rich_0m$Observed~sample_data(ps_0m_rare)$Abts_ativiral, p.adjust.method="fdr")
wilcox.test(rich_0m$Shannon~sample_data(ps_0m_rare)$Abts_ativiral,p.adjust.method="fdr")

kruskal.test(rich_0m$Observed~sample_data(ps_0m_rare)$Abts_ativiral)
pairwise.wilcox.test(rich_0m$Observed,sample_data(ps_0m_rare)$Abts_ativiral,p.adjust.method="fdr")

kruskal.test(rich_0m$Shannon~sample_data(ps_0m_rare)$Abts_ativiral)
pairwise.wilcox.test(rich_0m$Shannon,sample_data(ps_0m_rare)$Abts_ativiral,p.adjust.method="fdr")

sample_colors <- c("0"="#2874C5", "1"="#EABF00")
sample_types <- c("0","1")
sample_labels <- c("0"="No Abts + Antiviral","1"="Abts + Antiviral")

p<-plot_richness(ps_0m_rare,
                 x = "Abts_ativiral",   #compare diversity of datatype
                 measures = c("Observed", "Shannon")) +  #choose diversity measures
  geom_violin(aes(fill = "Abts_ativiral"), show.legend = FALSE)+ #make violin plot, set fill aes to Time
  geom_boxplot(width=0.1) + #add boxplot, set width
  theme_classic()+       #change theme to classic
  xlab(NULL)+    #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),  #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 15))+  #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+ #adjust headings
  scale_fill_manual(values= sample_colors) +   #set fill colors
  scale_x_discrete(          #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +        #add title
  theme(plot.title=element_text(size = 15, face = "bold", hjust = 0.5)) #change title size, face and position
p


# Rejection_1YA

(p <- plot_richness(ps_0m_rare, x = "Rejection_1YA", color = "Rejection_1YA", measures = c("Observed", "Shannon")))
p+geom_boxplot()
p

wilcox.test(rich_0m$Observed~sample_data(ps_0m_rare)$Rejection_1YA, p.adjust.method="fdr")
wilcox.test(rich_0m$Shannon~sample_data(ps_0m_rare)$Rejection_1YA,p.adjust.method="fdr")

# Violin plot of alph diversity
sample_colors <- c("No rejection"="#2874C5", "Rejection"="#EABF00")
sample_types <- c("No rejection","Rejection")
sample_labels <- c("No rejection"="No rejection","Rejection"="Rejection")


p<-plot_richness(ps_0m_rare,
                 x = "Rejection_1YA",   #compare diversity of datatype
                 measures = c("Observed", "Shannon")) +  #choose diversity measures
  geom_violin(aes(fill = "Rejection_1YA"), show.legend = FALSE)+ #make violin plot, set fill aes to Time
  geom_boxplot(width=0.1) + #add boxplot, set width
  theme_classic()+       #change theme to classic
  xlab(NULL)+    #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),  #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 15))+  #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+ #adjust headings
  scale_fill_manual(values= sample_colors) +   #set fill colors
  scale_x_discrete(          #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +        #add title
  theme(plot.title=element_text(size = 15, face = "bold", hjust = 0.5)) #change title size, face and position
p




# beta diversity: ???PCoA???on bray, plot type 2
#ordination
pcoa_0m <- ordinate(
  physeq = ps_0m_rare_new, 
  method = "PCoA", 
  distance = "bray")

# G_Age

sample_colors <- c("1"="coral","2"="cadetblue1","3"="darkviolet","4"="firebrick","5"="goldenrod")
sample_types <- c("1","2","3","4","5")
sample_labels <- c("1"="20's","2"="30's","3"="40's","4"="50's","5"="60's-70's")

p<-plot_ordination(
  physeq = ps_0m_rare,           #phyloseq object
  ordination = pcoa_0m)+     #ordination
  geom_point(aes(fill = G_Age, shape = G_Age), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22, 23,24,25))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

# Calculate bray curtis distance matrix: permanova
set.seed(1)

# Calculate bray curtis distance matrix
ps_bray_0m<- phyloseq::distance(ps_0m_rare, method = "bray")

# Adonis test
adonis(ps_bray_0m ~ sample_data(ps_0m_rare)$G_Age) # p 0.765

# Homogeneity of dispersion test
beta <- betadisper(ps_bray_0m, sample_data(ps_0m_rare)$G_Age)
permutest(beta) # p 0.016

 ## G_BMI

sample_colors <- c("Underweight"="coral","Normal"="cadetblue1","Overweight"="darkviolet","Obese"="firebrick")
sample_types <- c("Underweight", "Normal", "Overweight","Obese")


p<-plot_ordination(
  physeq = ps_0m_rare,           #phyloseq object
  ordination = pcoa_0m)+     #ordination
  geom_point(aes(fill = G_BMI, shape = G_BMI), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22, 23,24,25))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

# Calculate bray curtis distance matrix: permanova
set.seed(1)

# Calculate bray curtis distance matrix
ps_bray_0m<- phyloseq::distance(ps_0m_rare, method = "bray")

# Adonis test
adonis(ps_bray_0m ~ sample_data(ps_0m_rare)$G_BMI) # p 0.471

# Homogeneity of dispersion test
beta <- betadisper(ps_bray_0m, sample_data(ps_0m_rare)$G_BMI)
permutest(beta) # p 0.111

### Sex

sample_colors <- c("Female"="goldenrod", "Male"="darkviolet")
sample_types <- c("Female","Male")
sample_labels <- c("Female"="Female","Male"="Male")

p<-plot_ordination(
  physeq = ps_0m_rare,           #phyloseq object
  ordination = pcoa_0m)+     #ordination
  geom_point(aes(fill = Sex, shape = Sex), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

# Calculate bray curtis distance matrix: permanova
set.seed(1)

# Calculate bray curtis distance matrix
ps_bray_0m<- phyloseq::distance(ps_0m_rare, method = "bray")

# Adonis test
adonis(ps_bray_0m ~ sample_data(ps_0m_rare)$Sex) # p 0.029

# Homogeneity of dispersion test
beta <- betadisper(ps_bray_0m, sample_data(ps_0m_rare)$Sex)
permutest(beta) # p 0.2

### ISA

sample_colors <- c("No ISA"="goldenrod", "ISA"="darkviolet")
sample_types <- c("Female","Male")


p<-plot_ordination(
  physeq = ps_0m_rare,           #phyloseq object
  ordination = pcoa_0m)+     #ordination
  geom_point(aes(fill = ISA_before, shape = ISA_before), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

# Calculate bray curtis distance matrix: permanova
set.seed(1)

# Calculate bray curtis distance matrix
ps_bray_0m<- phyloseq::distance(ps_0m_rare, method = "bray")

# Adonis test
adonis(ps_bray_0m ~ sample_data(ps_0m_rare)$ISA_before) # p 0.679

# Homogeneity of dispersion test
beta <- betadisper(ps_bray_0m, sample_data(ps_0m_rare)$ISA_before)
permutest(beta) # p 0.927

### Abts

sample_colors <- c("No Abts"="goldenrod", "Abts"="darkviolet")


p<-plot_ordination(
  physeq = ps_0m_rare,           #phyloseq object
  ordination = pcoa_0m)+     #ordination
  geom_point(aes(fill = Abts_before, shape = Abts_before), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

# Calculate bray curtis distance matrix: permanova
set.seed(1)

# Calculate bray curtis distance matrix
ps_bray_0m<- phyloseq::distance(ps_0m_rare, method = "bray")

# Adonis test
adonis(ps_bray_0m ~ sample_data(ps_0m_rare)$Abts_before) # p 0.435

# Homogeneity of dispersion test
beta <- betadisper(ps_bray_0m, sample_data(ps_0m_rare)$Abts_before)
permutest(beta) # p 0.546

### Antiviral

sample_colors <- c("No Antiviral"="goldenrod", "Antiviral"="darkviolet")


p<-plot_ordination(
  physeq = ps_0m_rare,           #phyloseq object
  ordination = pcoa_0m)+     #ordination
  geom_point(aes(fill = Antiviral, shape = Antiviral), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

# Calculate bray curtis distance matrix: permanova
set.seed(1)

# Calculate bray curtis distance matrix
ps_bray_0m<- phyloseq::distance(ps_0m_rare, method = "bray")

# Adonis test
adonis(ps_bray_0m ~ sample_data(ps_0m_rare)$Antiviral) # p 0.35

# Homogeneity of dispersion test
beta <- betadisper(ps_bray_0m, sample_data(ps_0m_rare)$Antiviral)
permutest(beta) # p 0.001***


### Kt type

sample_colors <- c("DDKT"="goldenrod", "LDKT"="darkviolet")


p<-plot_ordination(
  physeq = ps_0m_rare,           #phyloseq object
  ordination = pcoa_0m)+     #ordination
  geom_point(aes(fill = Donor, shape = Donor), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

# Calculate bray curtis distance matrix: permanova
set.seed(1)

# Calculate bray curtis distance matrix
ps_bray_0m<- phyloseq::distance(ps_0m_rare, method = "bray")

# Adonis test
adonis(ps_bray_0m ~ sample_data(ps_0m_rare)$Donor) # p 0.184

# Homogeneity of dispersion test
beta <- betadisper(ps_bray_0m, sample_data(ps_0m_rare)$Donor)
permutest(beta) # p 0.062


### Modality

sample_colors <- c("Preemptive"="coral","HD"="cadetblue1","PD"="darkviolet","HDPD"="firebrick")

p<-plot_ordination(
  physeq = ps_0m_rare,           #phyloseq object
  ordination = pcoa_0m)+     #ordination
  geom_point(aes(fill = Modality, shape = Modality), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22,23,24))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

# Calculate bray curtis distance matrix: permanova
set.seed(1)

# Calculate bray curtis distance matrix
ps_bray_0m<- phyloseq::distance(ps_0m_rare, method = "bray")

# Adonis test
adonis(ps_bray_0m ~ sample_data(ps_0m_rare)$Modality) # p 0.187

# Homogeneity of dispersion test
beta <- betadisper(ps_bray_0m, sample_data(ps_0m_rare)$Modality)
permutest(beta) # p 0.001

### Induction

sample_colors <- c("No"="coral","Basiliximab"="cadetblue1","ATG"="darkviolet","Both"="firebrick")


p<-plot_ordination(
  physeq = ps_0m_rare,           #phyloseq object
  ordination = pcoa_0m)+     #ordination
  geom_point(aes(fill = Induction, shape = Induction), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22,23,24))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

# Calculate bray curtis distance matrix: permanova
set.seed(1)

# Calculate bray curtis distance matrix
ps_bray_0m<- phyloseq::distance(ps_0m_rare, method = "bray")

# Adonis test
adonis(ps_bray_0m ~ sample_data(ps_0m_rare)$Induction) # p 0.708

# Homogeneity of dispersion test
beta <- betadisper(ps_bray_0m, sample_data(ps_0m_rare)$Induction)
permutest(beta) # p 0.001


### DM

sample_colors <- c("No DM"="goldenrod", "DM"="darkviolet")
sample_types <- c("No DM","DM")
sample_labels <- c("No DM"="No DM","DM"="DM")

p<-plot_ordination(
  physeq = ps_0m_rare,           #phyloseq object
  ordination = pcoa_0m)+     #ordination
  geom_point(aes(fill = DM, shape = DM), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p


### Desentiazation

sample_colors <- c("0"="goldenrod", "1"="darkviolet")
sample_types <- c("0","1")
sample_labels <- c("0"="No Desensitization","1"="Desensitization")


p<-plot_ordination(
  physeq = ps_0m_rare,           #phyloseq object
  ordination = pcoa_0m)+     #ordination
  geom_point(aes(fill = Desensitization, shape = Desensitization), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

# Calculate bray curtis distance matrix: permanova
set.seed(1)

# Calculate bray curtis distance matrix
ps_bray_0m<- phyloseq::distance(ps_0m_rare, method = "bray")

# Adonis test
adonis(ps_bray_0m ~ sample_data(ps_0m_rare)$Desensitization) # p 0.019

# Homogeneity of dispersion test
beta <- betadisper(ps_bray_0m, sample_data(ps_0m_rare)$Desensitization)
permutest(beta) # p 0.531


# Abts_ativiral

sample_colors <- c("0"="goldenrod", "1"="darkviolet")
sample_types <- c("0","1")
sample_labels <- c("0"="No Abts + Antiviral","1"="Abts + Antiviral")


p<-plot_ordination(
  physeq = ps_0m_rare,           #phyloseq object
  ordination = pcoa_0m)+     #ordination
  geom_point(aes(fill = Abts_ativiral, shape = Abts_ativiral), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

# Calculate bray curtis distance matrix: permanova
set.seed(1)

# Calculate bray curtis distance matrix
ps_bray_0m<- phyloseq::distance(ps_0m_rare, method = "bray")

# Adonis test
adonis(ps_bray_0m ~ sample_data(ps_0m_rare)$Abts_ativiral) # ?

# Homogeneity of dispersion test
beta <- betadisper(ps_bray_0m, sample_data(ps_0m_rare)$Abts_ativiral)
permutest(beta) # p ?

# Rejection_1YA

sample_colors <- c("No Rejection"="goldenrod", "Rejection"="darkviolet")


p<-plot_ordination(
  physeq = ps_0m_rare,           #phyloseq object
  ordination = pcoa_0m)+     #ordination
  geom_point(aes(fill = Rejection_1YA, shape = Rejection_1YA), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

# Calculate bray curtis distance matrix: permanova
set.seed(1)

# Calculate bray curtis distance matrix
ps_bray_0m<- phyloseq::distance(ps_0m_rare, method = "bray")

# Adonis test
adonis(ps_bray_0m ~ sample_data(ps_0m_rare)$Rejection_1YA) # p 0.005

# Homogeneity of dispersion test
beta <- betadisper(ps_bray_0m, sample_data(ps_0m_rare)$Rejection_1YA)
permutest(beta) # p 0.415

# Rejection_1YA_B

sample_colors <- c("No Rejection"="goldenrod", "Rejection"="darkviolet")


p<-plot_ordination(
  physeq = ps_0m_rare,           #phyloseq object
  ordination = pcoa_0m)+     #ordination
  geom_point(aes(fill = Rejection_1YA_B, shape = Rejection_1YA_B), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

# Calculate bray curtis distance matrix: permanova
set.seed(1)

# Calculate bray curtis distance matrix
ps_bray_0m<- phyloseq::distance(ps_0m_rare, method = "bray")

# Adonis test
adonis(ps_bray_0m ~ sample_data(ps_0m_rare)$Rejection_1YA_B) # p 0.042

# Homogeneity of dispersion test
beta <- betadisper(ps_bray_0m, sample_data(ps_0m_rare)$Rejection_1YA_B)
permutest(beta) # p 0.728

# Rejection_1YA_B

sample_colors <- c("No Rejection"="goldenrod", "Rejection"="darkviolet")


p<-plot_ordination(
  physeq = ps_0m_rare,           #phyloseq object
  ordination = pcoa_0m)+     #ordination
  geom_point(aes(fill = Rejection_1YA_all.B, shape = Rejection_1YA_all.B), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

# Calculate bray curtis distance matrix: permanova
set.seed(1)

# Calculate bray curtis distance matrix
ps_bray_0m<- phyloseq::distance(ps_0m_rare, method = "bray")

# Adonis test
adonis(ps_bray_0m ~ sample_data(ps_0m_rare)$Rejection_1YA_all.B) # p 0.157

# Homogeneity of dispersion test
beta <- betadisper(ps_bray_0m, sample_data(ps_0m_rare)$Rejection_1YA_all.B)
permutest(beta) # p 0.141



### Differential analysis
#1. Preprocessing

ps_0m_filter
ps_0m_genus<-subset_taxa(ps_0m_filter, !is.na(Genus)) # remove na in genus
ps_0m_genus

# examining the taxonomy
rank_names(ps_0m_genus)
head(tax_table(ps_0m_genus))
table(tax_table(ps_0m_genus)[,2])
tax_tab_0m<-data.frame(tax_table(ps_0m_genus))

# Aggregate to genus level
library(microbiome)

genus_data_0m = aggregate_taxa(ps_0m_genus, "Genus")
genus_data_0m


(genus_data_0m_10 <- filter_taxa(genus_data_0m, function(x) sum(x > 0) > (0.1*length(x)), TRUE))   #removing species not seen > 10% of samples

head(sort(sample_sums(genus_data_0m_10))) 

genus_name = unlist(taxa_names(genus_data_0m_10))

# Output data
library(readr)
library(tibble)

write_tsv(data.frame(abundances(genus_data_0m_10), check.names = FALSE) %>%
            rownames_to_column("genus"), "genus_data_0m_10.txt")

write_tsv(meta(genus_data_0m_10), path = "meta_0m.txt") 
tax_0m = data.frame(tax_table(genus_data_0m_10)@.Data, check.names = FALSE)

library(tidyr)
tax_0m = tax_0m %>% rownames_to_column("Feature ID") %>%
  unite(col = "Taxon", Kingdom:Genus, sep = ";") %>%
  dplyr::select(-unique)

write_tsv(tax_0m, path = "tax_0m.txt")

sample_data(genus_data_0m_10)$Rejection_1YA <- factor(sample_data(genus_data_0m_10)$Rejection_1YA, levels = c("No Rejection", "Rejection"))
table(sample_data(genus_data_0m_10)$Rejection_1YA) 

tax_mat_0m = as(tax_table(genus_data_0m_10), "matrix")


# plot bar-relative abundance

#1.abundance 0.01
genus_0m_relative<-transform_sample_counts(genus_data_0m_10 , function(OTU) OTU/sum(OTU) )

genus_0m_relative= filter_taxa(genus_0m_relative, function(x) mean(x) > 1e-2, TRUE)
# abundance 0.01 

ps1<-merge_samples(genus_0m_relative,"Rejection_1YA")
warnings()
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Genus")


#2. abundance 0.001
genus_0m_relative<-transform_sample_counts(genus_data_0m_10 , function(OTU) OTU/sum(OTU) )

genus_0m_relative1= filter_taxa(genus_0m_relative, function(x) mean(x) > 1e-3, TRUE)


ps1<-merge_samples(genus_0m_relative1,"Rejection_1YA")
ps2 <- transform_sample_counts(ps1, function(x) x / sum(x))
plot_bar(ps2, fill="Genus")

# Running MaAsLin2
library(Maaslin2)


# data trasnformation to CRL
library(microbiome)
?microbiome::transform
ps_relative_0m<-microbiome::transform(genus_data_0m_10,"compositional")
otu_relative_0m<-t(data.frame(otu_table(ps_relative_0m))) # relative otu table (******)

ps_clr_0m<-microbiome::transform(genus_data_0m_10,"clr")
otu_clr_0m<-t(data.frame(otu_table(ps_clr_0m))) # clr transformed otu table

metadata_0m<-data.frame(sample_data(genus_data_0m_10))

write.csv(metadata_0m,"metadata_0m.csv")
write.csv(otu_relative_0m,"otu_relative_0m.csv")

getwd()
input_metadata_0m<-read.csv("/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/metadata_0m.csv",header=TRUE,sep=",",row.names = 1)


input_data_r_0m <- data.frame(otu_relative_0m)
input_data_clr_0m<-data.frame(otu_clr_0m)

input_metadata_0m$Desensitization<-factor(input_metadata_0m$Desensitization)
input_metadata_0m$Desensitization=revalue(input_metadata_0m$Desensitization,replace=c("0"="No desensitization","1"="Desensitization"))

# Rejection_1YA

mas_0m <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/Maaslin2_output_r_0m_R1Y",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = "Rejection_1YA",
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/Maaslin2_output_r_0m_R1Y_G_Agesex",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/Maaslin2_output_r_0m_R1Y_G_AgesexG_BMI",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/Maaslin2_output_r_0m_R1Y_G_AgesexG_BMI_DM",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/Maaslin2_output_r_0m_R1Y_G_AgesexG_BMI_DM_Abts",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/Maaslin2_output_r_0m_R1Y_G_AgesexG_BMI_DM_Abts_ISA",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before","ISA_before"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts","ISA_before,No ISA"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/Maaslin2_output_r_0m_R1Y_G_AgesexG_BMI_DM_Abts_ISA_desensitization",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before","ISA_before","desensitization"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts","ISA_before,No ISA","desensitization, No Desensitization"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

=============================

mas_0m_2 <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210727/Maaslin2_output_r_0m_agesexBMIDM",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","Age","Sex","BMI","DM"),
  reference=c("Rejection_1YA,No rejection", "Sex,Male","DM,No DM"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_2 <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210727/Maaslin2_output_r_0m_agesexBMIDM_desensi",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","Age","Sex","BMI","DM","Desensitization"),
  reference=c("Rejection_1YA,No rejection", "Sex,Male","DM,No DM","Desentization, No desensitization"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)



# crosstable
library(moonBook)

table0<-mytable(Rejection_1YA~., data=input_metadata_0m) 
str(input_metadata_0m)

mycsv(table1,file="table1.csv")

table1<-mytable(Rejection_1YA~., data=input_metadata_0m)
mycsv(table1,file="table1.csv")

table3<-mytable(Rejection_1YA_all.B~., data=input_metadata_0m)
mycsv(table3,file="table3.csv")




# Rejection_1YA_B

mas_0m <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/Maaslin2_output_r_0m_R1YB",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = "Rejection_1YA_B",
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/Maaslin2_output_r_0m_R1Y_Borderlibe_G_Agesex",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA_B","G_Age","sex"),
  reference=c("Rejection_1YA_B,No Rejection", "G_Age,Twenty","sex,Male"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/Maaslin2_output_r_0m_R1Y_G_AgesexG_BMI",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/Maaslin2_output_r_0m_R1Y_G_AgesexG_BMI_DM",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/Maaslin2_output_r_0m_R1Y_G_AgesexG_BMI_DM_Abts",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/Maaslin2_output_r_0m_R1Y_G_AgesexG_BMI_DM_Abts_ISA",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before","ISA_before"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts","ISA_before,No ISA"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/Maaslin2_output_r_0m_R1Y_G_AgesexG_BMI_DM_Abts_ISA_desensitization",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before","ISA_before","desensitization"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts","ISA_before,No ISA","desensitization, No Desensitization"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

=============================
  
  mas_0m_2 <- Maaslin2(
    input_data_r_0m,
    input_metadata_0m,
    output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210727/Maaslin2_output_r_0m_agesexBMIDM",
    min_abundance = 0.001,
    min_prevalence = 0.1,
    fixed_effects = c("Rejection_1YA","Age","Sex","BMI","DM"),
    reference=c("Rejection_1YA,No rejection", "Sex,Male","DM,No DM"),
    correction = "BH",
    max_significance = 0.25,
    normalization = 'NONE',
    standardize = FALSE)

mas_0m_2 <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210727/Maaslin2_output_r_0m_agesexBMIDM_desensi",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","Age","Sex","BMI","DM","Desensitization"),
  reference=c("Rejection_1YA,No rejection", "Sex,Male","DM,No DM","Desentization, No desensitization"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)



# Rejection_1YA_allB

mas_0m <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/Maaslin2_output_r_0m_R1YallB",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = "Rejection_1YA_all.B",
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_data_r_0m,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime-total/Phyloseq_total_20210818/Maaslin2_output_r_0m_R1Y_allB_G_Agesex",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA_all.B","G_Age","sex"),
  reference=c("Rejection_1YA_all.B,No Rejection", "G_Age,Twenty","sex,Male"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)


# PICRUST2 #############
 ##  preprocessing data for aldex
getwd()

library(ALDEx2)
library(readr)
library(dplyr)

df=read_tsv("~/Dropbox/KTPL_2021/Qiime_202108/Picrust2/picrust2_out_pipeline/pathways_out/path_abun_strat.tsv/path_abun_strat.tsv")
dfs=as.data.frame(df)
rownames(dfs)=paste(df$pathway,df$sequence,sep="|")
dfs=dfs[,-which(colnames(dfs) %in% c("pathway","sequence"))]
dfout=as.data.frame(df)
dfus=read_tsv("~/Dropbox/KTPL_2021/Qiime_202108/Picrust2/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv")
dfus=as.data.frame(dfus) # stratified metacycle file
rownames(dfus)=dfus$pathway
dfus=dfus[,-which(colnames(dfus)=="pathway")] # unstratified metacycle file

#nothing in the pathway, look in the kegg?
dfk=read_tsv("~/Dropbox/KTPL_2021/Qiime_202108/Picrust2/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv")
dfk=as.data.frame(dfk)
rownames(dfk)=dfk$"function"
dfk=dfk[,-which(colnames(dfk)=="function")]

rownames(dfout)=paste(dfout$pathway,dfout$sequence,sep="|")
dfout=dfout[,-which(colnames(dfout) %in% c("pathway","sequence"))]

df_tmp= df[, -which(colnames(df) == "sequence")]
df_abun=aggregate(.~pathway,data=df_tmp,FUN=sum) # sum pathway
rownames(df_abun)=df_abun$pathway
df_abun=df_abun[,-which(colnames(df_abun)=="pathway")]

#convert to relative abundance
df_rabun=data.frame(sweep(df_abun,2,colSums(df_abun),"/"),check.names=F) # relative abundance table
# asin transformation
asinT=function(x) asin(sqrt(x))
df_rabun2=df_rabun %>% mutate_all(asinT) # asin transformation table

#finally, round them so they
df_rrabun=round(df_rabun) # round of relative abundance of stratified pathway

library(phyloseq)
library(microbiome)
getwd()
mdata=read_tsv("/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/metadata_0m.tsv")
kotu=otu_table(dfk,taxa_are_rows=TRUE)
md=as.data.frame(mdata)
rownames(md)=mdata$SampleID
md=md[,-which(colnames(md)=="SampleID")]
md=sample_data(md)
kegg=phyloseq(kotu,md)
pwo=otu_table(dfout,taxa_are_rows=TRUE)
pw=phyloseq(pwo,md)

pwo_us=otu_table(dfus, taxa_are_rows=TRUE)
pw_us=phyloseq(pwo_us,md)

#kegg is the kegg table
#pw is the pathway table (much higher level)
pw_clr=microbiome::transform(pw,"clr")
kegg_clr=microbiome::transform(kegg,"clr")

phyloseq::otu_table(pw_clr)[1:5, 1:6]
phyloseq::otu_table(kegg_clr)[1:5, 1:6]

#rda is one of many ordination methods
ord_pw=phyloseq::ordinate(pw_clr,"RDA")
ord_kegg=phyloseq::ordinate(kegg_clr,"RDA")
#get the top eigenvalues of the first few PC axes
head(ord_pw$CA$eig)

#transform to percent
sapply(ord_pw$CA$eig[1:8], function(x) x / sum(ord_pw$CA$eig))

sapply(ord_kegg$CA$eig[1:8], function(x) x / sum(ord_kegg$CA$eig))

?plot_ordination
p=plot_ordination(pw,ord_pw,type="samples",color="Rejection_1YA",shape="Rejection_1YA")
p=p+geom_polygon(aes(fill=Rejection_1YA))
p

# fancy-pants statistical modeling with Aldex2

conds=mdata$Rejection_1YA

#kegg
x.kf=aldex(round(dfk),conds,effects=T,denom="all")

#unstratified pathway
x.pwu=aldex(round(dfus),conds,effects=TRUE,denom="all")

#stratified pathway
x.pws=aldex(round(dfs),conds,effects=TRUE, denom="all")

# error:vector memory exhausted (limit reached?)
# setting:  Once open, the R_MAX_VSIZE var can be set. R_MAX_VSIZE=100Gb
install.packages("usethis")
usethis::edit_r_environ()

#x.all=aldex(round(dfs),conds,effects=TRUE)
#x.all <- aldex(selex.sub,conds,mc.samples=16,test="t",effect=TRUE,include.samples.summary=FALSE,demon="all",verbose=FALSE)

# aldex.plot
par(mfrow=(c(2,2)))
aldex.plot(x.kf,type="MA",test="welch",xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(x.kf,type="MW",test="welch",xlab="Dispersion",
           ylab="Difference")
aldex.plot(x.pwu,type="MA",test="welch",xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(x.pwu,type="MW",test="welch",xlab="Dispersion",
           ylab="Difference")
aldex.plot(x.pwu,type="MA",test="welch",xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(x.pwu,type="MW",test="welch",xlab="Dispersion",
           ylab="Difference")

# visualizing the most significant pathways
xsig=x.pwu[x.pwu$we.ep<0.05,]
xsig$path=rownames(xsig)
xsig2=arrange(xsig,desc(abs(effect)))
xsig3=xsig2[1:10,] #just the first 10



print(nrow(xsig)/nrow(x.pwu)) #percent significant pathways


# visualizing the most significant pathways
ksig=x.kf[x.kf$we.ep<0.05,]
ksig$path=rownames(ksig)
ksig2=arrange(ksig,desc(abs(effect)))
ksig3=ksig2[1:10,] #just the first 10



print(nrow(ksig)/nrow(x.kf)) #percent significant pathways

library(ggplot2)
p=ggplot(ksig3,aes(x=path,y=effect))+geom_bar(stat="identity")+coord_flip()+theme_classic()
p


###### Maaslin2 for picrust2 output


# data trasnformation to CRL
library(microbiome)
?microbiome::transform
kegg_relative<-microbiome::transform(kegg,"compositional")
otu_kegg_relative<-t(data.frame(otu_table(kegg_relative))) # relative otu table (******)

otu_kegg_clr<-t(data.frame(otu_table(kegg_clr))) # clr transformed otu table

pwus_relative<-microbiome::transform(pw_us,"compositional")
otu_pwus_relative<-t(data.frame(otu_table(pwus_relative))) # relative otu table (******)

pwus_clr=microbiome::transform(pw_us,"clr")
otu_pwus_clr<-t(data.frame(otu_table(pwus_clr))) # clr transformed otu table


metadata_0m<-data.frame(sample_data(pw_us))

write.csv(metadata_0m,"metadata_0m.csv")
write.csv(otu_pwus_relative,"otu_pwus_relative.csv")
write.csv(otu_kegg_relative,"otu_kegg_relative.csv")

getwd()
input_metadata_0m<-read.csv("/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/metadata_0m.csv",header=TRUE,sep=",",row.names = 1)


 ## pwus_relative
input_pwus_r <- data.frame(otu_pwus_relative)


library(Maaslin2)

# Rejection_1YA

mas_0m <- Maaslin2(
  input_pwus_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_pwus_r_R1Y",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = "Rejection_1YA",
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_pwus_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_pwus_r_R1Y_G_Agesex",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_pwus_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_pwus_r_R1Y_G_AgesexG_BMI",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_pwus_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_pwus_r_R1Y_G_AgesexG_BMI_DM",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_pwus_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_pwus_r_R1Y_G_AgesexG_BMI_DM_Abts",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_pwus_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_pwus_r_R1Y_G_AgesexG_BMI_DM_Abts_ISA",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before","ISA_before"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts","ISA_before,No ISA"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_pwus_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_pwus_r_R1Y_G_AgesexG_BMI_DM_Abts_ISA_desensitization",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before","ISA_before","desensitization"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts","ISA_before,No ISA","desensitization, No Desensitization"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

# transform = 'AST' (default=LOG)

mas_0m_1 <- Maaslin2(
  input_pwus_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_AST_pwus_r_R1Y_G_AgesexG_BMI_DM_Abts_ISA_desensitization",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  transform = 'AST',
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before","ISA_before","desensitization"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts","ISA_before,No ISA","desensitization, No Desensitization"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)


## pwus_clr (Result: No difference)
input_pwus_clr <- data.frame(otu_pwus_clr)

# Rejection_1YA

mas_0m <- Maaslin2(
  input_pwus_clr,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_pwus_clr_R1Y",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = "Rejection_1YA",
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_pwus_clr,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_pwus_clr_R1Y_G_Agesex",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_pwus_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_pwus_r_R1Y_G_AgesexG_BMI",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_pwus_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_pwus_r_R1Y_G_AgesexG_BMI_DM",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_pwus_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_pwus_r_R1Y_G_AgesexG_BMI_DM_Abts",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_pwus_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_pwus_r_R1Y_G_AgesexG_BMI_DM_Abts_ISA",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before","ISA_before"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts","ISA_before,No ISA"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_pwus_clr,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_pwus_clr_R1Y_G_AgesexG_BMI_DM_Abts_ISA_desensitization",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before","ISA_before","desensitization"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts","ISA_before,No ISA","desensitization, No Desensitization"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

## otu_kegg_relative
input_kegg_r <- data.frame(otu_kegg_relative)

# Rejection_1YA

mas_0m <- Maaslin2(
  input_kegg_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_kegg_r_R1Y",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = "Rejection_1YA",
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_kegg_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_kegg_r_R1Y_G_Agesex",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_kegg_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_kegg_r_R1Y_G_AgesexG_BMI",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_kegg_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_kegg_r_R1Y_G_AgesexG_BMI_DM",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_kegg_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_kegg_r_R1Y_G_AgesexG_BMI_DM_Abts",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_kegg_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_kegg_r_R1Y_G_AgesexG_BMI_DM_Abts_ISA",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before","ISA_before"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts","ISA_before,No ISA"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_kegg_r,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_kegg_r_R1Y_G_AgesexG_BMI_DM_Abts_ISA_desensitization",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before","ISA_before","desensitization"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts","ISA_before,No ISA","desensitization, No Desensitization"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)


# otu_kegg_clr (result: No difference)

input_kegg_clr <- data.frame(otu_kegg_clr)

mas_0m_1 <- Maaslin2(
  input_kegg_clr,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_kegg_clr_R1Y_G_AgesexG_BMI_DM_Abts",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_kegg_clr,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_kegg_clr_R1Y_G_AgesexG_BMI_DM_Abts_ISA",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before","ISA_before"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts","ISA_before,No ISA"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)

mas_0m_1 <- Maaslin2(
  input_kegg_clr,
  input_metadata_0m,
  output = "/Users/Chohyunjeong/Dropbox/KTPL_2021/Qiime_202108/Picrust2/Maaslin2_kegg_clr_R1Y_G_AgesexG_BMI_DM_Abts_desensitization",
  min_abundance = 0.001,
  min_prevalence = 0.1,
  fixed_effects = c("Rejection_1YA","G_Age","sex","G_BMI","DM","Abts_before","ISA_before","desensitization"),
  reference=c("Rejection_1YA,No Rejection", "G_Age,Twenty","sex,Male","G_BMI,Normal","DM,No DM","Abts_before, No Abts","ISA_before,No ISA","desensitization, No Desensitization"),
  correction = "BH",
  max_significance = 0.25,
  normalization = 'NONE',
  standardize = FALSE)



??masslin2

class(rich_T)
summary(rich_T)
class(mdata)
 
merge()

ggplot(data=metadata_0m)+geom_histogram(aes(x=dialysis_duration_M))

summary(metadata_0m)

library(dplyr)
metadata_0m=metadata_0m %>% dplyr::mutate(G_dialysis = ifelse(dialysis_duration_M>=60,"> 5yrs",ifelse(dialysis_duration_M >=3,"< 5yrs","Preemptive")))


## G_dialysis

 ## Change metadata ****
library(phyloseq)
tax_0=tax_table(ps_0m_rare)
otu_0=otu_table(ps_0m_rare)
sam_0=metadata_0m
tree_0=phy_tree(ps_0m_rare)
otutax_0 = phyloseq(otu_0, tax_0)
ps_0m_rare_new=merge_phyloseq(otutax_0,sample_data(sam_0),tree_0)
ps_0m_rare
ps_0m_rare_new



sample_data(ps_0m_rare_new)$G_dialysis <- factor(
  sample_data(ps_0m_rare_new)$G_dialysis, 
  levels = c("Preemptive", "< 5yrs", "> 5yrs"))

(p <- plot_richness(ps_0m_rare_new, x = "G_dialysis", color = "G_dialysis", measures = c("Observed", "Shannon")))
p+geom_boxplot()
p

kruskal.test(rich_0m$Observed~sample_data(ps_0m_rare_new)$G_dialysis)
pairwise.wilcox.test(rich_0m$Observed,sample_data(ps_0m_rare_new)$G_dialysis,p.adjust.method="fdr")

kruskal.test(rich_0m$Shannon~sample_data(ps_0m_rare_new)$G_dialysis)
pairwise.wilcox.test(rich_0m$Shannon,sample_data(ps_0m_rare_new)$G_dialysis,p.adjust.method="fdr")

sample_colors <- c("No rejection"="#2874C5", "Rejection"="#EABF00")
sample_types <- c("Preemptive", "< 5yrs", "> 5yrs")
sample_labels <- c("Preemptive"="Preemptive", "< 5yrs"="< 5yrs", "> 5yrs"="> 5yrs")

p<-plot_richness(ps_0m_rare_new,
                 x = "G_dialysis",   #compare diversity of datatype
                 measures = c("Observed", "Shannon")) +  #choose diversity measures
  geom_violin(aes(fill = "G_dialysis"), show.legend = FALSE)+ #make violin plot, set fill aes to Time
  geom_boxplot(width=0.1) + #add boxplot, set width
  theme_classic()+       #change theme to classic
  xlab(NULL)+    #no label on x-axis
  theme(axis.text.y.left = element_text(size = 10),  #adjust y-axis text
        axis.text.x = element_text(size = 15, hjust = 0.5, vjust = 1, angle = 0),           #adjust x-axis label position
        axis.title.y = element_text(size = 15))+  #adjust y-axis title
  theme(strip.text = element_text(face = "bold", size = 15))+ #adjust headings
  scale_fill_manual(values= sample_colors) +   #set fill colors
  scale_x_discrete(          #change x-axis labels
    breaks = sample_types, 
    labels = sample_labels)+                   
  ggtitle("Alpha Diversity") +        #add title
  theme(plot.title=element_text(size = 15, face = "bold", hjust = 0.5)) #change title size, face and position
p


## G_dialysis

sample_colors <- c("Preemptive"="coral","< 5yrs"="cadetblue1","> 5yrs"="darkviolet")

p<-plot_ordination(
  physeq = ps_0m_rare_new,           #phyloseq object
  ordination = pcoa_0m)+     #ordination
  geom_point(aes(fill = G_dialysis, shape = G_dialysis), size = 3) +  #sets fill color to Time
  scale_shape_manual(values = c(21, 22, 23))+
  scale_fill_manual(values = sample_colors) +
  theme_classic() +     #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 15),    #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command
p

# Calculate bray curtis distance matrix: permanova
set.seed(1)

# Calculate bray curtis distance matrix
ps_bray_0m<- phyloseq::distance(ps_0m_rare_new, method = "bray")

# Adonis test
library(vegan)
adonis(ps_bray_0m ~ sample_data(ps_0m_rare_new)$G_dialysis) # p 0.328

# Homogeneity of dispersion test
beta <- betadisper(ps_bray_0m, sample_data(ps_0m_rare_new)$G_dialysis)
permutest(beta) # p 0.048


# crosstable
library(moonBook)

library(dplyr)
metadata_0m<-rename(metadata_0m, Rejection=Rejection_1YA,Borderline=Rejection_1YA_B, Al_B=Rejection_1YA_all.B )
metadata_0m<-rename(metadata_0m, Borderline=Rejection_1YA_B, Al_B=Rejection_1YA_all.B )

table_R<-mytable(Rejection~., data=metadata_0m) 
str(metadata_0m)

mycsv(table_R,file="table_R.csv")

table1<-mytable(Rejection_1YA~., data=input_metadata_0m)
mycsv(table1,file="table1.csv")

table3<-mytable(Rejection_1YA_all.B~., data=input_metadata_0m)
mycsv(table3,file="table3.csv")
