# load libraries
library("microseq")
library(readxl)
library(dplyr)
library(stringr)
library(dada2)
library("Biostrings")
library(phyloseq)
library(tidyr)
library(EcolUtils)
library(vegan)

# -------------
# 1) 18S data
# -------------

## Data import
# Fasta:
seqs<-readDNAStringSet("~/Desktop/oneblue/DATA/METAB/oneblue_18S.fasta")

# Raw file from dada2 by Ramiro:
ab<-as.data.frame(readRDS("~/Desktop/oneblue/DATA/METAB/oneblue_18S_seqtab_final.rds"))
colnames(ab)<-names(seqs)
rownames(ab)<-gsub("-18S_S1_L001_trimmed","",rownames(ab))

## Taxonomy file for Eukaryotes
euk_tax<-as.data.frame(readRDS("Desktop/oneblue/DATA/METAB/oneblue_18S_taxo_MB80.RData"));colnames(euk_tax)<-c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")

# OTU names
# last taxonomic rank available
euk_tax$Rank<-apply(euk_tax, 1, function(x) {last_col <- tail(names(x)[!is.na(x)], 1)
if (length(last_col) == 0) return(NA)
return(last_col)
})
euk_tax$Taxa <- apply(euk_tax[,1:9], 1, function(x) {
  non_na_vals <- x[!is.na(x)]
  if (length(non_na_vals) == 0) return(NA)
  return(tail(non_na_vals, 1))
})
euk_tax$asv<-names(seqs)
rownames(euk_tax)<-euk_tax$asv

## Phylogeny
euk_tree<-read_tree("~/Desktop/oneblue/DATA/METAB/oneblue_18S.tree")

## Samples info
info <- read.csv("~/Desktop/oneblue/DATA/METADATA/metadata_oneblue.csv", sep = ",", header = TRUE)
info$Date<-janitor::excel_numeric_to_date(info$Date)
info_18S<-info[info$MB_18S_ID %in% rownames(ab),-c(1,2,9)]
rownames(info_18S)<-info_18S$MB_18S_ID

# Final Phyloseq
ps18S<-phyloseq(otu_table(ab, taxa_are_rows = FALSE),
                sample_data(info_18S),
                phy_tree(euk_tree),
                tax_table(as.matrix(euk_tax)),
                refseq(readDNAStringSet("Desktop/oneblue/DATA/METAB/oneblue_18S.fasta"))
)

## Clean Phyloseq based on the control
ab <- as.data.frame(otu_table(ps18S))
asv.ctrl<-colnames(ab)[which(ab["BPCR5",]>0)] # which ASVs are present in the PCR control?
euk_tax[asv.ctrl,] # Their taxonomy
# beeeehehehe
# => there's just one, we remove it in the following cleaning steps based on taxonomy (Metazoa)

## Clean Phyloseq based on taxonomy
# Only Eukaryotes
table(euk_tax$Domain)
ps18S.c <- subset_taxa(ps18S, Domain %in% c("Eukaryota"))

# remove nuclear DNA
table(euk_tax$Supergroup)
ps18S.c <- subset_taxa(ps18S.c, !Supergroup %in% c("Cryptista:nucl", "TSAR:nucl"))

# remove Metazoan? (multicellular animals) remove fungi?
table(euk_tax$Subdivision)
euk_tax[euk_tax$Subdivision %in% "Metazoa",]
# There's some Homo sapiens in here, have a look maybe this is your DNA!
# There's also marine specific signal with Sponges, Zooplankton, worms and annelides, Coral, Molluscs (we usually remove them)
ps18S.c <- subset_taxa(ps18S.c, !Subdivision %in% c("Metazoa"))
# Removing Fungi is up to you I like to keep them
#ps18S.c <- subset_taxa(ps18S.c, !Subdivision %in% c("Fungi")) 

# remove multicellular plants and macro-algae?
table(euk_tax$Division)
ps18S.c <- subset_taxa(ps18S.c, !Division %in% c("Streptophyta")) # terrestrial plants
euk_tax[euk_tax$Division %in% "Rhodophyta",] # Rhodophyta may contain macro red-algae but Oltmannsiellopsis_viridis is a colonial protist
euk_tax[euk_tax$Class %in% "Ulvophyceae",]   # Ulvophyceae may contain macro green-algae but Rhodella and Porphyridiales are unicellular
# "Phaeophyceae" are brown-algae absent

## Remove the PCR controls
ps18S.c <- prune_samples(!sample_names(ps18S.c) %in% c("BPCR5"), ps18S.c)
ps18S.c <- prune_taxa(taxa_sums(ps18S.c) > 0, ps18S.c)

## Rarefying: normalize the total number of reads per samples
# convert the table in the required format 
sp18S<-as.data.frame(otu_table(ps18S.c))
# order samples by total number of reads
rowSums(sp18S)[order(rowSums(sp18S))]
# => the number of reads per sample is wuite homogeneous (47,000 to 103,000), 
# if it was up to me I wouldn't rarefy at all (not to lose rare ASVs)
sp18S_rare<-rrarefy.perm(sp18S, round.out = TRUE, sample = 47822, n = 99)
rowSums(sp18S_rare)

## Final Cleaned phyloseq object
rare_otu_phy <- otu_table(sp18S_rare, taxa_are_rows = FALSE)
ps18S.rarefied <- ps18S.c
otu_table(ps18S.rarefied) <- rare_otu_phy
ps18S.rarefied <- prune_taxa(taxa_sums(ps18S.rarefied) > 0, ps18S.rarefied)
sample_sums(ps18S.rarefied)
sample_data(ps18S.rarefied)
saveRDS(ps18S.rarefied, "~/Desktop/oneblue/DATA/METAB/phyloseq_oneblue_18S.rds")

# -------------
# 2) 16S data
# -------------

## Data import
# Fasta:
seqs<-readDNAStringSet("~/Desktop/oneblue/DATA/METAB/oneblue_16S.fasta")

# Raw file from dada2 by Ramiro:
ab<-as.data.frame(readRDS("~/Desktop/oneblue/DATA/METAB/oneblue_16S_seqtab_final.rds"))
colnames(ab)<-names(seqs)
rownames(ab)<-gsub("-16S_S1_L001_trimmed","",rownames(ab))

## Taxonomy file for Prokaryotes
prok_tax<-as.data.frame(readRDS("~/Desktop/oneblue/DATA/METAB/oneblue_16S_taxo_MB80.RData"))
prok_tax$Species<-ifelse(is.na(prok_tax$Species) == FALSE, paste(prok_tax$Genus,prok_tax$Species, sep = " "), NA)

# OTU names
# last taxonomic rank available
prok_tax$Rank<-apply(prok_tax, 1, function(x) {last_col <- tail(names(x)[!is.na(x)], 1)
if (length(last_col) == 0) return(NA)
return(last_col)
})
prok_tax$Taxa <- apply(prok_tax[,1:7], 1, function(x) {
  non_na_vals <- x[!is.na(x)]
  if (length(non_na_vals) == 0) return(NA)
  return(tail(non_na_vals, 1))
})
prok_tax$asv<-names(seqs)
rownames(prok_tax)<-prok_tax$asv

## Phylogeny
prok_tree<-read_tree("~/Desktop/oneblue/DATA/METAB/oneblue_16S.tree")

## Samples info
info <- read.csv("~/Desktop/oneblue/DATA/METADATA/metadata_oneblue.csv", sep = ",", header = TRUE)
info$Date<-janitor::excel_numeric_to_date(info$Date)
info_16S<-info[info$MB_16S_ID %in% rownames(ab),-c(2,3,9)]
rownames(info_16S)<-info_16S$MB_16S_ID

# Final Phyloseq
ps16S<-phyloseq(otu_table(ab, taxa_are_rows = FALSE),
                sample_data(info_16S),
                phy_tree(prok_tree),
                tax_table(as.matrix(prok_tax)),
                refseq(readDNAStringSet("~/Desktop/oneblue/DATA/METAB/oneblue_16S.fasta"))
)

## Clean Phyloseq based on the control
ab <- as.data.frame(otu_table(ps16S))
bpcr5<-colnames(ab)[which(ab[rownames(ab) %in% c("BPCR5"),]>0)] # which ASVs are present in the PCR controls?
bpcr6<-colnames(ab)[which(ab[rownames(ab) %in% c("BPCR6"),]>0)] # which ASVs are present in the PCR controls?
blank<-colnames(ab)[which(ab[rownames(ab) %in% c("BLANK"),]>0)] # which ASVs are present in the PCR controls?
asv.ctrl<-unique(c(bpcr5, bpcr6, blank))
prok_tax[asv.ctrl,] # Their taxonomy
# This time we have more than one asv in the controls, and also very abundant asvs

# let's create a file to explore the abundance of ASVs present in controls
ab.16S<-as.data.frame(otu_table(ps16S))
ab.16S[,asv.ctrl]
ctrl<-c("BPCR5", "BPCR6", "BLANK")
# mean abudances in the controls vs the environment
ctrl.ab<-data.frame(ctrl = apply(ab.16S[ctrl,asv.ctrl], 2, mean),
           env = apply(ab.16S[!rownames(ab.16S) %in% ctrl, asv.ctrl], 2, mean))
ctrl.ab$ratio<-ctrl.ab$ctrl/ctrl.ab$env
ctrl.ab<-merge(ctrl.ab, prok_tax, by = "row.names")

# These are taxa known to come from the human skin microbiome
# We could decide to remove them if they are not super abundant 
# but they could also be an environmental signal and not lab-contamination
skin.genus<-c("Enhydrobacter", "Cutibacterium", "Staphylococcus", "Pseudomonas", "Acinetobacter", "Corynebacterium")
skin.fam<-c("Neisseriaceae")
skin.asv<-prok_tax[prok_tax$Genus %in% c(skin.genus,skin.fam),"asv"]
skin.ab<-data.frame(ctrl = apply(ab.16S[ctrl,skin.asv], 2, max),
                    env = apply(ab.16S[!rownames(ab.16S) %in% ctrl, skin.asv], 2, max))

# Okay let's now filter thoss asv more abundant in the controls than in the env. samples
ps16S.c <- subset_taxa(ps16S, !asv %in% ctrl.ab[ctrl.ab$ratio > 0.5,"asv"] )

## Clean Phyloseq based on taxonomy
# Only Prokaryotes
table(prok_tax$Kingdom)
ps16S.c <- subset_taxa(ps16S.c, !Kingdom %in% c("Eukaryota"))

# remove organelles DNA
table(prok_tax$Order)
ps16S.c <- subset_taxa(ps16S.c, !Order %in% c("Chloroplast"))
table(prok_tax$Family)
ps16S.c <- subset_taxa(ps16S.c, !Family %in% c("Mitochondria"))

## Remove the PCR controls
ps16S.c <- prune_samples(!sample_names(ps16S.c) %in% c("BPCR5", "BPCR6", "BLANK"), ps16S.c)
ps16S.c <- prune_taxa(taxa_sums(ps16S.c) > 0, ps16S.c)

## Rarefying: normalize the total number of reads per samples
# convert the table in the required format 
sp16S<-as.data.frame(otu_table(ps16S.c))
# order samples by total number of reads
rowSums(sp16S)[order(rowSums(sp16S))]
# => the number of reads per sample is really heterogeneous (12,000 to 161,000), 
# if it was up to me I wouldn't rarefy at all (not to lose rare ASVs)
sp16S_rare<-rrarefy.perm(sp16S, round.out = TRUE, sample = 12318, n = 99)
rowSums(sp16S_rare)

## Final Cleaned phyloseq object
rare_otu_phy <- otu_table(sp16S_rare, taxa_are_rows = FALSE)
ps16S.rarefied <- ps16S.c
otu_table(ps16S.rarefied) <- rare_otu_phy
ps16S.rarefied <- prune_taxa(taxa_sums(ps16S.rarefied) > 0, ps16S.rarefied)
saveRDS(ps16S.rarefied, "~/Desktop/oneblue/DATA/METAB/phyloseq_oneblue_16S.rds")

