library("phyloseq")
library("ggplot2")
setwd("~/Google Drive/Hofmockel-agg-qiime-analysis/max-length/")
abundance_data <- read.delim(sep='\t', file="./summary2-orfs.txt",header=TRUE, strip.white=TRUE, row.names=1)
head(abundance_data)
summary(abundance_data)

#getting the core out of the raw data
min_count = 1
in_all_samples <- subset(abundance_data, abundance_data$CC_MI_H30 >=min_count & abundance_data$CC_MI_H39 >=min_count & abundance_data$PF_LM_H08 >=min_count &
                           abundance_data$PF_LM_H14 >=min_count & abundance_data$PF_LM_H16 >=min_count &
                           abundance_data$PF_LM_H03 >=min_count & abundance_data$PF_MI_H01 >=min_count &
                           abundance_data$PF_MI_H06 >=min_count & abundance_data$PF_MI_H12 >=min_count &
                           abundance_data$PF_MI_H13 >=min_count & abundance_data$PF_MM_H17 >=min_count &
                           abundance_data$PF_MM_H19 >=min_count & abundance_data$PF_MM_H20 >=min_count &
                           abundance_data$PF_SM_H02 >=min_count & abundance_data$PF_SM_H10 >=min_count &
                           abundance_data$PF_SM_H11 >=min_count & abundance_data$PF_WS_H04 >=min_count &
                           abundance_data$PF_WS_H07 >=min_count & abundance_data$PF_WS_H09 >=min_count &
                           abundance_data$PF_WS_H15 >=min_count & abundance_data$UP_MI_H41 >=min_count &
                           abundance_data$UP_MI_H43 >=min_count & abundance_data$UP_MI_H57 >=min_count &
                           abundance_data$UP_MI_H59 >=min_count)

dim(in_all_samples)
head(in_all_samples)
in_all_samples <- in_all_samples[-c(2,3,14,18)]
head(in_all_samples)

metadata = read.delim(file="./hof_metadata.txt", row.names = 1, sep="\t", header=TRUE)
ann_data <- read.delim(sep='\t', file="./observation_map.onecazy_filtered_cleaned2.txt",header=TRUE, strip.white=TRUE, row.names=1)

write.table(x=in_all_samples, file="./temp_to_write.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
abundance_data_norm <- read.delim(sep='\t', file="./summary_core_norm.txt",header=TRUE, strip.white=TRUE, row.names=1)


abundance_data_norm_matrix <- as.matrix(abundance_data_norm)
abundance <- otu_table(abundance_data_norm_matrix, taxa_are_rows=TRUE)
ann_data_matrix <- as.matrix(ann_data)
annotation <- tax_table(ann_data_matrix)
metadata <- sample_data(metadata)
all_agg <- phyloseq(metadata, annotation, abundance)
all_agg

#Looking at all core functions

#All soils
plot_bar(all_agg, "sample_name", fill="Cazy_fam2", facet_grid=~Cazy_fam2)

all_agg_just_pf = subset_samples(all_agg, crop_system == "prairieFert")
plot_bar(all_agg_just_pf, "sample_name", fill="Cazy_fam2", facet_grid=~Cazy_fam2)
tophits <- names(sort(taxa_sums(all_agg_just_pf),TRUE)[1:20])
all_agg_just_pf_top50 <- prune_species(tophits, all_agg_just_pf) 
plot_bar(all_agg_just_pf_top50, "agg_frac", fill="taxonomy", facet_grid=~Cazy_fam2)


all_agg_just_pf_just_micro = subset_samples(all_agg_just_pf, agg_frac == "micro")
plot_bar(all_agg_just_pf_just_micro, "sample_name", fill="Cazy_fam2", facet_grid=~Cazy_fam2)
plot_bar(all_agg_just_pf_just_micro, "Cazy_fam2", fill="Cazy_fam2")

tophits <- names(sort(taxa_sums(all_agg_just_pf_just_micro),TRUE)[1:20])
all_agg_just_pf_just_micro_top50 <- prune_species(tophits, all_agg_just_pf_just_micro)                 
plot_bar(all_agg_just_pf_just_micro_top50, "sample_name", fill="taxonomy", facet_grid=~Cazy_fam2)

all_agg_just_pf_just_micro = subset_samples(all_agg_just_pf, agg_frac == "micro")
plot_bar(all_agg_just_pf_just_micro, "sample_name", fill="Cazy_fam2", facet_grid=~Cazy_fam2)
plot_bar(all_agg_just_pf_just_micro, "sample_name", fill="taxonomy", facet_grid=~Cazy_fam2)

all_agg_just_pf_just_LM = subset_samples(all_agg_just_pf, agg_frac == "LM")
plot_bar(all_agg_just_pf_just_LM, "sample_name", fill="Cazy_fam2", facet_grid=~Cazy_fam2)
tophits <- names(sort(taxa_sums(all_agg_just_pf_just_LM),TRUE)[1:20])
all_agg_just_pf_just_LM_top50 <- prune_species(tophits, all_agg_just_pf_just_LM)                 
plot_bar(all_agg_just_pf_just_LM_top50, "sample_name", fill="taxonomy", facet_grid=~Cazy_fam2)

all_agg_just_pf_just_LM = subset_samples(all_agg_just_pf, agg_frac == "LM")
plot_bar(all_agg_just_pf_just_LM, "sample_name", fill="Cazy_fam2", facet_grid=~Cazy_fam2)
plot_bar(all_agg_just_pf_just_LM, "sample_name", fill="taxonomy", facet_grid=~Cazy_fam2)

all_agg_just_pf_just_LM = subset_samples(all_agg_just_pf, agg_frac == "LM")
plot_bar(all_agg_just_pf_just_LM, "sample_name", fill="Cazy_fam2", facet_grid=~Cazy_fam2)

tophits <- names(sort(taxa_sums(all_agg_just_pf_just_WS),TRUE)[1:20])
all_agg_just_pf_just_WS_top50 <- prune_species(tophits, all_agg_just_pf_just_WS)                 
plot_bar(all_agg_just_pf_just_WS_top50, "sample_name", fill="taxonomy", facet_grid=~Cazy_fam2)

all_agg_just_pf_just_WS = subset_samples(all_agg_just_pf, agg_frac == "WS")
plot_bar(all_agg_just_pf_just_WS, "sample_name", fill="Cazy_fam2", facet_grid=~Cazy_fam2)
plot_bar(all_agg_just_pf_just_WS, "sample_name", fill="taxonomy", facet_grid=~Cazy_fam2)




plot_bar(all_agg, "crop_system", fill="Cazy_fam2", facet_grid=~Cazy_fam2)
plot_bar(all_agg, "agg_frac", fill="Cazy_fam2", facet_grid=~Cazy_fam2)

tophits <- names(sort(taxa_sums(all_agg),TRUE)[1:20])
all_agg_top20 <- prune_species(tophits, all_agg)
all_agg_just_pf = subset_samples(all_agg, crop_system == "prairieFert")
tophits_just_pf <- names(sort(taxa_sums(all_agg_just_pf),TRUE)[1:20])
all_agg_just_pf_top20 <- prune_species(tophits_just_pf, all_agg_just_pf)
plot_bar(all_agg_top20, "crop_system", fill="Cazy_fam2", facet_grid=~Cazy_fam2)
plot_bar(all_agg_just_pf_top20, "agg_frac", fill="Cazy_fam2", facet_grid=~Cazy_fam2)
plot_bar(all_agg_top20, "crop_system", fill="taxonomy", facet_grid=~Cazy_fam2)
plot_bar(all_agg_just_pf_top20, "agg_frac", fill="Cazy_fam2")
plot_bar(all_agg_just_pf_top20, "agg_frac", fill = "Cazy_fam2", facet_grid = ~sample_block) +  theme(strip.text.x = element_text(size = 8, colour = "red", angle = 90))
plot_bar(all_agg_just_pf_top20, "agg_frac", fill="Cazy_fam2", facet_grid=~taxonomy)
plot_bar(all_agg_just_pf_top20, "agg_frac", fill="taxonomy", facet_grid=~Cazy_fam2)
plot_bar(all_agg_just_pf_top20, "agg_frac", fill="Cazy_fam2", facet_grid=~taxonomy)
plot_bar(all_agg_just_pf_top20, "agg_frac", fill="Cazy_fam2", facet_grid=agg_frac~taxonomy)
plot_bar(all_agg_just_pf_top20, "agg_frac", fill="Cazy_fam2", facet_grid=~taxonomy)