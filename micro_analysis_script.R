install.packages("tidyMicro")
library(tidyMicro)
library(magrittr)

atlas_otu <- read.csv("otu_family_input.csv", header=T, check.names=F)

df_atlas_otu <- as.data.frame(atlas_otu)

pheno <- read.table("clinical_data.txt", header = T)

clinical.pheno <- na.omit(pheno)

micro.set <- tidy_micro(otu_tabs = df_atlas_otu, clinical = clinical.pheno, tab_names = "Family")

taxanomy_file <- read.csv("taxanomy_table.csv",header=T)

unique_sub <- unique(micro.set$subject)

# number of subjects in micro.set 
length(unique_sub)
cat("The number of subjects in the micro.set is:", length(unique_sub))

# number of otu assayed in microbiome analysis
otu <- length(taxanomy_file$otu_id)
cat("The number of total assayed otu in the microbiome analysis is:", otu)

# filter out male samples
micro.male <- micro.set%<>%
  filter(sex == "male")

# number of male subjects
male_uni <- unique(micro.male$subject)
length(male_uni)
cat("Number of male subjects in micro.set:", length(male_uni))

# remove and reload micro.set
remove(micro.set)
micro.set <- tidy_micro(otu_tabs = df_atlas_otu, clinical = clinical.pheno, tab_names = "Family")

# filter out female samples
micro.female <- micro.set%<>%
  filter(sex == "female")

# number of female subjects
female_uni <- unique(micro.female$subject)
length(female_uni)
cat("Number of female subjects in micro.set:", length(female_uni))

# remove and reload micro.set
remove(micro.set)
micro.set <- tidy_micro(otu_tabs = df_atlas_otu, clinical = clinical.pheno, tab_names = "Family")

#filter out scandinavian subjects
micro.scan <- micro.set%<>%
  filter(nationality == "Scandinavia")

# number of scandinavian subjects
uni_scan <- unique(micro.scan$subject)
length(uni_scan)
cat("The number of scandinavian subjects in this study are:", length(uni_scan))

# remove and reload micro.set
remove(micro.set)
micro.set <- tidy_micro(otu_tabs = df_atlas_otu, clinical = clinical.pheno, tab_names = "Family")

#PCA plot
pdf("A5_PCA_BMI_Plot_CM.pdf")
micro.set%>% micro_pca(table = "Family", grp_var = bmi_group, legend_title = "BMI groups")
dev.off()

#Relative Abundance 

pdf("A5_Relative_Abundance_Bar_Plot_CM.pdf")
ra_bars(
  micro.set,
  table = "Family",
  bmi_group,
  top_taxa = 5,
  ylab = "% RA", xlab = "BMI Groups", main = "Stacked Bar Charts of Taxa by BMI Group") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

cat("The 5 taxa that have the highest relative abundance based on BMI group are:\n 1.Actinobacteria\n 2.Bacteroidetes\n 3.Clostridium cluster IV\n 4.Clostridium cluster XIVa\n 5.Uncultured Clostridiales\n")


# Box Plot
pdf("A5_Taxa_BoxPlot_BMI_CM.pdf")
taxa <- "Bacteroidetes"

taxa_boxplot(
  micro.set,
  taxa = taxa,
  bmi_group,
  xlab = "BMI Groups",
  ylab = "Relative Abundance",
  main = "Box Plot of Relative Abundance of Taxa Based on BMI"
)

dev.off()

cat("The morbidly obese BMI group demonstrates high abundance of the Bacteroidetes.")

# Univariate Comparisons

micro.filt <- micro.set %>%
  otu_filter(prev_cutoff = 1, ra_cutoff = 0.6)

nb_ord <- micro.filt%>%
  nb_mods(table = 'Family', bmi_group)

pdf("A5_Rocky_Mountain_Plot_CM.pdf")
nb_ord%>%
  micro_rocky_mtn(bmi_group, xlab = "Taxa", main = "Rocky Mountain Plot", subtitle = "Direction of bar indicates direction of relationship")
dev.off()

cat("Bacteroidetes and Proteobacteria are positively associated in the severe obese BMI group.")

# Two group comparison
lst_t <- list()
lst_w <- list()
candidates <- c("Actinobacteria", "Bacteroidetes","Fusobacteria")

# t-test
for (microbes in candidates) {
  micro.candidate <- micro.filt[micro.filt$Taxa == microbes, ]
  micro.under <- micro.candidate[micro.candidate$bmi_group == "underweight", ]
  micro.sev <- micro.candidate[micro.candidate$bmi_group == "severeobese", ]
  
  lst_t[[microbes]] <- t.test(micro.under$ra, micro.sev$ra, var.equal=T)
}

# Wilcox non-parametric test
for (microbes in candidates) {
  micro.candidate <- micro.filt[micro.filt$Taxa == microbes, ]
  micro.under <- micro.candidate[micro.candidate$bmi_group == "underweight", ]
  micro.sev <- micro.candidate[micro.candidate$bmi_group == "severeobese", ]
  lst_w[[microbes]] <- wilcox.test(micro.under$ra, micro.sev$ra)
  
}

cat("Bacteroidetes showed a significant difference between underweight and severely obese BMI groups with a P value of: 0.0178")