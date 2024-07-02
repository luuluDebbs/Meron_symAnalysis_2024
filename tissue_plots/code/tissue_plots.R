library(qiime2R)
setwd("~/OneDrive/Desktop/minimalR/tissue_plots")

SVs <- read_qza("table.qza")
head(SVs$type)
SVs$uuid
SVs$contents

metadata <- read_q2metadata("sample-metadata.txt")

taxonomy <- read_qza("blast-taxonomy-results.qza")

physeq <- qza_to_phyloseq(
  features = "table.qza",
  tree = "rooted_tree.qza",
  "blast-taxonomy-results.qza",
  metadata = "sample-metadata.txt")

library(tidyverse)

shannon <- read_qza("shannon_vector.qza")

shannon <- shannon$data %>% rownames_to_column("SampleID")

install.packages("ggplot")

ggplots::venn(list(metadata=metadata$SampleID, shannon=shannon$SampleID))

metadata <- 
  metadata %>% 
  left_join(shannon)
head(metadata)

target_genus <- c("Acropora", "Alveopora", "Galaxea", "Leptoseris", "Montipora",
                  "Pocillopora")
metadata %>% 
  group_by(Genus) %>% 
  count()


metadata <- metadata %>%
  mutate(alpha_genus= Genus == target_genus)

metadata1 <- droplevels(metadata[metadata$Genus == "Acropora"|metadata$Genus == "Alveopora"|
                                   metadata$Genus == "Galaxea"|metadata$Genus == "Leptoseris"|
                                   metadata$Genus == "Montipora"|metadata$Genus == "Pocillopora",])

summary(metadata1)

Meronthemealpha <-  theme(panel.border = element_rect(colour = "black", fill = NA),
                     panel.background = element_blank(),
                     panel.grid = element_blank(),
                     axis.text  = element_text(size = 14), 
                     axis.text.x = element_text(angle = 0, face = "italic"),
                     axis.title = element_text(size = 20),
                     strip.background = element_blank(),
                     strip.text = element_text(size = 10),
                     axis.line.x = element_line(size = 0.5),
                     axis.line.y = element_line(size = 0.5),
                     legend.position = "right",
                     legend.text =element_text(size =10))


metadata1 %>% 
  filter(!is.na(shannon_entropy), Type != "Planula") %>%
  ggplot(aes(Genus, shannon_entropy)) +
  geom_boxplot(position = "dodge")+
  geom_jitter(aes(fill = Depth), shape=21, width=0.2, height = 0, size = 2.5)+
  xlab("Coral Genera")+
  ylab("Shannon Diversity")+
  Meronthemealpha+
  scale_fill_gradient(low = "white", high = "darkblue")

#Beta diverisy PCoA

metadata<- read_q2metadata("sample-metadata.txt")
uwunifrac<-read_qza("unweighted_unifrac_pcoa_results.qza")
shannon<- read_qza("shannon_vector.qza")$data %>% 
  rownames_to_column("SampleID")
unifrac_modify<- uwunifrac$data$Vectors %>% 
  select(SampleID, PC1, PC2) %>% 
  left_join(metadata) %>% 
  left_join(shannon) 


unifrac_modify %>% 
  ggplot(aes(PC1, PC2, color= `Type`)) +
  geom_point(alpha=0.5, size=5) +
  stat_ellipse() +
  xlab(paste("PC1: ", round(100*uwunifrac$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*uwunifrac$data$ProportionExplained[2]), "%")) +
  theme_bw() +
  scale_shape_manual(values = c(17,16), name = "Region" ) +
  scale_color_discrete(name="Type")

unifrac_modify %>% 
  ggplot(aes(PC1, PC2, color= `Type`,shape= 'Region')) +
  geom_point(aes(shape= Type),alpha=0.5, size=5) +
  stat_ellipse() +
  xlab(paste("PC1: ", round(100*uwunifrac$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*uwunifrac$data$ProportionExplained[2]), "%")) +
  theme_bw() +
  #scale_size_continuous(name = "Shannon Diversity") +
  scale_shape_discrete(name ="Region", palette() ) +
  scale_color_discrete(name="Type")

unifrac_modify %>% 
  ggplot(aes(x=PC1, y=PC2, color= Depth))+
  geom_point(alpha=0.5, size = 5) +
  xlab(paste("PC1: ", round(100*uwunifrac$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*uwunifrac$data$ProportionExplained[2]), "%"))+
  #stat_ellipse()+
  #scale_shape_discrete(name ="Region")+
  #scale_size_continuous(name = "Shannon Diversity")+
  scale_color_continuous(name = "Depth in meters")+
  Merontheme

#phyloseq
library(phyloseq)

plot_bar(physeq, fill = "Genus") 

#build tree

plot_tree(physeq, color = "environment", label.tips = "taxa_names", ladderize = "left",
          plot.margin = 0.3)

plot_heatmap(physeq, taxa.label="Genus")


#alpha diversity

richness<- plot_richness(physeq,x="Type", color = "Region", measures = c("Chao1", "Shannon"))

richness <- richness + geom_point(size=4, alpha=0.7)

richness

ggsave("richenss.svg,")

plot_richness(physeq,x="sample_type", color = "depth", measures = c("Chao1", "Shannon"))

plot_richness(physeq,x="environment", color = "sample_type",measures = c("Chao1", "Shannon"))

plot_richness(physeq,x="depth", measures = c("Chao1", "Shannon"))


library(microbiome)
library(phyloseq)
library(knitr)
library(tidyverse)
#Get tibbles

otu.tib <- otu_tibble(physeq, column.id = "FeatureID")
tax.tib <- tax_tibble(physeq, column.id = "FeatureID")
sample_tib <- sample_tibble(physeq, column.id = "SampleID")

dat <- otu.tib %>% 
  pivot_longer(-FeatureID, names_to = "SampleID", values_to = "count")
dat

dat <- dat %>% 
  left_join(tax.tib, by="FeatureID")
dat

dat <- dat %>% 
  left_join(sample_tib, by = "SampleID")
dat

dat

class(dat)
write_csv(dat, file = "tissue_metadat.csv")
write_csv(tax.tib, file = "taxa_metadat.csv")

dat %>% 
  ggplot(aes(SampleID, count)) +
  facet_grid(~Type, scales = "free_x", space = "free_x")+
  geom_bar(aes (fill= Genus.x), stat = "identity", position = "fill", width = 0.9)

Merontheme <-  theme(panel.border = element_rect(colour = "black", fill = NA),
                     panel.background = element_blank(),
                     panel.grid = element_blank(),
                     axis.text  = element_text(size = 10), 
                     axis.text.x = element_text(angle = 90, face = "italic"),
                     axis.text.x.top = element_text(size=10),
                     axis.title = element_text(size = 20),
                     strip.background = element_blank(),
                     strip.text = element_text(size = 10),
                     axis.line.x = element_line(size = 0.5),
                     axis.line.y = element_line(size = 0.5),
                     legend.position = "right",
                     legend.text =element_text(size =8, face = "italic"))


dat <- dat %>% dat <- dat %>% NULL
  mutate(depth2= case_when(Depth <= 3 ~ "Shallow (1-3m)",
                           Depth >= 4 & Depth <= 10 ~ "Deeper Shallow (4-10m)",
                           Depth >= 11 & Depth <= 29 ~ "11-29m",
                           Depth >= 30 & Depth <= 50 ~ "Upper Mesophotic (30-50m)",
                           Depth >= 51 & Depth <= 70 ~ "51-70m",
                           Depth >70 ~ "Lower Mesophotic(>70)"))

dat$depth2 <- as_factor(dat$depth1)

levels(dat$depth1)

string <- ("C1-AB778611.1")
dat <- dat %>% 
  mutate(species_name= (str_replace_all(pattern = '-.*', replacement = "", Species.x)))

dat <- dat %>% 
  mutate(species_name2= (sub("\\..*", "", species_name)))

dat <- dat %>% 
  mutate(species_name3= (gsub("(\\d)[a-zA-Z]", "\\1", species_name2)))


dat$depth2 <- factor(dat$depth2, levels = c("Shallow (1-3m)", "Deeper Shallow (4-10m)",
                                            "11-29m", "Upper Mesophotic (30-50m)",
                                            "51-70m", "Lower Mesophotic(>70)"))


#merging rare taxa
  
dat %>% 
  filter(Type != "Planula") %>% 
  ggplot(aes(Genus.y, count)) +
  facet_grid(~depth2 , scales = "free_x", space = "free_x")+
  geom_bar(aes (fill= Species_name3), stat = "identity", position = "fill", width = 0.9)+
  scale_fill_discrete(name = NULL)+
  scale_fill_manual(values= c("#771155","#AA4488", "#774411", "#114477", "#777711", "#DDDD77",
                              "#4477AA", "#77AADD", "#117777","#44AAAA", "#77CCCC", "#117744",
                              "#44AA77", "#88CCAA","#CC99BB", "#AA7744", "#DDAA77", "#771122",
                              "#AAAA44", "#AA4455","#DD7788", "#00008F", "#0000FF", "#0070FF"), name = "ITS2 Types" )+
  scale_y_continuous(name = "Relative abundance",
                     labels = scales::percent) +
  labs(x = NULL)+
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 12))+
  Merontheme  
  
  
  
shannon <- read_qza("shannon_vector.qza")

shannon <- shannon$data %>% rownames_to_column("SampleID")

dat <- dat %>% 
  left_join(shannon, by = "SampleID")

dat %>%
  ggplot(aes(x = Genus.y, y = shannon_entropy, fill = Genus.y)) +
  geom_boxplot()+
  #geom_point()+
  labs(title = "Shannon Diversity Plot for Coral Tissue Samples", y = "Shannon Diversity", x = "Genera")+
  theme(axis.text.x = element_text(angle = 90)) +
  theme_classic()


uwunifrac<-read_qza("unweighted_unifrac_pcoa_results.qza")

unifrac_modify<- uwunifrac$data$Vectors %>% 
  select(SampleID, PC1, PC2)
dat <- dat %>% 
  left_join(unifrac_modify, by = "SampleID")

dat %>%
  filter(Type!= "Planula") %>% 
  ggplot(aes(PC1, PC2, color= `Genus.y`)) +
  geom_point(alpha=0.5, size=2) +
  #stat_ellipse() +
  xlab(paste("PC1: ", round(100*uwunifrac$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*uwunifrac$data$ProportionExplained[2]), "%")) +
  theme_bw() +
  #scale_size_continuous(name = "Shannon Diversity") +
  scale_shape_discrete(name ="Region", palette() ) +
  scale_color_discrete(name="Type")



library(vegan)

dat %>% group_by(SampleID) %>% 
  mutate(total = sum(count)) %>% 
  filter(total > 230) %>% 
  group_by(FeatureID) %>% 
  mutate(total = sum(count)) %>% 
  filter(total != 0)

