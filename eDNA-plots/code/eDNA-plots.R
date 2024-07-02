if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
setwd("~/OneDrive/Desktop/minimalR/eDNA-plots")

SVs <- read_qza("table.qza")
head(SVs$type)
SVs$uuid
SVs$contents

metadata <- read_q2metadata("sample-metadata.txt")

taxonomy <- read_qza("blast-taxonomy-results-97.qza")

physeq <- qza_to_phyloseq(
  features = "table.qza",
  tree = "rooted_tree.qza",
  "blast-taxonomy-results-97.qza",
  metadata = "sample-metadata.txt")

library(tidyverse)

#ALPHA diversity - Shannon Diversity

shannon <- read_qza("shannon_vector.qza")

shannon <- shannon$data %>% rownames_to_column("SampleID")

install.packages("ggplot")

gplots::venn(list(metadata=metadata$SampleID, shannon=shannon$SampleID))

metadata <- 
  metadata %>% 
  left_join(shannon)
head(metadata)

Meronthemealpha <-  theme(panel.border = element_rect(colour = "black", fill = NA),
                     panel.background = element_blank(),
                     panel.grid = element_blank(),
                     axis.text  = element_text(size = 12), 
                     axis.text.x = element_text(angle = 0),
                     axis.title = element_text(size = 20),
                     strip.background = element_blank(),
                     strip.text = element_text(size = 10),
                     axis.line.x = element_line(size = 0.5),
                     axis.line.y = element_line(size = 0.5),
                     legend.position = "right",
                     legend.text =element_text(size = 18))

metadata %>% 
  filter(!is.na(shannon_entropy)) %>%
  ggplot(aes(Type, shannon_entropy)) +
  geom_boxplot(position = "dodge")+
  geom_jitter(aes(fill = Depth), shape=21, width=0.2, height = 0, size = 2.5)+
  xlab(NULL)+
  ylab("Shannon Diversity")+
  Meronthemealpha+
  scale_fill_gradient(low = "white", high = "darkblue")

##Pielou's Evenness

evenness <- read_qza("evenness_vector.qza")

evenness <- evenness$data %>% rownames_to_column("SampleID")

metadata <- metadata %>% 
  left_join(evenness)


metadata %>% 
  #filter(!is.na(shannon_entropy)) %>%
  ggplot(aes(Type, pielou_evenness)) +
  geom_boxplot(position = "dodge")+
  geom_jitter(aes(fill = Depth), shape=21, width=0.2, height = 0, size = 2.5)+
  xlab(NULL)+
  ylab("Shannon Diversity")+
  Meronthemealpha+
  scale_fill_gradient(low = "white", high = "darkblue")




#Beta diverisy PCoA#Beta diverisy PCoAType

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
  ggplot(aes(PC1, PC2, color= Type,shape=Region)) +
  geom_point(alpha=0.5, size=5) +
  stat_ellipse() +
  xlab(paste("PC1: ", round(100*uwunifrac$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*uwunifrac$data$ProportionExplained[2]), "%")) +
  theme_bw() +
  #scale_size_continuous(name = "Shannon Diversity") +
  scale_shape_discrete(name ="Region") +
  scale_color_discrete(name="Type")

unifrac_modify %>% 
  ggplot(aes(x=PC1, y=PC2, color= Type, shape = Region, size = shannon_entropy))+
  geom_point(alpha=0.5) +
  xlab(paste("PC1: ", round(100*uwunifrac$data$ProportionExplained[1]), "%")) +
  ylab(paste("PC2: ", round(100*uwunifrac$data$ProportionExplained[2]), "%"))+
  #stat_ellipse()+
  scale_shape_discrete(name ="Region")+
  scale_size_continuous(name = "Shannon Diversity")+
  scale_color_discrete(name="Sample Type")+
  theme_bw()
  
  ggsave("PCoA.pdf", height = 4, width = 5, device= "pdf")



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



data <- read_csv("genus-level-composition.csv")


clado <- data %>% 
  group_by(Type, .drop = T) %>% 
  count(Cladocopium)



Fugacium, Clade I, Gerakladium, Durusdinium, Symbiodinium, Halluxium)

genus_order <- c("Cladocopium", "Durusdinium", "Fugacium", "Gerakladium",
                 "Clade I", "Halluxium", "Effrenium", "Woloszynskia",
                 "Biecheleria", "Ansanella", "Others")
ordered_genus <- factor(dat$Genus.x, levels = genus_order)


c("#771155","#AA4488", "#774411", "#114477", "#777711", "#DDDD77",
  "#4477AA", "#77AADD", "#117777","#44AAAA", "#77CCCC", "#117744",
  "#44AA77", "#88CCAA","#CC99BB", "#AA7744", "#DDAA77", "#771122",
  "#AAAA44", "#AA4455","#DD7788", )


"#00008F", "#774411", "#771122", "#CC99BB", "#117744", "#DDAA77", "#0000FF", "#DDDD77" , "#0070FF"




