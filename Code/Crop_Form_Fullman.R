##### Libraries #####
library(ggplot2) #plotting
library(dplyr)#data wrangling 
library(reshape2) #data wrangling
library(brms) #modeling
library(ggbreak) #plotting
library(ggsci) #plotting
library(ggpubr) #plots
library(tidybayes) #plotting
library(ape) #phylogeny
library(TreeTools) #phylogeny
library(rcartocolor)
library(taxize)

##### Import traits #####
form<-read.csv("~/Dropbox/other_projects/crop_architecture/crop_architecture/Data/Crop_Form_Fullman.csv", header=T) %>% na.omit()
form$Crop <- gsub("_$", "", form$Crop)
form$Crop <- gsub("__$", "", form$Crop)
form$Crop <- gsub("_x_", "_", form$Crop)
form$Crop <- gsub("_×_", "_", form$Crop)
form$Crop <- gsub("__", "_", form$Crop)
trait<-read.csv("~/Dropbox/other_projects/crop_architecture/crop_architecture/Data/20250211_LongFormat_croparchitecture_V3.csv", header=T)
trait<-trait %>%
  mutate(Herbaceous.woody = case_when(
    Herbaceous.woody == "Graminoid" ~ "Herbaceous",
    TRUE ~ Herbaceous.woody  # keeps all other values the same
  ))
trait <- trait[which(trait$wild.cultivated == "Cultivated"),]
transitioned<-trait %>% filter (WR2toCrop == 1 & Both == 0)

# check overlap
# trait[!trait$Species.ssp.var %in% form$Crop,]$Species.ssp.var

missing <- unique(trait[!trait$Species.ssp.var %in% form$Crop,]$Species.ssp.var)

missing_df <- data.frame(og_name = missing)
missing_df$name_str <- sapply(strsplit(missing_df$og_name, "_"), FUN = function(x) { paste(x, collapse = " ")})
missing_df$id <- rep("", nrow(missing_df))
missing_df$is_synonym <- rep(FALSE, nrow(missing_df))
missing_df$accepted_name <- rep("", nrow(missing_df))
for (i in 1:nrow(missing_df)) {
  missing_df$id[i] <- get_pow(missing_df$name_str[i], ask = FALSE)[1]
  Sys.sleep(0.2)
}

for (i in 1:nrow(missing_df)) {
  if (is.na(missing_df$id[i])) next
  result <- pow_lookup(id = missing_df$id[i])
  missing_df$is_synonym[i] <- result$meta$synonym
  if (result$meta$synonym) missing_df$accepted_name[i] <- result$meta$accepted$name
  Sys.sleep(0.2)
}

missing_df$accepted_name <- gsub(" ", "_", missing_df$accepted_name)

missing_df[which(missing_df$is_synonym & missing_df$accepted_name %in% form$Crop),]
#                   og_name              name_str                                 id is_synonym            accepted_name
# 42    Chionachne_gigantea   Chionachne gigantea urn:lsid:ipni.org:names:20010084-1       TRUE        Polytoca_gigantea
# 71      Hyptis_suaveolens     Hyptis suaveolens   urn:lsid:ipni.org:names:448418-1       TRUE  Mesosphaerum_suaveolens
# 106         Pisum_sativum         Pisum sativum urn:lsid:ipni.org:names:60454055-2       TRUE       Lathyrus_oleraceus
# 112 Pueraria_phaseoloides Pueraria phaseoloides   urn:lsid:ipni.org:names:516712-1       TRUE Neustanthus_phaseoloides
# 113    Quararibea_cordata    Quararibea cordata   urn:lsid:ipni.org:names:215373-2       TRUE          Matisia_cordata

# account for synonymy where possible
trait[which(trait$Species.ssp.var == "Chionachne_gigantea"),]$Species.ssp.var <- "Polytoca_gigantea"
trait[which(trait$Species.ssp.var == "Hyptis_suaveolens"),]$Species.ssp.var <- "Mesosphaerum_suaveolens"
trait[which(trait$Species.ssp.var == "Pisum_sativum"),]$Species.ssp.var <- "Lathyrus_oleraceus"
trait[which(trait$Species.ssp.var == "Pueraria_phaseoloides"),]$Species.ssp.var <- "Neustanthus_phaseoloides"
trait[which(trait$Species.ssp.var == "Quararibea_cordata"),]$Species.ssp.var <- "Matisia_cordata"

trait[which(trait$Species.ssp.var == "Vigna_unguiculata_Cowpea"),]$Species.ssp.var <- "Vigna_unguiculata"
trait[which(trait$Species.ssp.var == "Vigna_unguiculata_aparagus"),]$Species.ssp.var <- "Vigna_unguiculata_subsp._unguiculata"

# remove ruderals and low sample sizes
form <- form[form$Crop %in% trait$Species.ssp.var, ] %>%
  filter(Pathway_Fuller_Cleaned != "Fruit (not tree)") %>%
  filter(Pathway_Fuller_Cleaned != "Fodder") %>%
  filter(Pathway_Fuller_Cleaned != "Ornamental")
                                                               
form<-na.omit(form)
dim(form)
unique(form$Pathway_Fuller_Cleaned)
table(form$Pathway_Fuller_Cleaned)

# group several and 
form[which(form$Pathway_Fuller_Cleaned == "Fruit tree / Ecosystem engineering"),]$Pathway_Fuller_Cleaned <- "Fruit tree"
form[which(form$Pathway_Fuller_Cleaned == "Ruderal/Tuber"),]$Pathway_Fuller_Cleaned <- "Ruderal"
form[which(form$Pathway_Fuller_Cleaned == "Ruderal/Ornamental"),]$Pathway_Fuller_Cleaned <- "Ruderal"
form[which(form$Pathway_Fuller_Cleaned == "Tuber (but young shoots are eaten)"),]$Pathway_Fuller_Cleaned <- "Tuber"


table(form$Pathway_Fuller_Cleaned)

# get traits dataset for species in form database
trait <- trait[trait$Species.ssp.var %in% form$Crop,] %>%
  filter(Trait != "Branching") %>%
  filter(Trait != "Reiteration") %>%
  filter(Trait != "BasalBranching")

trait$form <- form[match(trait$Species.ssp.var, form$Crop),]$Pathway_Fuller_Cleaned

trait$crop_transition <- ifelse(trait$WR2toCrop == 1 & trait$Both == 0, TRUE, FALSE)

tab <- table(trait$Trait, trait$form, trait$crop_transition)

df <- as.data.frame(tab)
colnames(df) <- c("Trait", "Pathway", "Transitions", "Freq")

# df$Trait <- factor(df$Trait, levels = rev(c("BranchingPosition", "GrowthDirection", "BranchingMechanism",
#                                         "Flowering", "Phyllotaxis", "MeristemFunctioning", "BranchingType",
#                                         "ShortShoots")))

df$Trait <- factor(df$Trait, levels = c("BranchingPosition", "GrowthDirection", "BranchingMechanism",
                                            "Flowering", "Phyllotaxis", "MeristemFunctioning", "BranchingType",
                                            "ShortShoots"))

xlabs <- c("Branching Position", "Growth Direction", "Branching Mechanism",
           "Flowering Axis", "Phyllotaxis", "Meristem Functioning", "Branching Type",
           "Short Shoots")

p <- ggplot(data=df, aes(x=Transitions, fill=Pathway, y=Freq)) + 
  geom_bar(position = "fill", stat = "identity") + 
  geom_text(data=subset(df, Freq != 0, ), aes(label = Freq), stat='identity',
            position = position_fill(vjust = 0.5), col="white") +
  facet_grid(Trait ~ ., scales = "free", space="free", switch = "x") +
  coord_flip() +
  # scale_x_discrete(labels = function(x) stringr::str_wrap(rev(xlabs), width = 10)) +
  scale_fill_carto_d(type = "qualitative", palette = "Safe", direction = -1,
                     name = "Domestication\npathway",
                     labels = c("Ecosystem\nengineering",
                                "Fibre",
                                "Fruit tree",
                                "Grain",
                                "Ruderal",
                                "Segetal",
                                "Tuber")) + 
  ylab("Proportion") +
  xlab("Transitions") +
  theme_bw()

ggsave("~/Dropbox/other_projects/crop_architecture/figures/FigX_pathway_by_trait_v4.svg", p, width = 7, height = 9)

tests <- data.frame(trait=levels(df$Trait),
                    p=rep(0., 8))
for (i in levels(df$Trait)) {
  tmp <- df[which(df$Trait == i),]
  ss <- xtabs(Freq ~ Transitions + Pathway, data = tmp)
  res <- chisq.test(ss, simulate.p.value = T)
  tests[which(tests$trait == i),]$p <- res$p.value
}

# residuals for GrowthDirection
# Pathway
# Transitions     Ecosystem engineering       Fibre  Fruit tree       Grain     Ruderal     Segetal       Tuber
#           FALSE            0.15882071 -0.53473207  0.18021266 -0.14755616  0.22208556 -0.63378492  0.09788927
#           TRUE            -0.84935677  2.85969191 -0.96375869  0.78911514 -1.18769067  3.38941632 -0.52350172

# residuals for phyllotaxis
#            Pathway
# Transitions Ecosystem engineering       Fibre  Fruit tree       Grain     Ruderal     Segetal       Tuber
#       FALSE           -0.13675531 -0.56241759 -0.09933903  0.10707019  0.19760568  0.06532923  0.08492466
#       TRUE             0.86611695  3.56197806  0.62914722 -0.67811123 -1.25150262 -0.41375179 -0.53785620

# # to break up by transition terminal state, we need to compress states as in architectural convergence
# # these are either compressible OR have 1 observation
# transitioned[which(transitioned$VALUE == "Axis 1 Orthotropic"),]$VALUE <- "Orthotropic"
# # transitioned <- transitioned[-which(transitioned$VALUE == "Orthotropic and Plagiotropic"),]
# # 
# transitioned[which(transitioned$VALUE == "Axis 1 Monopodial"),]$VALUE <- "Monopodial"
# # transitioned <- transitioned[-which(transitioned$VALUE == "Monopodial or Sympodial"),]
# # 
# transitioned[which(transitioned$VALUE == "Axis 1 Indeterminate"),]$VALUE <- "Indeterminate"
# # transitioned <- transitioned[-which(transitioned$VALUE == "Determinate or Indeterminate"),]
# 
# 
# transitioned$Trait <- factor(transitioned$Trait, levels = c("BranchingPosition", "GrowthDirection", "BranchingMechanism",
#                                                                 "Flowering", "Phyllotaxis", "MeristemFunctioning", "BranchingType",
#                                                                 "ShortShoots"))
# 
# xlabs <- c("Branching Position", "Growth Direction", "Branching Mechanism",
#            "Flowering Axis", "Phyllotaxis", "Meristem Functioning", "Branching Type",
#            "Short Shoots")
# 
# 
# p <- ggplot(data=transitioned, aes(x=VALUE, fill=form)) +
#   geom_bar(position = "fill", stat = "count", width=0.8) +
#   geom_text(aes(label=after_stat(count)), stat="count", position = position_fill(vjust = 0.5), col="white") +
#   scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) +
#   coord_flip() +
#   facet_grid(Trait ~ ., scales = "free", space="free", switch = "x") +
#   scale_fill_carto_d(type = "qualitative", palette = "Safe", direction = -1,
#                      name = "Domestication\npathway",
#                      labels = c("Ecosystem\nengineering",
#                                 "Fruit tree",
#                                 "Grain",
#                                 "Ruderal",
#                                 "Segetal",
#                                 "Tuber")) +
#   xlab("Terminal state") +
#   ylab("Proportion of transitions") +
#   theme_bw()
# 
# ggsave("../figures/FigX_pathway_by_trait_value_v2.svg", p, width = 7, height = 8)
# 
# # # Create a contingency table
# # contingency_table <- table(form$Herbaceous.woody, form$Pathway_Fuller_Cleaned)
# # 
# # # Perform chi-squared test of independence
# # chi_square_result <- chisq.test(contingency_table,simulate.p.value = TRUE)
# # 
# # # Print the test result
# # print(chi_square_result)
# # 
# # plot_data <- form %>%
# #   group_by(Pathway_Fuller_Cleaned, Herbaceous.woody) %>%
# #   summarise(Count = n(), .groups = 'drop')
# # 
# # # Plot
# # plot<-ggplot(plot_data, aes(x = reorder(Pathway_Fuller_Cleaned,-Count), y = Count, fill = Herbaceous.woody)) +
# #   geom_bar(stat = "identity", position = ) +
# #   scale_fill_manual(values = c("Herbaceous" = "green", "Woody" = "purple")) +
# #   theme_minimal() +
# #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotating x-axis labels
# #   labs(title = "Crop target use",
# #        x = "",
# #        y = "Count",
# #        fill = "Type")
# # 
# # ggsave(plot = plot,filename = "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/Figures/TargetUse.png",scale = 0.7)
