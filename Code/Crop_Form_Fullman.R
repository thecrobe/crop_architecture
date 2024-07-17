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

##### Import traits #####
form<-read.csv("/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/Crop_Form_Fullman.csv", header=T) %>% na.omit()
trait<-read.csv("/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/LongFormat_croparchitecture_V2.csv", header=T)
transitioned<-trait %>% filter (WR2toCrop == 1)

form <- form[!form$Crop %in% transitioned$Tree, ] %>%
  filter(Pathway_Fuller_Cleaned != "Ruderal") %>% 
  filter(Pathway_Fuller_Cleaned != "Ruderal/Ornamental") %>% 
  filter(Pathway_Fuller_Cleaned != "Ruderal/Tuber")
                                                               
form<-na.omit(form)
dim(form)
unique(form$Pathway_Fuller_Cleaned)

# Create a contingency table
contingency_table <- table(form$Herbaceous.woody, form$Pathway_Fuller_Cleaned)

# Perform chi-squared test of independence
chi_square_result <- chisq.test(contingency_table,simulate.p.value = TRUE)

# Print the test result
print(chi_square_result)

plot_data <- form %>%
  group_by(Pathway_Fuller_Cleaned, Herbaceous.woody) %>%
  summarise(Count = n(), .groups = 'drop')

# Plot
plot<-ggplot(plot_data, aes(x = reorder(Pathway_Fuller_Cleaned,-Count), y = Count, fill = Herbaceous.woody)) +
  geom_bar(stat = "identity", position = ) +
  scale_fill_manual(values = c("Herbaceous" = "green", "Woody" = "purple")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotating x-axis labels
  labs(title = "Crop target use",
       x = "",
       y = "Count",
       fill = "Type")

ggsave(plot = plot,filename = "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/Figures/TargetUse.png",scale = 0.7)
