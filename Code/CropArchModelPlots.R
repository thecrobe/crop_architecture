# Load necessary libraries
library(ggplot2)
library(ggsci)
library(brms)
library(reshape2)
library(tidybayes)

# Set the directory containing the models and load models
directory <- "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/models_base"
files <- list.files(directory, pattern = "\\.RDS$", full.names = TRUE)
models <- lapply(files, readRDS)

##### Branching Mechanism #####
ps_woody <- posterior_samples(models$woody_BranchingMechanism)[,1:3]
names(ps_woody)
# Get chain from intercept
ps_exp_woody <- exp(ps_woody) # exponentiating all chains get get chain 
ps_sums_woody <- rowSums(exp(ps_woody))+1 # summing them to get the denominator add the 1 for the baseline case
ps_sums_woody
df_woody<-ps_exp_woody/ps_sums_woody
#df_woody$WR2toCrop<-1/ps_sums_woody #return baseline to df
dfm_woody <- df_woody %>%  # organize for plotting
  select( -b_muBoth_Intercept)  %>% 
  rename(
    "Wild" = b_muWR2toWR1_Intercept,
    "Crop" = b_muWR2toCrop_Intercept) %>% 
  mutate(form = "woody") %>%
  mutate(Trait = "Branching Mechanism") %>%
  melt()

ps_herb <- posterior_samples(models$herb_BranchingMechanism)[,1:3]
names(ps_herb)
ps_exp_herb <- exp(ps_herb) # exponentiating all chains
ps_sums_herb <- rowSums(exp(ps_herb))+1 # summing them to get the denominator add the 1 for the baseline case
ps_sums_herb
df_herb<-ps_exp_herb/ps_sums_herb
#df_herb$WR2toCrop<-1/ps_sums_herb #return baseline to df
dfm_herb <- df_herb %>%  # organize for plotting
  select(-b_muBoth_Intercept)  %>% 
  rename(
    "Wild" = b_muWR2toWR1_Intercept,
    "Crop" = b_muWR2toCrop_Intercept) %>% 
  mutate(form = "herb") %>% 
  mutate(Trait = "Branching Mechanism") %>%
  melt()

##### Hypothesis testing 
dat1 <- as.data.frame(models$woody_BranchingMechanism)
hypothesis(dat1, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

dat2 <- as.data.frame(models$herb_BranchingMechanism)
hypothesis(dat2, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

# Plot
df_plot<-rbind(dfm_woody,dfm_herb) 
estimates_BranchingMechanism<-df_plot

phylo_plot<-ggplot(df_plot, aes(x=value*100,y=form,fill=variable)) +
  stat_slabinterval() +
  theme_bw() + 
  xlab("Probability of change (%)") + 
  ylab("Density") +
  ggsci::scale_fill_aaas(name="Direction") + 
  ggsci::scale_color_aaas(name="Direction") + 
  ggtitle(label = "Branching Mechanism",subtitle = "P woody = 84% | P herb = 89%")
phylo_plot
ggsave("/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/Figures/BranchingMechanism_base.png",scale = 0.75)


##### Branching Position ######
ps_woody <- posterior_samples(models$woody_BranchingPosition)[,1:3]
names(ps_woody)
# Get chain from intercept
ps_exp_woody <- exp(ps_woody) # exponentiating all chains get get chain 
ps_sums_woody <- rowSums(exp(ps_woody)) + 1 # summing them to get the denominator add the 1 for the baseline case
ps_sums_woody
df_woody <- ps_exp_woody / ps_sums_woody
dfm_woody <- df_woody %>%
  dplyr::select(-b_muBoth_Intercept) %>% 
  dplyr::rename(
    "Wild" = b_muWR2toWR1_Intercept,
    "Crop" = b_muWR2toCrop_Intercept
  ) %>% 
  dplyr::mutate(form = "woody", Trait = "Branching Position") %>%
  tidyr::pivot_longer(cols = c(Wild, Crop), names_to = "variable", values_to = "value")

ps_herb <- posterior_samples(models$herb_BranchingPosition)[,1:3]
names(ps_herb)
ps_exp_herb <- exp(ps_herb) # exponentiating all chains
ps_sums_herb <- rowSums(exp(ps_herb)) + 1 # summing them to get the denominator add the 1 for the baseline case
ps_sums_herb
df_herb <- ps_exp_herb / ps_sums_herb
dfm_herb <- df_herb %>%
  dplyr::select(-b_muBoth_Intercept) %>%
  dplyr::rename(
    "Wild" = b_muWR2toWR1_Intercept,
    "Crop" = b_muWR2toCrop_Intercept
  ) %>%
  dplyr::mutate(form = "herb", Trait = "Branching Position") %>%
  tidyr::pivot_longer(cols = c(Wild, Crop), names_to = "variable", values_to = "value")

# Hypothesis testing 
dat1 <- as.data.frame(models$woody_BranchingPosition)
hypothesis(dat1, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

dat2 <- as.data.frame(models$herb_BranchingPosition)
hypothesis(dat2, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

# Combine and plot
df_plot <- dplyr::bind_rows(dfm_woody, dfm_herb)
phylo_plot <- ggplot(df_plot, aes(x = value * 100, y = form, fill = variable)) +
  stat_slabinterval() +
  theme_bw() + 
  xlab("Probability of change (%)") + 
  ylab("Density") +
  ggsci::scale_fill_aaas(name = "Direction") + 
  ggsci::scale_color_aaas(name = "Direction") + 
  ggtitle(label = "Branching Position", subtitle = "P woody = 99% | P herb = 98%")

print(phylo_plot)

# Save the plot
ggsave("/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/Figures/BranchingPosition_base.png",scale = 0.75)

##### Branching Type ######
ps_woody <- posterior_samples(models$woody_BranchingType)[, 1:3]
names(ps_woody)
# Get chain from intercept
ps_exp_woody <- exp(ps_woody)  # exponentiating all chains to get the chain
ps_sums_woody <- rowSums(exp(ps_woody)) + 1  # summing them to get the denominator, add the 1 for the baseline case
ps_sums_woody
df_woody <- ps_exp_woody / ps_sums_woody
dfm_woody <- df_woody %>%
  select(-b_muBoth_Intercept) %>% 
  rename(
    "Wild" = b_muWR2toWR1_Intercept,
    "Crop" = b_muWR2toCrop_Intercept
  ) %>% 
  mutate(form = "woody", Trait = "Branching Type") %>%
  pivot_longer(cols = c(Wild, Crop), names_to = "variable", values_to = "value")

ps_herb <- posterior_samples(models$herb_BranchingType)[, 1:3]
names(ps_herb)
ps_exp_herb <- exp(ps_herb)  # exponentiating all chains
ps_sums_herb <- rowSums(exp(ps_herb)) + 1  # summing them to get the denominator, add the 1 for the baseline case
ps_sums_herb
df_herb <- ps_exp_herb / ps_sums_herb
dfm_herb <- df_herb %>%
  select(-b_muBoth_Intercept) %>%
  rename(
    "Wild" = b_muWR2toWR1_Intercept,
    "Crop" = b_muWR2toCrop_Intercept
  ) %>%
  mutate(form = "herb", Trait = "Branching Type") %>%
  pivot_longer(cols = c(Wild, Crop), names_to = "variable", values_to = "value")

# Hypothesis testing 
dat1 <- as.data.frame(models$woody_BranchingType)
hypothesis(dat1, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

dat2 <- as.data.frame(models$herb_BranchingType)
hypothesis(dat2, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

# Combine and plot
df_plot <- bind_rows(dfm_woody, dfm_herb)
phylo_plot <- ggplot(df_plot, aes(x = value * 100, y = form, fill = variable)) +
  stat_slabinterval() +
  theme_bw() + 
  xlab("Probability of change (%)") + 
  ylab("Density") +
  ggsci::scale_fill_aaas(name = "Direction") + 
  ggsci::scale_color_aaas(name = "Direction") + 
  ggtitle(label = "Branching Type", subtitle = "P woody = 90% | P herb = 96%")

print(phylo_plot)

# Save the plot
ggsave("/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/Figures/BranchingType_base.png",scale = 0.75)

##### Flowering Axis ######
ps_woody <- posterior_samples(models$woody_FloweringAxis)[, 1:3]
names(ps_woody)
# Get chain from intercept
ps_exp_woody <- exp(ps_woody)  # exponentiating all chains to get the chain
ps_sums_woody <- rowSums(exp(ps_woody)) + 1  # summing them to get the denominator, add the 1 for the baseline case
ps_sums_woody
df_woody <- ps_exp_woody / ps_sums_woody
dfm_woody <- df_woody %>%
  select(-b_muBoth_Intercept) %>% 
  rename(
    "Wild" = b_muWR2toWR1_Intercept,
    "Crop" = b_muWR2toCrop_Intercept
  ) %>% 
  mutate(form = "woody", Trait = "Flowering Axis") %>%
  pivot_longer(cols = c(Wild, Crop), names_to = "variable", values_to = "value")

ps_herb <- posterior_samples(models$herb_FloweringAxis)[, 1:3]
names(ps_herb)
ps_exp_herb <- exp(ps_herb)  # exponentiating all chains
ps_sums_herb <- rowSums(exp(ps_herb)) + 1  # summing them to get the denominator, add the 1 for the baseline case
ps_sums_herb
df_herb <- ps_exp_herb / ps_sums_herb
dfm_herb <- df_herb %>%
  select(-b_muBoth_Intercept) %>%
  rename(
    "Wild" = b_muWR2toWR1_Intercept,
    "Crop" = b_muWR2toCrop_Intercept
  ) %>%
  mutate(form = "herb", Trait = "Flowering Axis") %>%
  pivot_longer(cols = c(Wild, Crop), names_to = "variable", values_to = "value")

# Hypothesis testing 
dat1 <- as.data.frame(models$woody_FloweringAxis)
hypothesis(dat1, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

dat2 <- as.data.frame(models$herb_FloweringAxis)
hypothesis(dat2, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

# Combine and plot
df_plot <- bind_rows(dfm_woody, dfm_herb)
phylo_plot <- ggplot(df_plot, aes(x = value * 100, y = form, fill = variable)) +
  stat_slabinterval() +
  theme_bw() + 
  xlab("Probability of change (%)") + 
  ylab("Density") +
  ggsci::scale_fill_aaas(name = "Direction") + 
  ggsci::scale_color_aaas(name = "Direction") + 
  ggtitle(label = "Flowering Axis", subtitle = "P woody = 79% | P herb = 90%")

print(phylo_plot)

# Save the plot
ggsave("/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/Figures/FloweringAxis_base.png",scale = 0.75)

##### Growth Direction #####
ps_woody <- posterior_samples(models$woody_GrowthDirection)[, 1:3]
names(ps_woody)
# Get chain from intercept
ps_exp_woody <- exp(ps_woody)  # exponentiating all chains to get the chain
ps_sums_woody <- rowSums(exp(ps_woody)) + 1  # summing them to get the denominator, add the 1 for the baseline case
ps_sums_woody
df_woody <- ps_exp_woody / ps_sums_woody
dfm_woody <- df_woody %>%
  select(-b_muBoth_Intercept) %>% 
  rename(
    "Wild" = b_muWR2toWR1_Intercept,
    "Crop" = b_muWR2toCrop_Intercept
  ) %>% 
  mutate(form = "woody", Trait = "Growth Direction") %>%
  pivot_longer(cols = c(Wild, Crop), names_to = "variable", values_to = "value")

ps_herb <- posterior_samples(models$herb_GrowthDirection)[, 1:3]
names(ps_herb)
ps_exp_herb <- exp(ps_herb)  # exponentiating all chains
ps_sums_herb <- rowSums(exp(ps_herb)) + 1  # summing them to get the denominator, add the 1 for the baseline case
ps_sums_herb
df_herb <- ps_exp_herb / ps_sums_herb
dfm_herb <- df_herb %>%
  select(-b_muBoth_Intercept) %>%
  rename(
    "Wild" = b_muWR2toWR1_Intercept,
    "Crop" = b_muWR2toCrop_Intercept
  ) %>%
  mutate(form = "herb", Trait = "Growth Direction") %>%
  pivot_longer(cols = c(Wild, Crop), names_to = "variable", values_to = "value")

# Hypothesis testing 
dat1 <- as.data.frame(models$woody_GrowthDirection)
hypothesis(dat1, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

dat2 <- as.data.frame(models$herb_GrowthDirection)
hypothesis(dat2, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

# Combine and plot
df_plot <- bind_rows(dfm_woody, dfm_herb)
phylo_plot <- ggplot(df_plot, aes(x = value * 100, y = form, fill = variable)) +
  stat_slabinterval() +
  theme_bw() + 
  xlab("Probability of change (%)") + 
  ylab("Density") +
  ggsci::scale_fill_aaas(name = "Direction") + 
  ggsci::scale_color_aaas(name = "Direction") + 
  ggtitle(label = "Growth Direction", subtitle = "P woody = 84% | P herb = 79%")

print(phylo_plot)

# Save the plot
ggsave("/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/Figures/GrowthDirection_base.png",scale = 0.75)

##### Meristem Function ######
ps_woody <- posterior_samples(models$woody_MeristemFunction)[, 1:3]
names(ps_woody)
# Get chain from intercept
ps_exp_woody <- exp(ps_woody)  # exponentiating all chains to get the chain
ps_sums_woody <- rowSums(exp(ps_woody)) + 1  # summing them to get the denominator, add the 1 for the baseline case
ps_sums_woody
df_woody <- ps_exp_woody / ps_sums_woody
dfm_woody <- df_woody %>%
  select(-b_muBoth_Intercept) %>% 
  rename(
    "Wild" = b_muWR2toWR1_Intercept,
    "Crop" = b_muWR2toCrop_Intercept
  ) %>% 
  mutate(form = "woody", Trait = "Meristem Function") %>%
  pivot_longer(cols = c(Wild, Crop), names_to = "variable", values_to = "value")

ps_herb <- posterior_samples(models$herb_MeristemFunction)[, 1:3]
names(ps_herb)
ps_exp_herb <- exp(ps_herb)  # exponentiating all chains
ps_sums_herb <- rowSums(exp(ps_herb)) + 1  # summing them to get the denominator, add the 1 for the baseline case
ps_sums_herb
df_herb <- ps_exp_herb / ps_sums_herb
dfm_herb <- df_herb %>%
  select(-b_muBoth_Intercept) %>%
  rename(
    "Wild" = b_muWR2toWR1_Intercept,
    "Crop" = b_muWR2toCrop_Intercept
  ) %>%
  mutate(form = "herb", Trait = "Meristem Function") %>%
  pivot_longer(cols = c(Wild, Crop), names_to = "variable", values_to = "value")

# Hypothesis testing 
dat1 <- as.data.frame(models$woody_MeristemFunction)
hypothesis(dat1, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

dat2 <- as.data.frame(models$herb_MeristemFunction)
hypothesis(dat2, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

# Combine and plot
df_plot <- bind_rows(dfm_woody, dfm_herb)
phylo_plot <- ggplot(df_plot, aes(x = value * 100, y = form, fill = variable)) +
  stat_slabinterval() +
  theme_bw() + 
  xlab("Probability of change (%)") + 
  ylab("Density") +
  ggsci::scale_fill_aaas(name = "Direction") + 
  ggsci::scale_color_aaas(name = "Direction") + 
  ggtitle(label = "Meristem Function", subtitle = "P woody = 55% | P herb = 97%")

print(phylo_plot)

# Save the plot
ggsave("/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/Figures/MeristemFunction_base.png",scale = 0.75)

##### Phyllotaxis ######
ps_woody <- posterior_samples(models$woody_Phyllotaxis)[, 1:3]
names(ps_woody)
# Get chain from intercept
ps_exp_woody <- exp(ps_woody)  # exponentiating all chains to get the chain
ps_sums_woody <- rowSums(exp(ps_woody)) + 1  # summing them to get the denominator, add the 1 for the baseline case
ps_sums_woody
df_woody <- ps_exp_woody / ps_sums_woody
dfm_woody <- df_woody %>%
  select(-b_muBoth_Intercept) %>% 
  rename(
    "Wild" = b_muWR2toWR1_Intercept,
    "Crop" = b_muWR2toCrop_Intercept
  ) %>% 
  mutate(form = "woody", Trait = "Phyllotaxis") %>%
  pivot_longer(cols = c(Wild, Crop), names_to = "variable", values_to = "value")

ps_herb <- posterior_samples(models$herb_Phyllotaxis)[, 1:3]
names(ps_herb)
ps_exp_herb <- exp(ps_herb)  # exponentiating all chains
ps_sums_herb <- rowSums(exp(ps_herb)) + 1  # summing them to get the denominator, add the 1 for the baseline case
ps_sums_herb
df_herb <- ps_exp_herb / ps_sums_herb
dfm_herb <- df_herb %>%
  select(-b_muBoth_Intercept) %>%
  rename(
    "Wild" = b_muWR2toWR1_Intercept,
    "Crop" = b_muWR2toCrop_Intercept
  ) %>%
  mutate(form = "herb", Trait = "Phyllotaxis") %>%
  pivot_longer(cols = c(Wild, Crop), names_to = "variable", values_to = "value")

# Hypothesis testing 
dat1 <- as.data.frame(models$woody_Phyllotaxis)
hypothesis(dat1, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

dat2 <- as.data.frame(models$herb_Phyllotaxis)
hypothesis(dat2, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

# Combine and plot
df_plot <- bind_rows(dfm_woody, dfm_herb)
phylo_plot <- ggplot(df_plot, aes(x = value * 100, y = form, fill = variable)) +
  stat_slabinterval() +
  theme_bw() + 
  xlab("Probability of change (%)") + 
  ylab("Density") +
  ggsci::scale_fill_aaas(name = "Direction") + 
  ggsci::scale_color_aaas(name = "Direction") + 
  ggtitle(label = "Phyllotaxis", subtitle = "P woody = 98% | P herb = 96%")

print(phylo_plot)

# Save the plot
ggsave("/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/Figures/Phyllotaxis_base.png",scale = 0.75)

##### Short Shoots #####
ps_woody <- posterior_samples(models$woody_ShortShoots)[, 1:3]
names(ps_woody)
# Get chain from intercept
ps_exp_woody <- exp(ps_woody)  # exponentiating all chains to get the chain
ps_sums_woody <- rowSums(exp(ps_woody)) + 1  # summing them to get the denominator, add the 1 for the baseline case
ps_sums_woody
df_woody <- ps_exp_woody / ps_sums_woody
dfm_woody <- df_woody %>%
  select(-b_muBoth_Intercept) %>% 
  rename(
    "Wild" = b_muWR2toWR1_Intercept,
    "Crop" = b_muWR2toCrop_Intercept
  ) %>% 
  mutate(form = "woody", Trait = "Short Shoots") %>%
  pivot_longer(cols = c(Wild, Crop), names_to = "variable", values_to = "value")

ps_herb <- posterior_samples(models$herb_ShortShoots)[, 1:3]
names(ps_herb)
ps_exp_herb <- exp(ps_herb)  # exponentiating all chains
ps_sums_herb <- rowSums(exp(ps_herb)) + 1  # summing them to get the denominator, add the 1 for the baseline case
ps_sums_herb
df_herb <- ps_exp_herb / ps_sums_herb
dfm_herb <- df_herb %>%
  select(-b_muBoth_Intercept) %>%
  rename(
    "Wild" = b_muWR2toWR1_Intercept,
    "Crop" = b_muWR2toCrop_Intercept
  ) %>%
  mutate(form = "herb", Trait = "Short Shoots") %>%
  pivot_longer(cols = c(Wild, Crop), names_to = "variable", values_to = "value")

# Hypothesis testing 
dat1 <- as.data.frame(models$woody_ShortShoots)
hypothesis(dat1, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

dat2 <- as.data.frame(models$herb_ShortShoots)
hypothesis(dat2, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

# Combine and plot
df_plot <- bind_rows(dfm_woody, dfm_herb)
phylo_plot <- ggplot(df_plot, aes(x = value * 100, y = form, fill = variable)) +
  stat_slabinterval() +
  theme_bw() + 
  xlab("Probability of change (%)") + 
  ylab("Density") +
  ggsci::scale_fill_aaas(name = "Direction") + 
  ggsci::scale_color_aaas(name = "Direction") + 
  ggtitle(label = "Short Shoots", subtitle = "P woody = 75% | P herb = 81%")

print(phylo_plot)

# Save the plot
ggsave("/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/Figures/ShortShoots_base.png",scale = 0.75)


