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
setwd("~/Dropbox/other_projects/crop_architecture/crop_architecture/")
trait<-read.csv("Data/20250211_LongFormat_croparchitecture_V3.csv", header=T)

c<-trait %>% filter(wild.cultivated =="Cultivated")

trait$WR2toCrop<-as.numeric(trait$WR2toCrop)
trait$WR2toWR1<-as.numeric(trait$WR2toWR1)
trait$Both<-as.numeric(trait$Both)
trait$Neither<-as.numeric(trait$Neither)

trait<-trait %>%
  mutate(Herbaceous.woody = case_when(
    Herbaceous.woody == "Graminoid" ~ "Herbaceous",
    TRUE ~ Herbaceous.woody  # keeps all other values the same
  ))
# IS TIP OR BASE THE RIGHT TREE TO USE? MODEL COEFFICIENTS ARE SIMILAR EITHER WAY, BUT NOT SURE WHAT ORIGINAL ANALYSIS USED
tree<-ape::read.tree("~/Dropbox/other_projects/crop_architecture/crop_architecture/Output/tree_pgls_base_clean.tre")

###### GrowthDirection Model #####
GrowthDirection <- trait %>% filter(Trait == "GrowthDirection", VALUE != "NoData")
GrowthDirection$Tree <- GrowthDirection$Species.ssp.var
GrowthDirection_drop<-trait %>% filter(Trait == "Growthdirection", VALUE == "NoData")
GrowthDirection_tree<-drop.tip(tree, GrowthDirection_drop$Tree)
labels<-data.frame(TipLabels(GrowthDirection_tree)) 
labels$Tree<-labels$TipLabels.GrowthDirection_tree.
d<-dplyr::inner_join(labels, GrowthDirection)
data<-d %>% filter(d$wild.cultivated == "Cultivated")
#Prepare tree for model 
A <- ape::vcv.phylo(GrowthDirection_tree)
B <-diag(diag(A + 0.00000000000000001)) #add to diagonal to make matrix positive definite.
row.names(B)<-row.names(A)
#Select only the cultivated crops, this is where the binary coding of changes lives in the data. Also seperate out and start with woody crops
woody<-data %>% filter(wild.cultivated == "Cultivated" & Herbaceous.woody == "Woody")
woody$Neither
#set response for model
woody$size <- with(woody,WR2toCrop+WR2toWR1+Both+Neither) 
woody$y <- with(woody, cbind(Neither,Both,WR2toCrop,WR2toWR1)) 
#Run model

woody_GrowthDirection <- brm(bf(y |trials(size) ~ 1 + (1|gr(Tree,cov=B))),
                             family = multinomial(),
                             data = data.frame(woody),
                             data2 = list(B=B),
                             chains = 4,
                             thin = 5,
                             iter = 3000,
                             warmup = 500,
                             sample_prior = TRUE,
                             cores = 4,
                             control = list(adapt_delta=0.85))
print(woody_GrowthDirection)
saveRDS(woody_GrowthDirection, file = "~/Dropbox/other_projects/crop_architecture/woody_GrowthDirection.RDS")
# woody_GrowthDirection <- readRDS("../woody_GrowthDirection.RDS")

#Repeat steps above for herbaceous plans
herb<-data %>% filter(wild.cultivated == "Cultivated" & Herbaceous.woody == "Herbaceous")
herb$size <- with(herb,WR2toCrop+WR2toWR1+Both+Neither)
herb$y <- with(herb, cbind(Neither,Both,WR2toCrop,WR2toWR1))

herb_GrowthDirection <- brm(bf(y |trials(size) ~ 1 + (1|gr(Tree,cov=B))),
                            family = multinomial(),
                            data = data.frame(herb),
                            data2 = list(B=B),
                            chains = 4,
                            thin = 5,
                            iter = 3000,
                            warmup = 500,
                            sample_prior = TRUE,
                            cores = 4)
print(herb_GrowthDirection)
saveRDS(herb_GrowthDirection, file = "~/Dropbox/other_projects/crop_architecture/herb_GrowthDirection.RDS")
# herb_GrowthDirection <- readRDS("../herb_GrowthDirection.RDS")

# Extract posterior distributions from the Models_tip,
# ps_woody <- posterior_samples(woody_GrowthDirection)[,1:3]
ps_woody <- as_draws_df(woody_GrowthDirection,
                        variable = c("b_muBoth_Intercept", "b_muWR2toCrop_Intercept", "b_muWR2toWR1_Intercept"))[,1:3]
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
  mutate(Trait = "Growth Direction") %>%
  melt()

# ps_herb <- posterior_samples(herb_GrowthDirection)[,1:3]
ps_herb <- as_draws_df(herb_GrowthDirection,
                       variable = c("b_muBoth_Intercept", "b_muWR2toCrop_Intercept", "b_muWR2toWR1_Intercept"))[,1:3]
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
  mutate(Trait = "Growth Direction") %>%
  melt()

##### Hypothesis testing 
dat1 <- as.data.frame(woody_GrowthDirection)
h_woody_GrowthDirection <- hypothesis(dat1, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")
dat2 <- as.data.frame(herb_GrowthDirection)
h_herb_GrowthDirection <- hypothesis(dat2, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

# Plot
df_plot<-rbind(dfm_woody,dfm_herb) 
estimates_GrowthDirection<-df_plot

phylo_plot<-ggplot(df_plot, aes(x=value*100,y=form,fill=variable)) +
  stat_slabinterval() +
  theme_bw() + 
  xlab("Probability of change (%)") + 
  ylab("Density") +
  ggsci::scale_fill_aaas(name="Direction") + 
  ggsci::scale_color_aaas(name="Direction")
phylo_plot



###### Phyllotaxis Model #####
Phyllotaxis<-trait %>% filter(Trait == "Phyllotaxis" & VALUE != "NoData")
Phyllotaxis$Tree <- Phyllotaxis$Species.ssp.var
Phyllotaxis_drop<-trait %>% filter(Trait == "Phyllotaxis", VALUE == "NoData")
Phyllotaxis_tree<-drop.tip(tree, Phyllotaxis_drop$Tree)
labels<-data.frame(TipLabels(Phyllotaxis_tree)) 
labels$Tree<-labels$TipLabels.Phyllotaxis_tree.
d<-dplyr::inner_join(labels, Phyllotaxis)
data<-d %>% filter(d$wild.cultivated == "Cultivated")

#Prepare tree for model 
A <- ape::vcv.phylo(Phyllotaxis_tree)
B <-diag(diag(A + 0.00000000000000001)) #add to diagonal to make matrix positive definite.
row.names(B)<-row.names(A)

#Select only the cultivated crops, this is where the binary coding of changes lives in the data. Also seperate out and start with woody crops
woody<-data %>% filter(wild.cultivated == "Cultivated" & Herbaceous.woody == "Woody")
#set response for model
woody$size <- with(woody,WR2toCrop+WR2toWR1+Both+Neither) 
woody$y <- with(woody, cbind(Neither,Both,WR2toCrop,WR2toWR1)) 
#Run model

woody_Phyllotaxis <- brm(bf(y |trials(size) ~ 1 + (1|gr(Tree,cov=B))),
                         family = multinomial(),
                         data = data.frame(woody)
                         ,data2 = list(B=B),
                         chains = 4,
                         thin = 5,
                         iter = 3000,
                         warmup = 500,
                         sample_prior = TRUE,
                         cores = 4,
                         control = list(adapt_delta=0.95))
print(woody_Phyllotaxis)
saveRDS(woody_Phyllotaxis, file = "~/Dropbox/other_projects/crop_architecture/woody_Phyllotaxis.RDS")
# woody_Phyllotaxis <- readRDS("../woody_Phyllotaxis.RDS")

#Repeat steps above for herbaceous plans
herb<-data %>% filter(wild.cultivated == "Cultivated" & Herbaceous.woody == "Herbaceous")
herb$size <- with(herb,WR2toCrop+WR2toWR1+Both+Neither)
herb$y <- with(herb, cbind(Neither,Both,WR2toCrop,WR2toWR1))

herb_Phyllotaxis <- brm(bf(y |trials(size) ~ 1 + (1|gr(Tree,cov=B))),
                        family = multinomial(),
                        data = data.frame(herb),
                        data2 = list(B=B),
                        chains = 4,
                        thin = 5,
                        iter = 3000,
                        warmup = 500,
                        sample_prior = TRUE,
                        cores = 4,
                        control = list(adapt_delta = 0.95))
print(herb_Phyllotaxis)
saveRDS(herb_Phyllotaxis, file = "~/Dropbox/other_projects/crop_architecture/herb_Phyllotaxis.RDS")
# herb_Phyllotaxis <- readRDS("../herb_Phyllotaxis.RDS")

# Extract posterior distributions from the Models_tip,
ps_woody <- as_draws_df(woody_Phyllotaxis,
                        variable = c("b_muBoth_Intercept", "b_muWR2toCrop_Intercept", "b_muWR2toWR1_Intercept"))[,1:3]
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
  mutate(Trait = "Phyllotaxis") %>%
  melt()

ps_herb <- as_draws_df(herb_Phyllotaxis,
                       variable = c("b_muBoth_Intercept", "b_muWR2toCrop_Intercept", "b_muWR2toWR1_Intercept"))[,1:3]
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
  mutate(Trait = "Phyllotaxis") %>%
  melt()

##### Hypothesis testing 
dat1 <- as.data.frame(woody_Phyllotaxis)
h_woody_Phyllotaxis <- hypothesis(dat1, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

dat2 <- as.data.frame(herb_Phyllotaxis)
h_herb_Phyllotaxis <- hypothesis(dat2, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

# Plot
df_plot<-rbind(dfm_woody,dfm_herb) 
estimates_Phyllotaxis<-df_plot

phylo_plot<-ggplot(df_plot, aes(x=value*100,y=form,fill=variable)) +
  stat_slabinterval() +
  theme_bw() + 
  xlab("Probability of change (%)") + 
  ylab("Density") +
  ggsci::scale_fill_aaas(name="Direction") + 
  ggsci::scale_color_aaas(name="Direction")
phylo_plot

###### BranchingType Model #####
BranchingType<-trait %>% filter(Trait == "BranchingType" & VALUE != "NoData")
BranchingType$Tree <- BranchingType$Species.ssp.var
BranchingType_drop<-trait %>% filter(Trait == "BranchingType", VALUE == "NoData")
BranchingType_tree<-drop.tip(tree, BranchingType_drop$Tree)
labels<-data.frame(TipLabels(BranchingType_tree)) 
labels$Tree<-labels$TipLabels.BranchingType_tree.
d<-dplyr::inner_join(labels, BranchingType)
data<-d %>% filter(d$wild.cultivated == "Cultivated")

#Prepare tree for model 
A <- ape::vcv.phylo(BranchingType_tree)
B <-diag(diag(A + 0.00000000000000001)) #add to diagonal to make matrix positive definite.
row.names(B)<-row.names(A)

#Select only the cultivated crops, this is where the binary coding of changes lives in the data. Also seperate out and start with woody crops
woody<-data %>% filter(wild.cultivated == "Cultivated" & Herbaceous.woody == "Woody")
#set response for model
woody$size <- with(woody,WR2toCrop+WR2toWR1+Both+Neither) 
woody$y <- with(woody, cbind(Neither,Both,WR2toCrop,WR2toWR1)) 
#Run model

woody_BranchingType <- brm(bf(y |trials(size) ~ 1 + (1|gr(Tree,cov=B))),
                           family = multinomial(),
                           data = data.frame(woody)
                           ,data2 = list(B=B),
                           chains = 4,
                           thin = 5,
                           iter = 3000,
                           warmup = 500,
                           sample_prior = TRUE,
                           cores = 4,
                           control = list(adapt_delta = 0.85))
print(woody_BranchingType)
saveRDS(woody_BranchingType, file = "~/Dropbox/other_projects/crop_architecture/woody_BranchingType.RDS")
# woody_BranchingType <- readRDS("../woody_BranchingType.RDS")

#Repeat steps above for herbaceous plans
herb<-data %>% filter(wild.cultivated == "Cultivated" & Herbaceous.woody == "Herbaceous")
herb$size <- with(herb,WR2toCrop+WR2toWR1+Both+Neither)
herb$y <- with(herb, cbind(Neither,Both,WR2toCrop,WR2toWR1))

herb_BranchingType <- brm(bf(y |trials(size) ~ 1 + (1|gr(Tree,cov=B))),
                          family = multinomial(),
                          data = data.frame(herb),
                          data2 = list(B=B),
                          chains = 4,
                          thin = 5,
                          iter = 3000,
                          warmup = 500,
                          sample_prior = TRUE,
                          cores = 4,
                          control = list(adapt_delta = 0.9))
print(herb_BranchingType)
saveRDS(herb_BranchingType, file = "~/Dropbox/other_projects/crop_architecture/herb_BranchingType.RDS")
# herb_BranchingType <- readRDS("../herb_BranchingType.RDS")

# Extract posterior distributions from the Models_tip,
ps_woody <- as_draws_df(woody_BranchingType,
                        variable = c("b_muBoth_Intercept", "b_muWR2toCrop_Intercept", "b_muWR2toWR1_Intercept"))[,1:3]
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
  mutate(Trait = "Branching Type") %>%
  melt()

ps_herb <- as_draws_df(herb_BranchingType,
                       variable = c("b_muBoth_Intercept", "b_muWR2toCrop_Intercept", "b_muWR2toWR1_Intercept"))[,1:3]
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
  mutate(Trait = "Branching Type") %>%
  melt()

##### Hypothesis testing 
dat1 <- as.data.frame(woody_BranchingType)
h_woody_BranchingType <- hypothesis(dat1, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

dat2 <- as.data.frame(herb_BranchingType)
h_herb_BranchingType <- hypothesis(dat2, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

# Plot
df_plot<-rbind(dfm_woody,dfm_herb) 
estimates_BranchingType<-df_plot

phylo_plot<-ggplot(df_plot, aes(x=value*100,y=form,fill=variable)) +
  stat_slabinterval() +
  theme_bw() + 
  xlab("Probability of change (%)") + 
  ylab("Density") +
  ggsci::scale_fill_aaas(name="Direction") + 
  ggsci::scale_color_aaas(name="Direction")
phylo_plot


###### BranchingPosition Model #####
BranchingPosition<-trait %>% filter(Trait == "BranchingPosition" & VALUE != "NoData")
BranchingPosition$Tree <- BranchingPosition$Species.ssp.var
BranchingPosition_drop<-trait %>% filter(Trait == "BranchingPosition", VALUE == "NoData")
BranchingPosition_tree<-drop.tip(tree, BranchingPosition_drop$Tree)
labels<-data.frame(TipLabels(BranchingPosition_tree)) 
labels$Tree<-labels$TipLabels.BranchingPosition_tree.
d<-dplyr::inner_join(labels, BranchingPosition)
data<-d %>% filter(d$wild.cultivated == "Cultivated")

#Prepare tree for model 
A <- ape::vcv.phylo(BranchingPosition_tree)
B <-diag(diag(A + 0.00000000000000001)) #add to diagonal to make matrix positive definite.
row.names(B)<-row.names(A)

#Select only the cultivated crops, this is where the binary coding of changes lives in the data. Also seperate out and start with woody crops
woody<-data %>% filter(wild.cultivated == "Cultivated" & Herbaceous.woody == "Woody")
#set response for model
woody$size <- with(woody,WR2toCrop+WR2toWR1+Both+Neither) 
woody$y <- with(woody, cbind(Neither,Both,WR2toCrop,WR2toWR1)) 
#Run model

woody_BranchingPosition <- brm(bf(y |trials(size) ~ 1 + (1|gr(Tree,cov=B))),
                               family = multinomial(),
                               data = data.frame(woody)
                               ,data2 = list(B=B),
                               chains = 4,
                               thin = 5,
                               iter = 3000,
                               warmup = 500,
                               sample_prior = TRUE,
                               cores = 4)

print(woody_BranchingPosition)
saveRDS(woody_BranchingPosition, file = "~/Dropbox/other_projects/crop_architecture/woody_BranchingPosition.RDS")
# woody_BranchingPosition <- readRDS("../woody_BranchingPosition.RDS")

#Repeat steps above for herbaceous plans
herb<-data %>% filter(wild.cultivated == "Cultivated" & Herbaceous.woody == "Herbaceous")
herb$size <- with(herb,WR2toCrop+WR2toWR1+Both+Neither)
herb$y <- with(herb, cbind(Neither,Both,WR2toCrop,WR2toWR1))

herb_BranchingPosition <- brm(bf(y |trials(size) ~ 1 + (1|gr(Tree,cov=B))),
                              family = multinomial(),
                              data = data.frame(herb),
                              data2 = list(B=B),
                              chains = 4,
                              thin = 5,
                              iter = 3000,
                              warmup = 500,
                              sample_prior = TRUE,
                              cores = 4)

print(herb_BranchingPosition)
saveRDS(herb_BranchingPosition, file = "~/Dropbox/other_projects/crop_architecture/herb_BranchingPosition.RDS")
# herb_BranchingPosition <- readRDS("../herb_BranchingPosition.RDS")

# Extract posterior distributions from the Models_tip,
ps_woody <- as_draws_df(woody_BranchingPosition,
                        variable = c("b_muBoth_Intercept", "b_muWR2toCrop_Intercept", "b_muWR2toWR1_Intercept"))[,1:3]
names(ps_herb)
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
  mutate(Trait = "Branching Position") %>%
  melt()

ps_herb <- as_draws_df(herb_BranchingPosition,
                       variable = c("b_muBoth_Intercept", "b_muWR2toCrop_Intercept", "b_muWR2toWR1_Intercept"))[,1:3]
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
  mutate(Trait = "Branching Position") %>%
  melt()

x<-dfm_herb  %>% filter(variable == "Crop")
mean(x$value)
##### Hypothesis testing 
dat1 <- as.data.frame(woody_BranchingPosition)
h_woody_BranchingPosition <- hypothesis(dat1, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

dat2 <- as.data.frame(herb_BranchingPosition)
h_herb_BranchingPosition <- hypothesis(dat2, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

# Plot
df_plot<-rbind(dfm_woody,dfm_herb) 
estimates_BranchingPosition<-df_plot

phylo_plot<-ggplot(df_plot, aes(x=value*100,y=form,fill=variable)) +
  stat_slabinterval() +
  theme_bw() + 
  xlab("Probability of change (%)") + 
  ylab("Density") +
  ggsci::scale_fill_aaas(name="Direction") + 
  ggsci::scale_color_aaas(name="Direction")
phylo_plot


###### FloweringAxis Model #####
Flowering<-trait %>% filter(Trait == "Flowering" & VALUE != "NoData")
Flowering$Tree <- Flowering$Species.ssp.var 
Flowering_drop<-trait %>% filter(Trait == "Flowering", VALUE == "NoData")
Flowering_tree<-drop.tip(tree, Flowering_drop$Tree)
labels<-data.frame(TipLabels(Flowering_tree)) 
labels$Tree<-labels$TipLabels.Flowering_tree.
d<-dplyr::inner_join(labels, Flowering)
data<-d %>% filter(d$wild.cultivated == "Cultivated")

#Prepare tree for model 
A <- ape::vcv.phylo(Flowering_tree)
B <-diag(diag(A + 0.00000000000000001)) #add to diagonal to make matrix positive definite.
row.names(B)<-row.names(A)

#Select only the cultivated crops, this is where the binary coding of changes lives in the data. Also seperate out and start with woody crops
woody<-data %>% filter(wild.cultivated == "Cultivated" & Herbaceous.woody == "Woody")
#set response for model
woody$size <- with(woody,WR2toCrop+WR2toWR1+Both+Neither) 
woody$y <- with(woody, cbind(Neither,Both,WR2toCrop,WR2toWR1)) 
#Run model

woody_Flowering <- brm(bf(y |trials(size) ~ 1 + (1|gr(Tree,cov=B))),
                           family = multinomial(),
                           data = data.frame(woody)
                           ,data2 = list(B=B),
                           chains = 4,
                           thin = 5,
                           iter = 3000,
                           warmup = 500,
                           sample_prior = TRUE,
                           cores = 4,
                       control = list(adapt_delta=0.85))
print(woody_Flowering)
saveRDS(woody_Flowering, file = "~/Dropbox/other_projects/crop_architecture/woody_Flowering.RDS")
# woody_Flowering <- readRDS("../woody_Flowering.RDS")

#Repeat steps above for herbaceous plans
herb<-data %>% filter(wild.cultivated == "Cultivated" & Herbaceous.woody == "Herbaceous")
herb$size <- with(herb,WR2toCrop+WR2toWR1+Both+Neither)
herb$y <- with(herb, cbind(Neither,Both,WR2toCrop,WR2toWR1))

herb_Flowering <- brm(bf(y |trials(size) ~ 1 + (1|gr(Tree,cov=B))),
                      family = multinomial(),
                      data = data.frame(herb),
                      data2 = list(B=B),
                      chains = 4,
                      thin = 5,
                      iter = 3000,
                      warmup = 500,
                      sample_prior = TRUE,
                      cores = 4,
                      control = list(adapt_delta=0.95)
                      )
print(herb_Flowering)
saveRDS(herb_Flowering, file = "~/Dropbox/other_projects/crop_architecture/herb_Flowering.RDS")
# herb_Flowering <- readRDS("../herb_Flowering.RDS")

# Extract posterior distributions from the Models_tip,
ps_woody <- as_draws_df(woody_Flowering,
                        variable = c("b_muBoth_Intercept", "b_muWR2toCrop_Intercept", "b_muWR2toWR1_Intercept"))[,1:3]
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
  mutate(Trait = "Flowering Axis") %>%
  melt()

ps_herb <- as_draws_df(herb_Flowering,
                       variable = c("b_muBoth_Intercept", "b_muWR2toCrop_Intercept", "b_muWR2toWR1_Intercept"))[,1:3]
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
  mutate(Trait = "Flowering Axis") %>%
  melt()

##### Hypothesis testing 
dat1 <- as.data.frame(woody_Flowering)
h_woody_Flowering <- hypothesis(dat1, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

dat2 <- as.data.frame(herb_Flowering)
h_herb_Flowering <- hypothesis(dat2, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

# Plot
df_plot<-rbind(dfm_woody,dfm_herb) 
estimates_Flowering<-df_plot

phylo_plot<-ggplot(df_plot, aes(x=value*100,y=form,fill=variable)) +
  stat_slabinterval() +
  theme_bw() + 
  xlab("Probability of change (%)") + 
  ylab("Density") +
  ggsci::scale_fill_aaas(name="Direction") + 
  ggsci::scale_color_aaas(name="Direction")
phylo_plot

###### MeristemFunction Model #####
MeristemFunction<-trait %>% filter(Trait == "MeristemFunctioning" & VALUE != "NoData")
MeristemFunction$Tree <- MeristemFunction$Species.ssp.var
MeristemFunction_drop<-trait %>% filter(Trait == "MeristemFunctioning", VALUE == "NoData")
MeristemFunction_tree<-drop.tip(tree, MeristemFunction_drop$Tree)
labels<-data.frame(TipLabels(MeristemFunction_tree)) 
labels$Tree<-labels$TipLabels.MeristemFunction_tree.
d<-dplyr::inner_join(labels, MeristemFunction)
data<-d %>% filter(d$wild.cultivated == "Cultivated")

#Prepare tree for model 
A <- ape::vcv.phylo(MeristemFunction_tree)
B <-diag(diag(A + 0.00000000000000001)) #add to diagonal to make matrix positive definite.
row.names(B)<-row.names(A)

#Select only the cultivated crops, this is where the binary coding of changes lives in the data. Also seperate out and start with woody crops
woody<-data %>% filter(wild.cultivated == "Cultivated" & Herbaceous.woody == "Woody")
#set response for model
woody$size <- with(woody,WR2toCrop+WR2toWR1+Both+Neither) 
woody$y <- with(woody, cbind(Neither,Both,WR2toCrop,WR2toWR1)) 
#Run model

woody_MeristemFunction <- brm(bf(y |trials(size) ~ 1 + (1|gr(Tree,cov=B))),
                              family = multinomial(),
                              data = data.frame(woody)
                              ,data2 = list(B=B),
                              chains = 4,
                              thin = 5,
                              iter = 3000,
                              warmup = 500,
                              sample_prior = TRUE,
                              cores = 4,
                              control = list(adapt_delta=0.85))
print(woody_MeristemFunction)
saveRDS(woody_MeristemFunction, file = "~/Dropbox/other_projects/crop_architecture/woody_MeristemFunction.RDS")
# woody_MeristemFunction <- readRDS("../woody_MeristemFunction.RDS")

#Repeat steps above for herbaceous plans
herb<-data %>% filter(wild.cultivated == "Cultivated" & Herbaceous.woody == "Herbaceous")
herb$size <- with(herb,WR2toCrop+WR2toWR1+Both+Neither)
herb$y <- with(herb, cbind(Neither,Both,WR2toCrop,WR2toWR1))

herb_MeristemFunction <- brm(bf(y |trials(size) ~ 1 + (1|gr(Tree,cov=B))),
                             family = multinomial(),
                             data = data.frame(herb),
                             data2 = list(B=B),
                             chains = 4,
                             thin = 5,
                             iter = 3000,
                             warmup = 500,
                             sample_prior = TRUE,
                             cores = 4,
                             control = list(adapt_delta = 0.9))
print(herb_MeristemFunction)
saveRDS(herb_MeristemFunction, file = "~/Dropbox/other_projects/crop_architecture/herb_MeristemFunction.RDS")
# herb_MeristemFunction <- readRDS("../herb_MeristemFunction.RDS")

# Extract posterior distributions from the Models_tip,
ps_woody <- as_draws_df(woody_MeristemFunction,
                        variable = c("b_muBoth_Intercept", "b_muWR2toCrop_Intercept", "b_muWR2toWR1_Intercept"))[,1:3]
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
  mutate(Trait = "Meristem Function") %>%
  melt()

ps_herb <- as_draws_df(herb_MeristemFunction,
                       variable = c("b_muBoth_Intercept", "b_muWR2toCrop_Intercept", "b_muWR2toWR1_Intercept"))[,1:3]
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
  mutate(Trait = "Meristem Function") %>%
  melt()

##### Hypothesis testing 
dat1 <- as.data.frame(woody_MeristemFunction)
h_woody_MeristemFunction <- hypothesis(dat1, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

dat2 <- as.data.frame(herb_MeristemFunction)
h_herb_MeristemFunction <- hypothesis(dat2, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

# Plot
df_plot<-rbind(dfm_woody,dfm_herb) 
estimates_MeristemFunction<-df_plot

phylo_plot<-ggplot(df_plot, aes(x=value*100,y=form,fill=variable)) +
  stat_slabinterval() +
  theme_bw() + 
  xlab("Probability of change (%)") + 
  ylab("Density") +
  ggsci::scale_fill_aaas(name="Direction") + 
  ggsci::scale_color_aaas(name="Direction")
phylo_plot


###### BranchingMechanism Model #####
BranchingMechanism<-trait %>% filter(Trait == "BranchingMechanism" & VALUE != "NoData")
BranchingMechanism$Tree <- BranchingMechanism$Species.ssp.var
BranchingMechanism_drop<-trait %>% filter(Trait == "BranchingMechanism", VALUE == "NoData")
BranchingMechanism_tree<-drop.tip(tree, BranchingMechanism_drop$Tree)
labels<-data.frame(TipLabels(BranchingMechanism_tree)) 
labels$Tree<-labels$TipLabels.BranchingMechanism_tree.
d<-dplyr::inner_join(labels, BranchingMechanism)
data<-d %>% filter(d$wild.cultivated == "Cultivated")

#Prepare tree for model 
A <- ape::vcv.phylo(BranchingMechanism_tree)
B <-diag(diag(A + 0.00000000000000001)) #add to diagonal to make matrix positive definite.
row.names(B)<-row.names(A)

#Select only the cultivated crops, this is where the binary coding of changes lives in the data. Also seperate out and start with woody crops
woody<-data %>% filter(wild.cultivated == "Cultivated" & Herbaceous.woody == "Woody")
#set response for model
woody$size <- with(woody,WR2toCrop+WR2toWR1+Both+Neither) 
woody$y <- with(woody, cbind(Neither,Both,WR2toCrop,WR2toWR1)) 
#Run model
woody_BranchingMechanism <- brm(bf(y |trials(size) ~ 1 + (1|gr(Tree,cov=B))),
                                family = multinomial(),
                                data = data.frame(woody)
                                ,data2 = list(B=B),
                                chains = 4,
                                thin = 5,
                                iter = 3000,
                                warmup = 500,
                                sample_prior = TRUE,
                                cores = 4,
                                control = list(adapt_delta = 0.85))
print(woody_BranchingMechanism)
saveRDS(woody_BranchingMechanism, file = "~/Dropbox/other_projects/crop_architecture/woody_BranchingMechanism.RDS")
# woody_BranchingMechanism <- readRDS("../woody_BranchingMechanism.RDS")

#Repeat steps above for herbaceous plans
herb<-data %>% filter(wild.cultivated == "Cultivated" & Herbaceous.woody == "Herbaceous")
herb$size <- with(herb,WR2toCrop+WR2toWR1+Both+Neither)
herb$y <- with(herb, cbind(Neither,Both,WR2toCrop,WR2toWR1))

herb_BranchingMechanism <- brm(bf(y |trials(size) ~ 1 + (1|gr(Tree,cov=B))),
                               family = multinomial(),
                               data = data.frame(herb),
                               data2 = list(B=B),
                               chains = 4,
                               thin = 5,
                               iter = 3000,
                               warmup = 500,
                               sample_prior = TRUE,
                               cores = 4,
                               control = list(adapt_delta = 0.95)
                               )
print(herb_BranchingMechanism)
saveRDS(herb_BranchingMechanism, file = "~/Dropbox/other_projects/crop_architecture/herb_BranchingMechanism.RDS")
# herb_BranchingMechanism <- readRDS("../herb_BranchingMechanism.RDS")

# Extract posterior distributions from the Models_tip,
ps_woody <- as_draws_df(woody_BranchingMechanism,
                        variable = c("b_muBoth_Intercept", "b_muWR2toCrop_Intercept", "b_muWR2toWR1_Intercept"))[,1:3]
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

ps_herb <- as_draws_df(herb_BranchingMechanism,
                       variable = c("b_muBoth_Intercept", "b_muWR2toCrop_Intercept", "b_muWR2toWR1_Intercept"))[,1:3]
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
dat1 <- as.data.frame(woody_BranchingMechanism)
h_woody_BranchingMechanism <- hypothesis(dat1, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

dat2 <- as.data.frame(herb_BranchingMechanism)
h_herb_BranchingMechanism <- hypothesis(dat2, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

# Plot
df_plot<-rbind(dfm_woody,dfm_herb) 
estimates_BranchingMechanism<-df_plot

phylo_plot<-ggplot(df_plot, aes(x=value*100,y=form,fill=variable)) +
  stat_slabinterval() +
  theme_bw() + 
  xlab("Probability of change (%)") + 
  ylab("Density") +
  ggsci::scale_fill_aaas(name="Direction") + 
  ggsci::scale_color_aaas(name="Direction")
phylo_plot


####### MAYBE REMOVE  ###### ShortShoots Model #####
ShortShoots<-trait %>% filter(Trait == "ShortShoots" & VALUE != "NoData")
ShortShoots$Tree <- ShortShoots$Species.ssp.var
ShortShoots_drop<-trait %>% filter(Trait == "ShortShoots", VALUE == "NoData")
ShortShoots_tree<-drop.tip(tree, ShortShoots_drop$Tree)
labels<-data.frame(TipLabels(ShortShoots_tree))
labels$Tree<-labels$TipLabels.ShortShoots_tree.
d<-dplyr::inner_join(labels, ShortShoots)
data<-d %>% filter(d$wild.cultivated == "Cultivated")

#Prepare tree for model
A <- ape::vcv.phylo(ShortShoots_tree)
B <-diag(diag(A + 0.00000000000000001)) #add to diagonal to make matrix positive definite.
row.names(B)<-row.names(A)

#Select only the cultivated crops, this is where the binary coding of changes lives in the data. Also seperate out and start with woody crops
woody<-data %>% filter(wild.cultivated == "Cultivated" & Herbaceous.woody == "Woody")
#set response for model
woody$size <- with(woody,WR2toCrop+WR2toWR1+Both+Neither)
woody$y <- with(woody, cbind(Neither,Both,WR2toCrop,WR2toWR1))
#Run model

woody_ShortShoots <- brm(bf(y |trials(size) ~ 1 + (1|gr(Tree,cov=B))),
                         family = multinomial(),
                         data = data.frame(woody)
                         ,data2 = list(B=B),
                         chains = 4,
                         thin = 5,
                         iter = 3000,
                         warmup = 500,
                         sample_prior = TRUE,
                         cores = 4,
                         control = list(adapt_delta = 0.9))
print(woody_ShortShoots)
saveRDS(woody_ShortShoots, file = "~/Dropbox/other_projects/crop_architecture/woody_ShortShoots.RDS")
# woody_ShortShoots <- readRDS("../woody_ShortShoots.RDS")

#Repeat steps above for herbaceous plans
herb<-data %>% filter(wild.cultivated == "Cultivated" & Herbaceous.woody == "Herbaceous")
herb$size <- with(herb,WR2toCrop+WR2toWR1+Both+Neither)
herb$y <- with(herb, cbind(Neither,Both,WR2toCrop,WR2toWR1))

herb_ShortShoots <- brm(bf(y |trials(size) ~ 1 + (1|gr(Tree,cov=B))),
                        family = multinomial(),
                        data = data.frame(herb),
                        data2 = list(B=B),
                        chains = 4,
                        thin = 5,
                        iter = 3000,
                        warmup = 500,
                        sample_prior = TRUE,
                        cores = 4,
                        control = list(adapt_delta = 0.9))
print(herb_ShortShoots)
saveRDS(herb_ShortShoots, file = "~/Dropbox/other_projects/crop_architecture/herb_ShortShoots.RDS")
# herb_ShortShoots <- readRDS("../herb_ShortShoots.RDS")

# Extract posterior distributions from the Models_tip,
ps_woody <- as_draws_df(woody_ShortShoots,
                        variable = c("b_muBoth_Intercept", "b_muWR2toCrop_Intercept", "b_muWR2toWR1_Intercept"))[,1:3]
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
  mutate(Trait = "Short Shoots") %>%
  melt()

ps_herb <- as_draws_df(herb_ShortShoots,
                       variable = c("b_muBoth_Intercept", "b_muWR2toCrop_Intercept", "b_muWR2toWR1_Intercept"))[,1:3]
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
  mutate(Trait = "Short Shoots") %>%
  melt()

##### Hypothesis testing
dat1 <- as.data.frame(woody_ShortShoots)
h_woody_ShortShoots <- hypothesis(dat1, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

dat2 <- as.data.frame(herb_ShortShoots)
h_herb_ShortShoots <- hypothesis(dat2, "b_muWR2toCrop_Intercept > b_muWR2toWR1_Intercept")

# Plot
df_plot<-rbind(dfm_woody,dfm_herb)
estimates_ShortShoots<-df_plot

phylo_plot<-ggplot(df_plot, aes(x=value*100,y=form,fill=variable)) +
  stat_slabinterval() +
  theme_bw() +
  xlab("Probability of change (%)") +
  ylab("Density") +
  ggsci::scale_fill_aaas(name="Direction") +
  ggsci::scale_color_aaas(name="Direction")
phylo_plot



###### Plotting Model ######
estimates<-rbind(estimates_GrowthDirection,estimates_Phyllotaxis,
                 estimates_BranchingType,estimates_BranchingPosition,
                 estimates_Flowering,estimates_MeristemFunction,estimates_BranchingMechanism,
                 estimates_ShortShoots) %>%
  mutate(form = ifelse(form == "herb", "Herbaceous", form)) %>%
  mutate(form = ifelse(form == "woody", "Woody", form))

d<-data.frame(estimates_Phyllotaxis  %>% filter(form=="herb"))
sd(d$value)

p.values <- c(h_woody_GrowthDirection$hypothesis$Post.Prob, h_herb_GrowthDirection$hypothesis$Post.Prob,
              h_woody_Phyllotaxis$hypothesis$Post.Prob, h_herb_Phyllotaxis$hypothesis$Post.Prob,
              h_woody_BranchingType$hypothesis$Post.Prob, h_herb_BranchingType$hypothesis$Post.Prob,
              h_woody_BranchingPosition$hypothesis$Post.Prob, h_herb_BranchingPosition$hypothesis$Post.Prob,
              h_woody_Flowering$hypothesis$Post.Prob, h_herb_Flowering$hypothesis$Post.Prob,
              h_woody_MeristemFunction$hypothesis$Post.Prob, h_herb_MeristemFunction$hypothesis$Post.Prob,
              h_woody_BranchingMechanism$hypothesis$Post.Prob, h_herb_BranchingMechanism$hypothesis$Post.Prob,
              h_woody_ShortShoots$hypothesis$Post.Prob, h_herb_ShortShoots$hypothesis$Post.Prob)

p.values <- data.frame(Trait = c(rep("Growth Direction", 2),
                                 rep("Phyllotaxis", 2),
                                 rep("Branching Type", 2),
                                 rep("Branching Position", 2),
                                 rep("Flowering Axis", 2),
                                 rep("Meristem Function", 2),
                                 rep("Branching Mechanism", 2),
                                 rep("Short Shoots", 2)),
                       form = rep(c("Woody", "Herbaceous"), 8),
                       p = p.values
                       )

groups <- unique(estimates$Trait) # get a list of unique groups
plot_list <- list() # initialize empty list to store plots

for (g in groups) {
  # create a ggplot for the current group
  plot_data <- estimates[estimates$Trait == g, ]
  ps <- p.values[p.values$Trait == g, ]$p
  ps <- signif(ps*100, 2)
  if (g %in% c("Branching Position", "Phyllotaxis")) {
    p <- ggplot(plot_data, aes(x=value*100,y=form, fill=variable)) +
      stat_slabinterval() +
      ggsci::scale_fill_aaas(name="Direction") + 
      ggsci::scale_color_aaas(name="Direction") +
      labs(title = g,
           subtitle = paste0("P woody = ", ps[1], "% | ", "P herb = ", ps[2], "%"),
             x = "Probability of change (%)") +
      theme_bw() +
      theme(plot.title = element_text(size=12),
            axis.title.y = element_blank(),
            plot.subtitle = element_text(size=10))
    
    # save the ggplot in the list
    plot_list[[g]] <- p
  } else {
    p <- ggplot(plot_data, aes(x=value*100,y=form, fill=variable)) +
      stat_slabinterval() +
      ggsci::scale_fill_aaas(name="Direction") + 
      ggsci::scale_color_aaas(name="Direction") +
      labs(title = g,
           subtitle = paste0("P woody = ", ps[1], "% | ", "P herb = ", ps[2], "%"),
             x = "Probability of change (%)") +
      theme_bw() +
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size=12),
            plot.subtitle = element_text(size=10))
    
    # save the ggplot in the list
    plot_list[[g]] <- p
  }
  
}

# print the ggplots from the list
for (i in seq_along(plot_list)) {
  print(plot_list[[i]])
}

ordered_list<-list(plot_list$`Branching Position`,plot_list$`Growth Direction`,
                   plot_list$`Branching Mechanism`,plot_list$`Flowering Axis`,
                   plot_list$Phyllotaxis,plot_list$`Meristem Function`,
                   plot_list$`Branching Type`,plot_list$`Short Shoots`
                   )

ggarrange(plotlist = ordered_list,ncol = 4,nrow = 2,common.legend = TRUE,
          legend="right", widths = c(1.25,1,1,1))



###### Plotting Hypothesis Testing #####
results<-read.csv(file = "Data/hypothesistesting.csv", header=T) %>%
  rename("Posterior [Crop > Wild]" = Estimate) %>%
  rename("Posterior probabability" = Posterior.probability) %>%
  select(-"Posterior [Crop > Wild]") %>%
  melt()

ggplot(results, aes(x=Form, y=Trait, fill=value)) + 
  geom_tile() + 
  facet_grid(~variable) + 
  scale_fill_viridis_c(begin = 0,end = 1,option = "mako")  +
  theme_bw() 

summary(results$value)
