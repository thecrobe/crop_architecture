library(ggplot2) 
library(dplyr)
library(stringr)
library(flipTables)
library(ggsci)
library(ggpubr)
library(tidyr)

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

# get the relevant end states - note that this will only pull rows with transitions
# not every row
trait$WR2State <- str_split_fixed(trait$Direction.WR2toCrop, " to ", 2)[,1]
trait$CropState <- str_split_fixed(trait$Direction.WR2toCrop, " to ", 2)[,2]
trait$WR1State <- str_split_fixed(trait$Direction.WR2toWR1, " to ", 2)[,2]

# prepare dataframe for plotting



###### GrowthDirection #####

GrowthDirection <- trait %>% filter(Trait == "GrowthDirection", VALUE != "NoData")
data <- GrowthDirection %>% filter(GrowthDirection$wild.cultivated == "Cultivated")
data <- data %>% filter(data$Neither == 0)

table(data$CropState)
table(data$WR1State)
# table(data$WR2State)

# now we need to extract the relevant transitions
data <- data %>% filter(data$Both == 0)

# collapse some states

# Axis 1 Orthotropic -> Orthotropic
data[which(data$CropState == "Axis 1 Orthotropic"),]$CropState <- "Orthotropic"
data[which(data$WR1State == "Axis 1 Orthotropic"),]$WR1State <- "Orthotropic"
data[which(data$WR2State == "Axis 1 Orthotropic"),]$WR2State <- "Orthotropic"

# Drop Orthotropic and Plagiotropic
data <- data[-which(data$CropState == "Orthotropic and Plagiotropic"),]

tab1 <- table(data$CropState)
tab1 <- tab1[2:length(tab1)]
tab2 <- table(data$WR1State)
tab2 <- tab2[2:length(tab2)]

merged <- t(Merge2Tables(tab1, tab2, direction = "S"))
rownames(merged) <- c("Crop", "Wild")
merged[is.na(merged)] <- 0
chisq.test(merged, simulate.p.value = T)

# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# 
# data:  merged
# X-squared = 5.781, df = NA, p-value = 0.1189

# 
# # do posthoc OR CI
# 
# # ageotropic
# 
# GrowthDirection_OR <- data.frame(trait = colnames(merged),
#                                  or = rep(0, ncol(merged)),
#                                  ci.l = rep(0, ncol(merged)),
#                                  ci.u = rep(0, ncol(merged)))
# 
# for (i in colnames(merged)) {
#   tab <- cbind(merged[,match(i, colnames(merged))],
#                rowSums(merged[,colnames(merged)[-match(i, colnames(merged))]]))
#   f <- fisher.test(tab)
#   or <- f$estimate
#   ci.l <- f$conf.int[1]
#   ci.u <- f$conf.int[2]
#   GrowthDirection_OR[which(GrowthDirection_OR$trait == i),2:4] <- c(or, ci.l, ci.u)
# }
# 
# tab <- cbind(merged[,1], rowSums(merged[,c(2, 3, 4)]))

tmp <- data.frame(trait = rep("GrowthDirection", 2*ncol(merged)),
                  transition = c(rep("Crop", ncol(merged)), rep("Wild", ncol(merged))),
                  state = rep(colnames(merged),2),
                  count = c(merged[1,], merged[2,]),
                  prop = c(merged[1,]/sum(merged[1,]), merged[2,]/sum(merged[2,]))
                  )

plotdf <- tmp

###### Phyllotaxis #####

Phyllotaxis <- trait %>% filter(Trait == "Phyllotaxis", VALUE != "NoData")
data <- Phyllotaxis %>% filter(Phyllotaxis$wild.cultivated == "Cultivated")
data <- data %>% filter(data$Neither == 0)

table(data$CropState)
table(data$WR1State)
# table(data$WR2State)

# now we need to extract the relevant transitions
data <- data %>% filter(data$Both == 0)

tab1 <- table(data$CropState)
tab1 <- tab1[2:length(tab1)]
tab2 <- table(data$WR1State)
tab2 <- tab2[2:length(tab2)]

merged <- t(Merge2Tables(tab1, tab2, direction = "S"))
rownames(merged) <- c("Crop", "Wild")
merged[is.na(merged)] <- 0
chisq.test(merged, simulate.p.value = T)

# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# 
# data:  merged
# X-squared = 6.619, df = NA, p-value = 0.2884

tmp <- data.frame(trait = rep("Phyllotaxis", 2*ncol(merged)),
                  transition = c(rep("Crop", ncol(merged)), rep("Wild", ncol(merged))),
                  state = rep(colnames(merged),2),
                  count = c(merged[1,], merged[2,]),
                  prop = c(merged[1,]/sum(merged[1,]), merged[2,]/sum(merged[2,]))
)

plotdf <- rbind(plotdf, tmp)

###### BranchingType #####

BranchingType <- trait %>% filter(Trait == "BranchingType", VALUE != "NoData")
data <- BranchingType %>% filter(BranchingType$wild.cultivated == "Cultivated")
data <- data %>% filter(data$Neither == 0)

table(data$CropState)
table(data$WR1State)
# table(data$WR2State)

# now we need to extract the relevant transitions
data <- data %>% filter(data$Both == 0)

# collapse some states

# Axis 1 Monopodial -> Monopodial
data[which(data$CropState == "Axis 1 Monopodial"),]$CropState <- "Monopodial"
data[which(data$WR1State == "Axis 1 Monopodial"),]$WR1State <- "Monopodial"
# data[which(data$WR2State == "Axis 1 Monopodial"),]$WR2State <- "Monopodial"

# drop Monopodial or Sympodial
data <- data[-which(data$CropState == "Monopodial or Sympodial"),]

tab1 <- table(data$CropState)
tab1 <- tab1[2:length(tab1)]
tab2 <- table(data$WR1State)
tab2 <- tab2[2:length(tab2)]

merged <- t(Merge2Tables(tab1, tab2, direction = "S"))
rownames(merged) <- c("Crop", "Wild")
merged[is.na(merged)] <- 0
chisq.test(merged, simulate.p.value = T)

# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# 
# data:  merged
# X-squared = 11.927, df = NA, p-value = 0.01299

round(residuals(chisq.test(merged, simulate.p.value = T)), 3)
#       Absent Axis 1 Monopodial Axis 2 Sympodial Monopodial and Sympodial Monopodial Sympodial
# Crop  1.354                             -0.012                   -1.160     -0.719    -0.461
# Wild -1.940                              0.017                    1.661      1.030     0.660

tmp <- data.frame(trait = rep("BranchingType", 2*ncol(merged)),
                  transition = c(rep("Crop", ncol(merged)), rep("Wild", ncol(merged))),
                  state = rep(colnames(merged),2),
                  count = c(merged[1,], merged[2,]),
                  prop = c(merged[1,]/sum(merged[1,]), merged[2,]/sum(merged[2,]))
)

plotdf <- rbind(plotdf, tmp)

###### BranchingPosition #####

BranchingPosition <- trait %>% filter(Trait == "BranchingPosition", VALUE != "NoData")
data <- BranchingPosition %>% filter(BranchingPosition$wild.cultivated == "Cultivated")
data <- data %>% filter(data$Neither == 0)

table(data$CropState)
table(data$WR1State)
# table(data$WR2State)

# now we need to extract the relevant transitions
data <- data %>% filter(data$Both == 0)

tab1 <- table(data$CropState)
tab1 <- tab1[2:length(tab1)]
tab2 <- table(data$WR1State)
tab2 <- tab2[2:length(tab2)]

merged <- t(Merge2Tables(tab1, tab2, direction = "S"))
rownames(merged) <- c("Crop", "Wild")
merged[is.na(merged)] <- 0
chisq.test(merged, simulate.p.value = T)
# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# 
# data:  merged
# X-squared = 12.13, df = NA, p-value = 0.02699
round(residuals(chisq.test(merged, simulate.p.value = T)), 3)
#       Absent Acrotonic Acrotonic and Basitonic Basitonic Basitonic and Mesotonic Mesotonic
# Crop  1.129     1.390                  -0.348    -0.725                   0.390    -0.954
# Wild -1.368    -1.684                   0.421     0.878                  -0.473     1.156

tmp <- data.frame(trait = rep("BranchingPosition", 2*ncol(merged)),
                  transition = c(rep("Crop", ncol(merged)), rep("Wild", ncol(merged))),
                  state = rep(colnames(merged),2),
                  count = c(merged[1,], merged[2,]),
                  prop = c(merged[1,]/sum(merged[1,]), merged[2,]/sum(merged[2,]))
)

plotdf <- rbind(plotdf, tmp)

###### FloweringAxis #####

Flowering <- trait %>% filter(Trait == "Flowering", VALUE != "NoData")
data <- Flowering %>% filter(Flowering$wild.cultivated == "Cultivated")
data <- data %>% filter(data$Neither == 0)

table(data$CropState)
table(data$WR1State)
# table(data$WR2State)

# now we need to extract the relevant transitions
data <- data %>% filter(data$Both == 0)

# drop Absent only present in a single crop
data <- data[-which(data$CropState == "Absent"),]


tab1 <- table(data$CropState)
tab1 <- tab1[2:length(tab1)]
tab2 <- table(data$WR1State)
tab2 <- tab2[2:length(tab2)]

merged <- t(Merge2Tables(tab1, tab2, direction = "S"))
rownames(merged) <- c("Crop", "Wild")
merged[is.na(merged)] <- 0
chisq.test(merged, simulate.p.value = T)

# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# 
# data:  merged
# X-squared = 3.5638, df = NA, p-value = 0.2964

tmp <- data.frame(trait = rep("FloweringAxis", 2*ncol(merged)),
                  transition = c(rep("Crop", ncol(merged)), rep("Wild", ncol(merged))),
                  state = rep(colnames(merged),2),
                  count = c(merged[1,], merged[2,]),
                  prop = c(merged[1,]/sum(merged[1,]), merged[2,]/sum(merged[2,]))
)

plotdf <- rbind(plotdf, tmp)


###### MeristemFunction #####

MeristemFunction <- trait %>% filter(Trait == "MeristemFunctioning", VALUE != "NoData")
data <- MeristemFunction %>% filter(MeristemFunction$wild.cultivated == "Cultivated")
data <- data %>% filter(data$Neither == 0)

table(data$CropState)
table(data$WR1State)
# table(data$WR2State)

# now we need to extract the relevant transitions
data <- data %>% filter(data$Both == 0)

# drop Axis 1 Indeterminate only present in a single crop
data <- data[-which(data$CropState == "Axis 1 Indeterminate"),]
# drop Determinate or Indeterminate only present in a single crop
data <- data[-which(data$CropState == "Determinate or Indeterminate"),]

tab1 <- table(data$CropState)
tab1 <- tab1[2:length(tab1)]
tab2 <- table(data$WR1State)
tab2 <- tab2[2:length(tab2)]

merged <- t(Merge2Tables(tab1, tab2, direction = "S"))
rownames(merged) <- c("Crop", "Wild")
merged[is.na(merged)] <- 0
chisq.test(merged, simulate.p.value = T)

# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# 
# data:  merged
# X-squared = 1.4961, df = NA, p-value = 0.3108

tmp <- data.frame(trait = rep("MeristemFunctioning", 2*ncol(merged)),
                  transition = c(rep("Crop", ncol(merged)), rep("Wild", ncol(merged))),
                  state = rep(colnames(merged),2),
                  count = c(merged[1,], merged[2,]),
                  prop = c(merged[1,]/sum(merged[1,]), merged[2,]/sum(merged[2,]))
)

plotdf <- rbind(plotdf, tmp)


###### BranchingMechanism #####

BranchingMechanism <- trait %>% filter(Trait == "BranchingMechanism", VALUE != "NoData")
data <- BranchingMechanism %>% filter(BranchingMechanism$wild.cultivated == "Cultivated")
data <- data %>% filter(data$Neither == 0)

table(data$CropState)
table(data$WR1State)
# table(data$WR2State)

# now we need to extract the relevant transitions
data <- data %>% filter(data$Both == 0)

tab1 <- table(data$CropState)
tab1 <- tab1[2:length(tab1)]
tab2 <- table(data$WR1State)
tab2 <- tab2[2:length(tab2)]

merged <- t(Merge2Tables(tab1, tab2, direction = "S"))
rownames(merged) <- c("Crop", "Wild")
merged[is.na(merged)] <- 0
chisq.test(merged, simulate.p.value = T)

# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# 
# data:  merged
# X-squared = 2.0449, df = NA, p-value = 0.6257

tmp <- data.frame(trait = rep("BranchingMechanism", 2*ncol(merged)),
                  transition = c(rep("Crop", ncol(merged)), rep("Wild", ncol(merged))),
                  state = rep(colnames(merged),2),
                  count = c(merged[1,], merged[2,]),
                  prop = c(merged[1,]/sum(merged[1,]), merged[2,]/sum(merged[2,]))
)

plotdf <- rbind(plotdf, tmp)

###### ShortShoots #####

ShortShoots <- trait %>% filter(Trait == "ShortShoots", VALUE != "NoData")
data <- ShortShoots %>% filter(ShortShoots$wild.cultivated == "Cultivated")
data <- data %>% filter(data$Neither == 0)

table(data$CropState)
table(data$WR1State)
# table(data$WR2State)

# now we need to extract the relevant transitions
data <- data %>% filter(data$Both == 0)

tab1 <- table(data$CropState)
tab1 <- tab1[2:length(tab1)]
tab2 <- table(data$WR1State)
tab2 <- tab2[2:length(tab2)]

merged <- t(Merge2Tables(tab1, tab2, direction = "S"))
rownames(merged) <- c("Crop", "Wild")
merged[is.na(merged)] <- 0
chisq.test(merged, simulate.p.value = T)

# Pearson's Chi-squared test with simulated p-value (based on 2000 replicates)
# 
# data:  merged
# X-squared = 0.12454, df = NA, p-value = 1

tmp <- data.frame(trait = rep("ShortShoots", 2*ncol(merged)),
                  transition = c(rep("Crop", ncol(merged)), rep("Wild", ncol(merged))),
                  state = rep(colnames(merged),2),
                  count = c(merged[1,], merged[2,]),
                  prop = c(merged[1,]/sum(merged[1,]), merged[2,]/sum(merged[2,]))
)

plotdf <- rbind(plotdf, tmp)

###### Plotting #####

plotdf$trait <- as.factor(plotdf$trait)

# difference in proportions version, uncomment for

# plotdf %>%
#   group_by(trait, state) %>%
#   reframe(diff = prop[transition == "Crop"] - prop[transition == "Wild"],
#           crop_n = count[transition == "Crop"],
#           wild_n = count[transition == "Wild"]) -> plotdf

# plotdf$trait <- factor(plotdf$trait, levels = c("BranchingPosition", "GrowthDirection", "BranchingMechanism",
#                                                 "FloweringAxis", "Phyllotaxis", "MeristemFunctioning", "BranchingType",
#                                                 "ShortShoots"))
# 
# p <- ggplot(plotdf, aes(x = state, y = diff)) + 
#   geom_bar(position = position_dodge(width = 0.8), stat = "identity", width = 0.7,
#            aes(fill = ifelse(diff > 0, "Crop", "Wild"))) + 
#   ggsci::scale_fill_aaas(name="Transition") + 
#   scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) + 
#   coord_flip() + 
#   facet_grid(trait ~ ., scales = "free", space="free", switch = "x") +
#   theme_bw() + 
#   ylim(-.5, .5) + 
#   labs(x = "Terminal state",
#        y = "Difference in proportion of transitions resulting in state") +
#   geom_hline(yintercept=0, linetype="dashed") +
#   geom_text(aes(label = crop_n, colour = ifelse(diff > 0, "white", "black")), y=0.01, hjust="left",
#             show.legend = F) +
#   geom_text(aes(label = wild_n,  colour = ifelse(diff < 0, "white", "black")),y=-0.01, hjust="right",
#             show.legend = F) + 
#   scale_color_manual(values = c("black", "white"))

# mirrored proportions version

plotdf$trait <- as.factor(plotdf$trait)

plotdf %>%
  group_by(trait, state) %>%
  reframe(prop_c = prop[transition == "Crop"], prop_w = -prop[transition == "Wild"],
          crop_n = count[transition == "Crop"],
          wild_n = count[transition == "Wild"]) -> plotdf

plotdf$trait <- factor(plotdf$trait, levels = c("BranchingPosition", "GrowthDirection", "BranchingMechanism",
                                                "FloweringAxis", "Phyllotaxis", "MeristemFunctioning", "BranchingType",
                                                "ShortShoots"))

p <- ggplot(plotdf, aes(x = state)) + 
    geom_bar(position = position_dodge(width = 0.8), stat = "identity", width = 0.7,
             aes(y = prop_c, fill = ifelse(prop_c >= 0, "Crop", "Wild"))) +
    geom_bar(position = position_dodge(width = 0.8), stat = "identity", width = 0.7,
           aes(y = prop_w, fill = ifelse(prop_w <= 0, "Wild", "Crop"))) + 
    ggsci::scale_fill_aaas(name="Transition") +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) +
    coord_flip() +
    facet_grid(trait ~ ., scales = "free", space="free", switch = "x") +
    theme_bw() +
    ylim(-.75, .75) +
    labs(x = "Terminal state",
         y = "Proportion of transitions resulting in state") +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_text(aes(label = crop_n, colour = "white"), y=0.01, hjust="left",
              show.legend = F) +
    geom_text(aes(label = wild_n, colour = "white"),y=-0.01, hjust="right",
              show.legend = F) +
    scale_colour_manual(values = c("white"))

p
