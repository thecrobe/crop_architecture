library("ape")
library("phytools")
library("WorldFlora")
here::i_am("Code/phylo_prep_for_pgls.R")

tree <- read.tree(here::here("Data/ALLMB.tre"))
tree <- multi2di(tree)

full.data <- read.csv(here::here("Data/20250727_LongFormat_croparchitecture_V4.csv"))

crops <- full.data[full.data$wild.cultivated == "Cultivated", ]

crops$Crop <- gsub("[ \t]+$", "", crops$Species.ssp.var, perl = TRUE)
crops$Crop <- gsub(" ", "_", crops$Crop)
crops$genus <- sapply(crops$Crop, function(x){strsplit(x, split = "_", fixed = TRUE)[[1]][1]})
crops$Family <- gsub("[ \t]+$", "", crops$Family, perl = TRUE)

## Fixing families
## crops$Family[crops$Tree == "Cyrtosperma_merkusii"] <- "Araceae"
## crops$Family[grep("Breynia", crops$Tree)] <- "Phyllanthaceae"
## crops$Family[grep("Elaeagnus", crops$Tree)] <- "Elaeagnaceae"
## crops$Family[crops$Tree == "Elymus_caninus"] <- "Gramineae"
## crops$Family[crops$Tree == "Myrica_rubra"] <- "Myricaceae"
## crops$Family[crops$Tree == "Myristica_fragrans"] <- "Myristicaceae"
## ## Problems with Sapindaceae - Allophylus_cobbe -> could be Rhus_cobbe, Anacardiaceae
## crops$Family[grep("Vitis", crops$Tree)] <- "Vitaceae"

genera.in.phylo <- lapply(unique(crops$genus), function(x){grep(x, tree$tip.label)})
genera.not.in.phylo <- unique(crops$genus)[which(sapply(genera.in.phylo, length) == 0)]

## genus.conv <- data.frame(not.phylo = genera.not.in.phylo,
##                          in.phylo = c("Uvaria",
##                                       "Apium", ## sister genus
##                                       "Wahlenbergia",
##                                       "Polymnia",
##                                       "Ambrosia", ## sister genus, along with Xanthium
##                                       "Merremia",
##                                       "Sedum",
##                                       "Dupontia", ## sister group to Poa
##                                       "Bambusa",
##                                       "Myrciaria",
##                                       "Oxyria", ## one of sister genera, along with Rheum
##                                       "Chaenomeles",
##                                       "Symphonia", ## from sister group Symphonieae
##                                       "Cardiospermum", ## sister genus
##                                       "Chrysophyllum",
##                                       "Musella" ## Sister genus to Musa
##                                       )
##                          )                                      

genus.conv <- data.frame(not.phylo = genera.not.in.phylo,
                         in.phylo = c("Bambusa", ## Synonym to Fimbribambusa
                                      "Chrysophyllum" ## Synonym to Gambeya
                                      )
                         )

tree$tip.label <- gsub("Bambusa_horsfieldii", "Fimbribambusa_horsfieldii", tree$tip.label)
# tree$tip.label <- gsub("Chaenomeles_sinensis", "Pseudocydonia_sinensis", tree$tip.label)
tree$tip.label <- gsub("Chrysophyllum_africanum", "Gambeya_africana", tree$tip.label)
tree$tip.label <- gsub("Chrysophyllum_beguei", "Gambeya_beguei", tree$tip.label)
tree$tip.label <- gsub("Chrysophyllum_boukokoense", "Gambeya_boukokoensis", tree$tip.label)

## crops$genus[!is.na(match(crops$genus, genus.conv$not.phylo))] <- genus.conv$in.phylo[na.omit(match(crops$genus, genus.conv$not.phylo))]


## ## Reincluding Vigna_unguiculata after Guillaume's instructions

## crops$Crop[419] <- "Vigna_unguiculata_cowpea"
## crops$Crop[420] <- "Vigna_unguiculata_asparagus"

tips.to.keep <- sapply(lapply(paste0(sort(unique(crops$genus)), "_"), function(x){grep(x, tree$tip.label, fixed = TRUE)}), function(y){y[1]})

tree.pruned <- drop.tip(tree, tree$tip.label[-tips.to.keep])

tree.gls <- phytools::force.ultrametric(tree.pruned, "nnls")
tree.gls$node.label <- NULL

for(i in 1:nrow(crops)){
    print(i)
    tree.gls <- bind.tip(tree.gls, tip.label = crops$Species.ssp.var[i], where = grep(paste0(crops$genus[i], "_"), tree.gls$tip.label, fixed = TRUE)[1], edge.length = 0)
}

tree.final.tip <- drop.tip(tree.gls, tree.gls$tip.label[!(tree.gls$tip.label %in% crops$Species.ssp.var)])
tree.final.tip <- drop.tip(tree.final.tip, which(duplicated(tree.final.tip$tip.label) == TRUE))

write.tree(tree.final.tip, here::here("Output/tree_pgls_tip_clean.tre"))

## Polytomies at the base of the genera

tree.base <- tree.pruned

tips.to.keep.full <- c()
for(i in 1:length(unique(crops$genus))){
    tips.to.keep.full <- c(tips.to.keep.full, grep(paste0(unique(crops$genus)[i], "_"), tree.base$tip.label, fixed = TRUE)[1])
}
tips.to.keep.full <- unique(tips.to.keep.full)

tree.base <- drop.tip(tree.base, tree.base$tip.label[-tips.to.keep.full])
height <- max(nodeHeights(tree.base))

for(i in 1:nrow(crops)){
    print(i)
    genus <- grep(paste0(crops$genus[i], "_"), tree.base$tip.label, fixed = TRUE)
    if(length(genus) == 1){
        focal.node <- nodepath(tree.base, genus, getMRCA(tree.base, tree.base$tip.label))[2]
    } else {
      focal.node <- getMRCA(tree.base, tree.base$tip.label[genus])
    }
    tree.base <- bind.tip(tree.base, tip.label = crops$Species.ssp.var[i], where = focal.node, edge.length = height - nodeheight(tree.base, focal.node))
}

tree.final.base <- drop.tip(tree.base, tree.base$tip.label[!(tree.base$tip.label %in% crops$Species.ssp.var)])
tree.final.base <- drop.tip(tree.final.base, which(duplicated(tree.final.base$tip.label) == TRUE))

write.tree(tree.final.base, here::here("Output/tree_pgls_base_clean.tre"))
