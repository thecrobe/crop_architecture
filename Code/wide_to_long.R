library("ape")
library("readxl")
library("stringr")
library("dplyr")
library("data.table")

setwd("~/Dropbox/other_projects/crop_architecture/crop_architecture/")

# load wide
wide <- read_excel("Data/Microp_plant species_20-05-2022_GC_NWH_20250727_GC_NWH.xlsx")

# form long
long <- melt(setDT(wide), id.vars = c(1:5), measure.vars = c(9:19),
             variable.name = "Trait", value.name = "VALUE")

# reorder to match previous long
setorder(long, "Species.ssp.var")

# replace spaces in species names
long[, "Species.ssp.var" := gsub(" ", "_", Species.ssp.var)]

# remove suffix from TRAIT column
long[, "Trait" := gsub("_CLEAN", "", Trait)]

# Make sure all common names have consistent formatting
long[, Common.name := str_to_sentence(Common.name)]

# 20250727 make sure all VALUE have consistent formatting
long[, VALUE := str_to_sentence(VALUE)]

# change NAs to 'NoData'
long[is.na(get("VALUE")), ("VALUE"):="NoData"]

# 20250707 make sure all Wild/Cultivated are sentence case
long[, Wild.cultivated := str_to_sentence(Wild.cultivated)]

# remove trailing spaces on common names
long[, Common.name := trimws(Common.name)]

# 20250727 remove trailing whitespace on traits and values
long[, Trait := trimws(Trait)]
long[, VALUE := trimws(VALUE)]

# for Wild/cultivated, we need to make them 'Wild' if field Common.name contains
# 'wild relative'
long[Wild.cultivated == "Wild/cultivated" & Common.name %flike% "wild relative", Wild.cultivated := "Wild"]
long[Wild.cultivated == "Wild/cultivated", Wild.cultivated := "Cultivated"]

# read in tree
tree <- read.tree("Output/tree_pgls_base_clean.tre")

# get tip labels
tips.in.tree <- tree$tip.label

# get species in db
tips.in.db <- long$Species.ssp.var[which(long$Wild.cultivated == "Cultivated")]

# find number of species in db that are not in tip labels
length(unique(tips.in.db[!tips.in.db %in% tips.in.tree]))
# [1] 99

# inspecting these reveals a few missing - concerning - but a few are just name mismatches
# e.g. Vasconcellea_x_pentagona != Vasconcellea_pentagona
# Solanum_lycopersicum__Heinz_1706 != Solanum_lycopersicum_Heinz_1706
# Citrus_x_taitensis != Citrus_taitensis
# so let's fix these 

# replace _x_ with _
long[, "Species.ssp.var" := gsub("_x_", "_", Species.ssp.var)]
# some Xs are actually ×, i.e. the multiplication symbol :/
long[, "Species.ssp.var" := gsub("_×_", "_", Species.ssp.var)]
# replace __ with _
long[, "Species.ssp.var" := gsub("__", "_", Species.ssp.var)]

# check matching again

# get species in db
tips.in.db <- long$Species.ssp.var[which(long$Wild.cultivated == "Cultivated")]

# find number of species in db that are not in tip labels
length(unique(tips.in.db[!tips.in.db %in% tips.in.tree]))
# [1] 85
unique(tips.in.db[!tips.in.db %in% tips.in.tree])

length(unique(tips.in.tree[!tips.in.tree %in% tips.in.db]))
# [1] 2
unique(tips.in.tree[!tips.in.tree %in% tips.in.db])
# [1] "Vigna_unguiculata_Cowpea"   "Vigna_unguiculata_aparagus"
# we'll switch these in the db
long[Common.name == "Cowpea", Species.ssp.var := "Vigna_unguiculata_Cowpea"]
long[Common.name == "Asparagus bean", Species.ssp.var := "Vigna_unguiculata_aparagus"]

# # European raspberry has no wild relatives in the database, so let's remove
# long <- long[Common.name != "European raspberry"]
# # ditto American raspberry
# long <- long[Common.name != "American raspberry"]
# # ditto ajwain

# for ease, we will filter everything with no wild relative
for (i in unique(long[Wild.cultivated == "Cultivated", Common.name])) {
  if (nrow(long[Common.name %flike% paste(i, "wild relative")]) == 0) {
    long <- long[Common.name != i]
  } 
}

# 20250727 figure out things with multiple wild relatives
cultivated <- long[Wild.cultivated == "Cultivated"]
for (i in cultivated$Common.name) {
  sub <- long[Common.name %flike% paste(i, "wild relative")]
  sub <- sub %>% group_by(Common.name)
  wr <- str_sub(sub$Common.name, -1)
  if ("3" %in% wr | "4" %in% wr) print(i)
}
# fix these ex situ

# let's prepare columns to hold the change coding
long$WR2toCrop <- rep(0, nrow(long))
long$WR2toWR1 <- rep(0, nrow(long))
long$Both <- rep(0, nrow(long))
long$Neither <- rep(1, nrow(long))
long$Direction.WR2toCrop <- rep("", nrow(long))
long$Direction.WR2toWR1 <- rep("", nrow(long))

# # to do the change coding it's easiest to do a loop
# # probably those smarter and better at R than me can vectorise but...
# for (i in unique(long[Wild.cultivated == "Cultivated", Common.name])) {
#   query1 <- paste(i, "wild relative 1")
#   query2 <- paste(i, "wild relative 2")
#   # sometimes we have wild relative 2 and 3
#   if (nrow(long[Common.name == query1]) == 0) {
#     query1 <- paste(i, "wild relative 2")
#     query2 <- paste(i, "wild relative 3")
#   }
#   # sometimes we have wild relative 1 and 3
#   if (nrow(long[Common.name == query2]) == 0) {
#     query1 <- paste(i, "wild relative 1")
#     query2 <- paste(i, "wild relative 3")
#   }
#   
#   for (j in unique(long[,Trait])) {
#     cult.val <- long[Common.name == i & Trait == j, VALUE]
#     wr1.val <- long[Common.name == query1 & Trait == j, VALUE]
#     wr2.val <- long[Common.name == query2 & Trait == j, VALUE]
#     if (wr2.val != cult.val & wr2.val != wr1.val) {
#       long[Common.name == i & Trait == j, Both := 1]
#     } else if (wr2.val != cult.val) {
#       long[Common.name == i & Trait == j, WR2toCrop := 1]
#     } else if (wr2.val != wr1.val) {
#       long[Common.name == i & Trait == j, WR2toWR1 := 1]
#     } else {
#       long[Common.name == i & Trait == j, Neither := 1]
#     }
#   }
# }

# turns out we avoid some problems of shared common names this way
for (i in 1:nrow(long)) {
  if (long[i, Wild.cultivated == "Cultivated"]) {
    sp <- long[i, Species.ssp.var]
    common <- long[i, Common.name]
    trait <- long[i, Trait]
    cult.val <- long[i, VALUE]
    query1 <- paste(common, "wild relative 1")
    query2 <- paste(common, "wild relative 2")
    # # sometimes we have wild relative 2 and 3
    # if (nrow(long[Common.name == query1]) == 0) {
    #   query1 <- paste(common, "wild relative 2")
    #   query2 <- paste(common, "wild relative 3")
    # }
    # # sometimes we have wild relative 1 and 3
    # if (nrow(long[Common.name == query2]) == 0) {
    #   query1 <- paste(common, "wild relative 1")
    #   query2 <- paste(common, "wild relative 3")
    # }
    wr1.val <- long[Common.name == query1 & Trait == trait,
                    VALUE]
    wr2.val <- long[Common.name == query2 & Trait == trait,
                    VALUE]
    # we skip anything which can't inform the comparison
    if (any(c(cult.val, wr1.val, wr2.val) == "NoData")) next
    # we test for differences
    if (wr2.val != cult.val & wr2.val != wr1.val) {
      long[i, WR2toCrop := 1] # this is to match Justin's coding
      long[i, WR2toWR1 := 1] # ditto
      long[i, Both := 1]
      long[i, Neither := 0]
      long[i, Direction.WR2toCrop := paste(wr2.val, cult.val, sep = " to ")]
      long[i, Direction.WR2toWR1 := paste(wr2.val, wr1.val, sep = " to ")]
    } else if (wr2.val != cult.val) {
      long[i, WR2toCrop := 1]
      long[i, Neither := 0]
      long[i, Direction.WR2toCrop := paste(wr2.val, cult.val, sep = " to ")]
    } else if (wr2.val != wr1.val) {
      long[i, WR2toWR1 := 1]
      long[i, Neither := 0]
      long[i, Direction.WR2toWR1 := paste(wr2.val, wr1.val, sep = " to ")]
    }
  }
}

colnames(long) <- c("Family", "Species.ssp.var", "common.name", "wild.cultivated",
                    "Herbaceous.woody", "Trait", "VALUE", "WR2toCrop", "WR2toWR1",
                    "Both", "Neither", "Direction.WR2toCrop", "Direction.WR2toWR1")


write.csv(long, "Data/20250727_LongFormat_croparchitecture_V4.csv", row.names = F)
# unique(long[Trait == "GrowthDirection", VALUE])
