# need development terra for horizontal continuous legend
# remotes::install_github("rspatial/terra")

library(tidyverse)
library(reshape2)
library(brms)
library(sf)
library(sp)
# library(raster)
library(terra)
library(ggpp)
library(ggplot2)

## Join  together cultivation data with age
cultivation<-(read.csv(file = "~/Dropbox/other_projects/crop_architecture/crop_architecture/Data/20250211_LongFormat_croparchitecture_V3.csv", header=T)) %>% filter(wild.cultivated == "Cultivated")
age<-read.csv(file = "~/Dropbox/other_projects/crop_architecture/crop_architecture/Data/milla_crop_origins_11march2020.csv", header=T)
age$Species.ssp.var<-age$species_name_
join<-inner_join(age,cultivation) %>% 
  group_by("species_name")

merged<-join %>% 
  filter(wild.cultivated == "Cultivated") %>%
  filter(!is.na(mode_ecoreg_centroid_lon) & !is.na(mode_ecoreg_centroid_lat) & !is.na(minimum_time_cultivation))
# Create a SpatialPoints object with latitude and longitude coordinates
points <- SpatialPointsDataFrame(coords = cbind(merged$mode_ecoreg_centroid_lon,merged$mode_ecoreg_centroid_lat), data = merged, proj4string = CRS("+proj=longlat +datum=WGS84"))
plot(points)

# Load your raster data
bio1 <- rast("~/Dropbox/other_projects/crop_architecture/crop_architecture/Data/CHELSA_bio1_1981-2010_V.2.1.tif")
bio12 <- rast("~/Dropbox/other_projects/crop_architecture/crop_architecture/Data/CHELSA_bio12_1981-2010_V.2.1.tif")
# bio1[bio1 == 0] <- NA # Set oceans/ice as NA


bio1_extract <- terra::extract(bio1, vect(points))
bio12_extract <- terra::extract(bio12, vect(points))
merged$MAT<-bio1_extract$`CHELSA_bio1_1981-2010_V.2.1`
merged$MAP<-bio12_extract$`CHELSA_bio12_1981-2010_V.2.1`
merged$cultivation<-as.numeric(merged$minimum_time_cultivation)
plot(bio1)

# here we need to try and work with bio1 to hack it into a projection

##### Map plot #####
library(rnaturalearth); library(rnaturalearthdata)
# library(rgdal) # depends rgdal but retired at end 2023 - install 1.6-7 from source https://cran.r-project.org/src/contrib/Archive/rgdal/
# install.packages('~/Downloads/rgdal_1.6-7.tar.gz', repos = NULL, type = 'source')

# Define Robinson projection
# crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
# these are the same, and also ESRI:53030
# robin.crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
crs <- "ESRI:53030"

# Create worldmap for grey background
Worldmap <- ne_countries(scale="medium", returnclass = "sv")
Worldmap <- subset(Worldmap, Worldmap$name != "Antarctica")
Worldmap <- project(Worldmap, crs)

ext(bio1) <- ext(-180, 180, -88, 88)
bio1.rob <- project(bio1, crs, method="average") # any losses don't matter because it's just for plotting
# ext(bio1.rob) <- c(-180, 180, -88, 88)
bio1.mask <- mask(bio1.rob, Worldmap)
# ext(bio1.mask) <- ext(-180, 180, -88, 88)
svg("~/Dropbox/other_projects/crop_architecture/figures/CropArchitectureMap.svg", width = 8.85, height = 5.24)
plot(bio1.mask, axes=F, buffer=F, box=F, clip=F, plg=list(title = "Mean annual temperature (°C)",
                                                          title.x = 7e+06, title.y = -5e+06,
                                                          title.cex=0.7, cex=0.7,
                                                          x=c(3e+06, 1.1e+07), y=c(-5.5e+06, -5e+06), horiz=T,
                                                          tic="out"),
     col=viridis::viridis(n=500,option = "magma",direction = 1), mar=.2, maxcell=1e8)
points <- project(vect(points), crs)
points(points, cex=0.5)
add_legend(-1e+06, -4.3e+06, legend="Crop origin", pch=20, col="black", cex=0.7, bty="n")
dev.off()

# Define lat/long lines for plotting
# graticule = st_graticule(lat=seq(-90, 90, 30), lon=seq(-180, 180, 30)) # Create graticules (lat/long lines on map)
# graticule = st_transform(graticule, crs) # Convert to robinson projection
# 
# points <- SpatialPointsDataFrame(coords = cbind(merged$mode_ecoreg_centroid_lon,merged$mode_ecoreg_centroid_lat), data = merged, proj4string = CRS("+proj=longlat +datum=WGS84"))
# points <- st_as_sf(points, coords = c("longitude", "latitude"), crs =4326)
# points.rob <- st_transform(points, crs = crs)
# plot(points.rob)
# length(unique(points$common_name_crop))
# extent(bio1)<- c(-180, 180, -88, 88) # Adjust extreme edges of extent for re-projection
# 
# # Export high res plot with grey background, lat/long graticule lines
# png('~/Dropbox/other_projects/crop_architecture/figures/CropArchitectureMap.png', width=(ncol(bio1)/10), height=(nrow(bio1)/10))
# par(bg=NA) # remove plot background
# Worldmap <- st_as_sf(Worldmap)
# Worldmap <- st_transform(Worldmap, crs = st_crs(bio1))
# raster::plot(bio1,  col=viridis::viridis(n=500,option = "magma",direction = 1), maxpixels=1e8, bg=NA, add=F, colNA="white")
# points <- st_as_sf(points)
# points <- st_transform(points, crs = st_crs(bio1))
# plot(points$geometry,col="black",cex=3.5,pch=19,add=T)
# #plot(graticule$geometry, col="black", lty = 5, lwd=5)
# dev.off()

###### Iterate through trait values to test for cultivation time and climate
unique_traits <- unique(merged$Trait) 
unique_traits <- unique_traits[unique_traits != "Branching" & unique_traits != "Reiteration"
                               & unique_traits != "BasalBranching"]


# Initialize an empty dataframe to store results
final_results <- data.frame()
hypothesis_results_list <- list()

# Iterate over the selected traits
for (selected_trait in unique_traits) {
  # Filter data for the selected trait
  trait_data <- merged %>% filter(Trait == selected_trait)
  
  # Fit the logistic regression model for WR2toCrop
  model_crop <- brm(
    formula = as.factor(WR2toCrop) ~ MAP + MAT + cultivation,
    data = trait_data,
    family = bernoulli(link = "logit"),
    warmup = 500,
    iter = 2000,
    chains = 4,
    cores = 4,
    seed = 123
  )
  
  # Fit the logistic regression model for WR2toWR1
  model_wild <- brm(
    formula = as.factor(WR2toWR1) ~ MAP + MAT + cultivation,
    data = trait_data,
    family = bernoulli(link = "logit"),
    warmup = 500,
    iter = 2000,
    chains = 4,
    cores = 4,
    seed = 123
  )
  
  # Extract coefficients, variable names, and estimates with errors for the crop model
  crop_summary <- summary(model_crop)
  
  # Extract coefficients, variable names, and estimates with errors for the wild model
  wild_summary <- summary(model_wild)
  
  # Access coefficients and standard errors for both models
  crop_coef <- as.data.frame(crop_summary$fixed)
  crop_coef$Variable <- row.names(crop_coef)
  wild_coef <- as.data.frame(wild_summary$fixed)
  wild_coef$Variable <- row.names(wild_coef)
  
  # Add the trait name to the dataframes
  crop_coef$Trait <- selected_trait
  wild_coef$Trait <- selected_trait
  
  # # Rename columns for clarity
  # colnames(crop_coef) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Rhat", "Bulk_ESS", "Tail_ESS", "row.names.crop_coef.", "Trait")
  # colnames(wild_coef) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Rhat", "Bulk_ESS", "Tail_ESS", "row.names.wild_coef.", "Trait")
  # 
  # # these are already the column names
  
  # Create dataframes and run hypothesis tests
  crop <- data.frame(model_crop)
  wild <- data.frame(model_wild)
  
  # Perform hypothesis tests for MAP and MAT
  test_map <- data.frame(cbind(crop$b_MAP, wild$b_MAP))
  hypothesis_map <- hypothesis(test_map, hypothesis = "X1 > X2")
  hypothesis_map <- hypothesis_map$hypothesis
  
  test_mat <- data.frame(cbind(crop$b_MAT, wild$b_MAT))
  hypothesis_mat <- hypothesis(test_mat, hypothesis = "X1 > X2")
  hypothesis_mat <- hypothesis_mat$hypothesis
  
  test_cultivation <- data.frame(cbind(crop$b_cultivation, wild$b_cultivation))
  hypothesis_cultivation <- hypothesis(test_cultivation, hypothesis = "X1 > X2")
  hypothesis_cultivation <- hypothesis_cultivation$hypothesis
  
  # Add the trait name and variable (MAP or MAT) to the hypothesis test results
  hypothesis_map$Trait <- selected_trait
  hypothesis_map$Variable <- "MAP"
  
  hypothesis_mat$Trait <- selected_trait
  hypothesis_mat$Variable <- "MAT"
  
  hypothesis_cultivation$Trait <- selected_trait
  hypothesis_cultivation$Variable <- "Cultivation"
  
  # Combine hypothesis test results into a single dataframe
  hypothesis_dataframe <- rbind(hypothesis_map, hypothesis_mat,hypothesis_cultivation )
  # Append the hypothesis test results for this trait to the list
  hypothesis_results_list[[selected_trait]] <- rbind(hypothesis_map, hypothesis_mat,hypothesis_cultivation)
  
  # Combine results into a single dataframe for this trait
  # results_dataframe <- bind_rows(
  #   crop_coef %>% mutate(Model = "Crop"),
  #   wild_coef %>% mutate(Model = "Wild")
  # )
  results_dataframe <- rbind(crop_coef, wild_coef)
  rownames(results_dataframe) <- NULL
  results_dataframe$Model <- c(rep("Crop",4), rep("Wild",4))
  
  # Append the results for this trait to the final results
  final_results <- rbind(final_results, results_dataframe)
  hypothesis_dataframe <- do.call(rbind, hypothesis_results_list)
  
}

# Optionally, you can save the final results to a CSV file
write.csv(final_results, "~/Dropbox/other_projects/crop_architecture/crop_architecture/Output/20250415_MAT_MAP_hypothesis_testing.csv",
          row.names = FALSE)

hypothesis_dataframe <- hypothesis_dataframe %>%
  mutate(
    Trait = case_when(
      Trait == "BranchingPosition" ~ "Branching position",
      Trait == "BranchingType" ~ "Branching type",
      Trait == "Growthdirection" ~ "Growth direction",
      Trait == "FloweringAxis" ~ "Flowering axis",
      Trait == "MeristemFunction" ~ "Meristem function",
      Trait == "BranchingMechanism" ~ "Branching mechanism",
      Trait == "ShortShoots" ~ "Short shoots",
      TRUE ~ Trait
    ),
    Variable = case_when(
      Variable == "Cultivation" ~ "Antiquity",
      TRUE ~ Variable
    )
  )

# Plot the model results
hypothesis_dataframe %>% filter(Trait !="Branching") %>%
ggplot(aes(x = Trait, y = Variable, fill = Estimate)) + 
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "mako") +
  # Add tiles with white borders
  geom_text(aes(label = sprintf("%.3f", Post.Prob)), size = 3.5, vjust = 1.5, color="red") +  # Add post. prob. labels
  labs(title = "Coefficient estimates and significance probabilities", x = "Trait", y = "Variable", fill = "Estimate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> p  # Adjust the angle of x-axis labels for better visibility

ggsave("~/Dropbox/other_projects/crop_architecture/figures/20250415_Fig_SX_climate_antiquity.svg", p, width = 8, height =4)

##### Plot of cultivation time 

#Get average to remove duplicate points b/c multiple traits
average_cultivation_df <- merged %>%
  group_by(common_name_crop, species_name) %>%
  summarise(AverageCultivationTime = mean(minimum_time_cultivation, na.rm = TRUE))

# # Extractr random crops to show on plot highlight
# random_crops <- average_cultivation_df[sample(nrow(average_cultivation_df), 5), ]
# 
# selected_rows <- merged[merged$common_name_crop %in% random_crops$common_name_crop, ] %>% 
#   dplyr::select(species_name, minimum_time_cultivation ) %>% group_by(species_name) %>%
#   summarise(AverageCultivationTime = mean(minimum_time_cultivation, na.rm = TRUE))
# 
# 
# 
# tt<-merged %>% filter(common_name_crop == "Lard_seed")
# 
# tt$species_name
# oldest_crop <- average_cultivation_df$common_name_crop[which.max(average_cultivation_df$AverageCultivationTime)]
# most_recent_crop <- average_cultivation_df$common_name_crop[which.min(average_cultivation_df$AverageCultivationTime)]
 
crops2label <- c("Turnip", "Bay_laurel", "Sunflower", "Almond", "White_yam")
names(crops2label) <- c("Brassica rapa", "Laurus nobilis", "Helianthus annuus",
                        "Prunus dulcis", "Dioscorea caynennensis")

plot<-average_cultivation_df %>% 
  ggplot(aes(x=(AverageCultivationTime/1000), y="")) +
  geom_point(position = position_jitter(width = 0.3, seed = 2),
             color = ifelse(average_cultivation_df$common_name_crop %in% crops2label,
                            "red", "black")) + 
  ggrepel::geom_label_repel(aes(label = ifelse(common_name_crop %in% crops2label,
                                species_name, "")),
                            position = ggpp::position_jitternudge(width = 0.3, seed = 2,
                                                                  y=1.3),
                            box.padding = 1.2, fontface = "italic",
                            segment.color = "red") +
  xlab("Cultivation age (1000s of years) ") +
  ylab("") + 
  theme_minimal()

plot

ggsave(plot = plot,filename = "~/Dropbox/other_projects/crop_architecture/figures/plot_cult.svg",
       width=9, height=2)
