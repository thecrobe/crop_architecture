library(tidyverse)
library(reshape2)
library(brms)
library(sf)
library(sp)
library(raster)
library(terra)

## Join  together cultivation data with age
cultivation<-(read.csv(file = "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/LongFormat_croparchitecture_V2.csv", header=T)) %>% filter(wild.cultivated == "Cultivated")
age<-read.csv(file = "/Users/justinstewart/Downloads/crop_origins_11march2020.csv", header=T)
age$Species_._ssp_._var<-age$species_name_
join<-inner_join(age,cultivation) %>% 
  group_by("species_name")

merged<-join %>% 
  filter(wild.cultivated == "Cultivated") %>%
  filter(!is.na(mode_ecoreg_centroid_lon) & !is.na(mode_ecoreg_centroid_lat) & !is.na(minimum_time_cultivation))
# Create a SpatialPoints object with latitude and longitude coordinates
points <- SpatialPointsDataFrame(coords = cbind(merged$mode_ecoreg_centroid_lon,merged$mode_ecoreg_centroid_lat), data = merged, proj4string = CRS("+proj=longlat +datum=WGS84"))
plot(points)

# Load your raster data
bio1 <- raster("/Users/justinstewart/Downloads/MAT_Masked.tif")
bio12 <- raster("/Users/justinstewart/CHELSA_bio8_1981-2010_V.2.1.tif")
#bio1[bio1 == 0] <- NA # Set oceans/ice as NA


bio1_extract <- raster::extract(bio1, points)
bio12_extract <- raster::extract(bio12, points)
merged$MAT<-(bio1_extract*.1)-273.15
merged$MAP<-bio12_extract*0.1
merged$cultivation<-as.numeric(merged$minimum_time_cultivation)
plot(bio1)



##### Map plot #####
library(sf)
library(rnaturalearth); library(rnaturalearthdata)
library(rgdal)
# Define Robinson projection
crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
robin.crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"

# Create worldmap for grey background
Worldmap <- ne_countries(scale="medium")
Worldmap <- as(Worldmap, "Spatial")
Worldmap.robin <- spTransform(Worldmap, CRSobj = crs)

# Define lat/long lines for plotting
graticule = st_graticule(lat=seq(-90, 90, 30), lon=seq(-180, 180, 30)) # Create graticules (lat/long lines on map)
graticule = st_transform(graticule, crs) # Convert to robinson projection

points <- SpatialPointsDataFrame(coords = cbind(merged$mode_ecoreg_centroid_lon,merged$mode_ecoreg_centroid_lat), data = merged, proj4string = CRS("+proj=longlat +datum=WGS84"))
points <- st_as_sf(points, coords = c("longitude", "latitude"), crs =4326)
points.rob <- st_transform(points, crs = crs)
plot(points.rob)
length(unique(points$common_name_crop))
ext(bio1) <- c(-180, 180, -88, 88) # Adjust extreme edges of extent for re-projection

# Export high res plot with grey background, lat/long graticule lines
png('/Users/justinstewart/Desktop/CropArchitectureMap.png', width=(ncol(bio1)), height=(nrow(bio1)))
par(bg=NA) # remove plot background
Worldmap <- st_as_sf(Worldmap)
Worldmap <- st_transform(Worldmap, crs = st_crs(bio1))
raster::plot(bio1,  col=viridis::viridis(n=500,option = "magma",direction = 1), maxpixels=1e8, bg=NA,add=TR)
points <- st_as_sf(points)
points <- st_transform(points, crs = st_crs(bio1))
plot(points$geometry,col="black",cex=35,pch=19,add=T)
#plot(graticule$geometry, col="black", lty = 5, lwd=5)
dev.off()

###### Iterate through trait values to test for cultivation time and climate
unique_traits <- unique(merged$Trait) 
unique_traits <- unique_traits[unique_traits != "Branching"]


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
  crop_coef$Variable <- data.frame(row.names(crop_coef))
  wild_coef <- as.data.frame(wild_summary$fixed)
  wild_coef$Variable <- data.frame(row.names(wild_coef))
  
  # Add the trait name to the dataframes
  crop_coef$Trait <- selected_trait
  wild_coef$Trait <- selected_trait
  
  # Rename columns for clarity
  colnames(crop_coef) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Rhat", "Bulk_ESS", "Tail_ESS", "row.names.crop_coef.", "Trait")
  colnames(wild_coef) <- c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Rhat", "Bulk_ESS", "Tail_ESS", "row.names.wild_coef.", "Trait")
  
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
  hypothesis_dataframe <- bind_rows(hypothesis_map, hypothesis_mat,hypothesis_cultivation )
  # Append the hypothesis test results for this trait to the list
  hypothesis_results_list[[selected_trait]] <- bind_rows(hypothesis_map, hypothesis_mat,hypothesis_cultivation)
  
  # Combine results into a single dataframe for this trait
  results_dataframe <- bind_rows(
    crop_coef %>% mutate(Model = "Crop"),
    wild_coef %>% mutate(Model = "Wild")
  )
  
  # Append the results for this trait to the final results
  final_results <- bind_rows(final_results, results_dataframe)
  hypothesis_dataframe <- do.call(rbind, hypothesis_results_list)
  
}

# Optionally, you can save the final results to a CSV file
write.csv(final_results, "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/results_crop_climate_antiquity.csv", row.names = FALSE)

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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust the angle of x-axis labels for better visibility

# Display the plot
print(p)

##### Plot of cultivation time 

#Get average to remove duplicate points b/c multiple traits
average_cultivation_df <- merged %>%
  group_by(common_name_crop) %>%
  summarise(AverageCultivationTime = mean(cultivation, na.rm = TRUE))

# Extractr random crops to show on plot highlight
random_crops <- average_cultivation_df[sample(nrow(average_cultivation_df), 5), ]

selected_rows <- merged[merged$common_name_crop %in% random_crops$common_name_crop, ] %>% 
  dplyr::select(species_name, cultivation ) %>% group_by(species_name) %>%
  summarise(AverageCultivationTime = mean(cultivation, na.rm = TRUE))



tt<-merged %>% filter(common_name_crop == "Lard_seed")

tt$species_name
oldest_crop <- average_cultivation_df$common_name_crop[which.max(average_cultivation_df$AverageCultivationTime)]
most_recent_crop <- average_cultivation_df$common_name_crop[which.min(average_cultivation_df$AverageCultivationTime)]
 
plot<-average_cultivation_df %>% 
  ggplot(aes(x=(AverageCultivationTime/1000), y="")) +
  geom_jitter(size=2, alpha=0.8) + 
  xlab("Cultivation age (per 1000 years ago) ") + 
  theme_minimal()

ggsave(plot = plot,filename = "/Users/justinstewart/Downloads/plot_cult.png",scale=0.6)
