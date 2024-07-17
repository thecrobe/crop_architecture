## Join  together cultivation data with age
cultivation<-(read.csv(file = "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/LongFormat_croparchitecture_V2.csv", header=T)) 

# Load necessary library
library(dplyr)

# Assuming 'df' is your dataframe
# Calculate the total number of unique plant species - 1/3 of total database used in analysis
total_unique_species <- cultivation %>%
  distinct(Species_._ssp_._var) %>%
  nrow()


# Count unique species per trait
species_count_per_trait <- cultivation %>%
  filter(VALUE != "NoData") %>%  # Exclude rows where VALUE is "NoData"
  group_by(Trait) %>% 
  summarise(Species_Count = n_distinct(Species_._ssp_._var))

# Print the results
print("Count of unique plant species per trait:")
print(species_count_per_trait)

##### 
cultivation<-(read.csv(file = "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/LongFormat_croparchitecture_V2.csv", header=T))

# Define the cleanup function
clean_value_function <- function(values) {
  values <- tolower(values)  # Convert all to lowercase
  values <- trimws(values)  # Remove leading/trailing whitespace
  values <- gsub("terminodatal", "terminal", values)  # Correct specific errors
  values <- gsub("indeterminodatate", "indeterminate", values)
  values <- gsub("determinodatate", "determinate", values)
  values <- gsub(" {2,}", " ", values)  # Replace multiple spaces
  values <- gsub("\\s+and\\s+", " & ", values)  # Standardize 'and' to '&'
  values <- gsub("/", " & ", values)  # Replace '/' with ' & '
  values <- gsub("spipral", "spiral", values)  # Correct typos
  values <- gsub("rhytmic", "rhythmic", values)  # Correct typos
  return(values)
}

# Clean up typos
clean_value_function <- function(x) {
  x <- gsub("axis 1 orthtropic axis 2 plagiotropic", "axis 1 - orthotropic axis 2 - plagiotropic", x, fixed = TRUE)
  x <- gsub("axis 1 orthotropic", "axis 1 - orthotropic", x)
  x <- gsub("axis 1- orthotropic", "axis 1 - orthotropic", x)
  x <- gsub("axis 1 - spiral axis 2 - distichous", "axis 1 - spiral axis 2 - distichous", x, fixed = TRUE)
  x <- gsub("\\bspira\\b(?!l)", "spiral", x, perl = TRUE)
  x <- gsub("axis 1 monopodial axis 2 sympodial", "axis 1 - monopodial axis 2 - sympodial", x, fixed = TRUE)
  x <- gsub("acrotonic & basitonic", "basitonic & acrotonic", x, fixed = TRUE)
  x <- gsub("basitonic & mesotonic", "mesotonic & basitonic", x, fixed = TRUE)
  x <- gsub("Monopodial +Axis 1", "Axis 1 - Monopodial", x, ignore.case = TRUE)
  x <- gsub("Sympodial +Axis 2", "Axis 2 - Sympodial", x, ignore.case = TRUE)
  x <- trimws(x) # Remove any spaces at the end of each entry
  x <- tolower(x) # Convert all characters to lowercase
  return(x)
}
# Apply the function to multiple columns in the 'cultivation' dataframe
cultivation <- cultivation %>%
  mutate(VALUE = clean_value_function(VALUE)) %>%
  mutate(Direction.WR2.to.Crop = clean_value_function(Direction.WR2.to.Crop)) %>%
  mutate(Direction.WR2.to.CWR1 = clean_value_function(Direction.WR2.to.CWR1))

# Optionally, print the updated dataframe to verify changes
print(cultivation)

# Get the list of unique traits from the cultivation dataset
unique_traits <- unique(cultivation$Trait)

# Initialize a list to store the combined data for each trait
all_traits_data <- list()

# Loop through each trait
for (trait in unique_traits) {
  # Filter data for the current trait where conditions are met
  test <- cultivation %>%
    filter(Trait == trait) %>%
    filter(WR2toCrop == 1 & WR2toWR1 == 0) %>%
    dplyr::select(Direction.WR2.to.Crop) %>%
    filter(!is.na(Direction.WR2.to.Crop), Direction.WR2.to.Crop != "") %>%
    # Handle rows that do not split correctly
    separate(Direction.WR2.to.Crop, into = c("Before", "After"), sep = " to ", extra = "drop", fill = "right") %>%
    # Drop rows where 'Before' or 'After' is NA
    filter(!is.na(Before), !is.na(After))
  
  # Calculate the proportion of each 'Before' value
  before_proportions <- test %>%
    count(Before) %>%
    mutate(Proportion = n / sum(n) * 100,
           State = "Before",
           Trait = trait) %>%
   dplyr:: select(Trait, State, Before, Proportion)
  
  # Calculate the proportion of each 'After' value
  after_proportions <- test %>%
    count(After) %>%
    mutate(Proportion = n / sum(n) * 100,
           State = "After",
           Trait = trait) %>%
    dplyr::select(Trait, State, Before = After, Proportion)
  
  # Combine the proportions and store in the list
  combined_proportions <- bind_rows(before_proportions, after_proportions)
  all_traits_data[[trait]] <- combined_proportions
}

# Combine all traits data into one DataFrame
final_data <- bind_rows(all_traits_data, .id = "Trait") %>%
  mutate(Before = gsub("axis 1 monopodial axis 2 sympodial", "axis 1 - monopodial and axis 2 - sympodial", Before, fixed = TRUE))

str(final_data)

write.csv(final_data, "/Users/justinstewart/Dropbox/Collaborations/TobyKiers/CropArchitecture/CropArchitectureDomestication/Data/NewDatabase_Analyszed_arch/TraitProportions.csv", row.names = FALSE)

ggplot(final_data, aes(x = State, y = Proportion, fill = Before)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Trait, scales = "free_y") + # Faceting by Trait with free y scales if needed
  scale_fill_manual(values = rainbow(length(unique(final_data$Before)))) + # Setting colors
  labs(x = "State", y = "Proportion", fill = "Before") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust text angle for clarity
