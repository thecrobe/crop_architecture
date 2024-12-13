library("here")
library("tidyr")
here::i_am("Code/ProportionChangeArchitecture.R")

## Join  together cultivation data with age
cultivation <- (read.csv(file = here::here("Data/LongFormat_croparchitecture_V2.csv"), header=T))

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
cultivation <- (read.csv(file = here::here("Data/LongFormat_croparchitecture_V2.csv"), header=T))
cultivation$VALUE <- trimws(tolower(cultivation$VALUE))

# Clean up typos
clean_value_function <- function(x) {
  x <- tolower(x)  # Convert all to lowercase
  x <- trimws(x)  # Remove leading/trailing whitespace
  x <- gsub(" {2,}", " ", x)  # Replace multiple spaces
  x <- gsub("\\s+\\&\\s+", " and ", x)  # Standardize 'and' to '&'
  x <- gsub("^$", "nodata", x, perl = TRUE) # Replace "" with "nodata"
  x <- gsub("/ ", "/", x, fixed = TRUE)  # Replace '/ ' with '/'
  #x <- gsub("/", " & ", x)  # Replace '/' with ' & '
  }

# Apply the cleaning function to multiple columns in the 'cultivation' dataframe
cultivation <- cultivation %>%
  mutate(VALUE = clean_value_function(VALUE)) %>%
  mutate(Direction.WR2.to.Crop = clean_value_function(Direction.WR2.to.Crop)) %>%
  mutate(Direction.WR2.to.CWR1 = clean_value_function(Direction.WR2.to.CWR1))

conversion.Growthdirection <-
    data.frame(
        original = sort(unique(c(levels(as.factor(cultivation[cultivation$Trait == "Growthdirection", "Direction.WR2.to.Crop"])),
                            levels(as.factor(cultivation[cultivation$Trait == "Growthdirection", "Direction.WR2.to.CWR1"]))
                            ))),
        grouped = c("ageotropic to orthotropic",
                    "ageotropic to plagiotropic",
                    "axis 1 - orthotropic axis 2 - plagiotropic and orthotropic to orthotropic",
                    "axis 1 - orthotropic axis 2 - plagiotropic to orthotropic",
                    "axis 1 - orthotropic axis 2 - plagiotropic to plagiotropic",
                    "axis 1 - orthotropic axis 2 - plagiotropic to orthotropic",
                    "axis 1 - orthotropic axis 2 - plagiotropic to orthotropic",
                    "orthotropic to orthotropic",
                    "axis 1 - orthotropic axis 2 - plagiotropic to orthotropic",
                    "axis 1 - orthotropic axis 2 - plagiotropic to plagiotropic",
                    "nodata",
                    "orthotropic to ageotropic",
                    "orthotropic to orthotropic",
                    "orthotropic to ageotropic",
                    "orthotropic and plagiotropic to orthotropic",
                    "orthotropic to ageotropic",
                    "orthotropic to axis 1 - orthotropic axis 2 - plagiotropic",
                    "orthotropic to axis 1 - orthotropic axis 2 - plagiotropic and orthotropic",
                    "orthotropic to orthotropic",
                    "orthotropic to axis 1 orthotropic and plagiotropic",
                    "orthotropic to plagiotropic",
                    "plagiotropic to ageotropic",
                    "plagiotropic to axis 1 - orthotropic axis 2 - plagiotropic",
                    "plagiotropic to axis 1 - orthotropic axis 2 - plagiotropic",
                    "plagiotropic to orthotropic")
        )



conversion.BranchingMechanism <-
    data.frame(
        original = sort(unique(c(levels(as.factor(cultivation[cultivation$Trait == "BranchingMechanism", "Direction.WR2.to.Crop"])),
                            levels(as.factor(cultivation[cultivation$Trait == "BranchingMechanism", "Direction.WR2.to.CWR1"]))
                            ))),
        grouped = c("absent to continuous",
                    "absent to rhythmic",
                    "continuous to absent",
                    "continuous to diffuse",
                    "continuous to rhythmic",
                    "diffuse to continuous",
                    "diffuse to rhythmic",
                    "nodata",
                    "rhythmic to absent",
                    "rhythmic to continuous",
                    "rhythmic to diffuse")
        )

conversion.BranchingPosition <-
    data.frame(
        original = sort(unique(c(levels(as.factor(cultivation[cultivation$Trait == "BranchingPosition", "Direction.WR2.to.Crop"])),
                                 levels(as.factor(cultivation[cultivation$Trait == "BranchingPosition", "Direction.WR2.to.CWR1"]))
                                 ))),
        grouped = c("basitonic and mesotonic to absent",
                    "absent to basitonic",                               
                    "absent to mesotonic",                               
                    "absent to mesotonic and basitonic",
                    "acrotonic and basitonic to acrotonic",
                    "acrotonic and basitonic to basitonic",              
                    "acrotonic and basitonic to acrotonic and basitonic",
                    "acrotonic and basitonic to mesotonic",              
                    "acrotonic to absent",                               
                    "acrotonic to basitonic",                            
                    "acrotonic to acrotonic and basitonic",              
                    "acrotonic to basitonic and mesotonic",              
                    "acrotonic to mesotonic",                            
                    "acrotonic and basitonic to acrotonic",              
                    "acrotonic and basitonic to acrotonic and mesotonic",
                    "acrotonic and basitonic to basitonic",              
                    "acrotonic and basitonic to mesotonic",              
                    "basitonic and mesotonic",                           
                    "basitonic and mesotonic to acrotonic",
                    "basitonic and mesotonic to basitonic",
                    "basitonic and mesotonic to acrotonic and basitonic",
                    "basitonic and mesotonic to mesotonic",              
                    "basitonic to absent",                               
                    "basitonic to acrotonic",                            
                    "basitonic to acrotonic and basitonic",              
                    "basitonic to basitonic and mesotonic",              
                    "basitonic to mesotonic",                            
                    "basitonic to basitonic and mesotonic",              
                    "acrotonic and mesotonic to basitonic",              
                    "acrotonic and mesotonic to mesotonic",              
                    "acrotonic and mesotonic",                           
                    "acrotonic and mesotonic to acrotonic",              
                    "acrotonic and mesotonic to acrotonic and basitonic",
                    "acrotonic and mesotonic to basitonic",              
                    "acrotonic and mesotonic to mesotonic",              
                    "acrotonic and mesotonic to mesotonic",             
                    "acrotonic and mesotonic to acrotonic",             
                    "mesotonic to mesotonic and basitonic",                 
                    "mesotonic to absent",                               
                    "mesotonic to acrotonic",                            
                    "mesotonic to acrotonic and basitonic",              
                    "mesotonic to basitonic",                            
                    "mesotonic to acrotonic and basitonic",              
                    "mesotonic to basitonic and mesotonic",              
                    "mesotonic to basitonic and mesotonic",              
                    "mesotonic to mesotonic and basitonic", 
                    "nodata"                                            )
    )


conversion.BranchingType <-
    data.frame(
        original = sort(unique(c(levels(as.factor(cultivation[cultivation$Trait == "BranchingType", "Direction.WR2.to.Crop"])),
                                 levels(as.factor(cultivation[cultivation$Trait == "BranchingType", "Direction.WR2.to.CWR1"]))
                                 ))),
        grouped = c("absent to monopodial",
                    "absent to monopodial and sympodial",
                    "absent to monopodial",
                    "absent to sympodial",
                    "sympodial to axis 1 - monopodial axis 2 - sympodial",
                    "axis 1 - monopodial axis 2 - sympodial to monopodial",
                    "axis 1 monopodial and axis 2 plagiotrop apposition to monopodial",
                    "axis 1 monopodial and axis 2 plagiotrop apposition to sympodial",
                    "axis 1 - monopodial axis 2 - sympodial to monopodial",
                    "axis 1 - monopodial axis 2 - sympodial to monopodial",
                    "axis 1 - monopodial axis 2 - sympodial to sympodial",
                    "monopodial and sympodial to monopodial",
                    "monopodial and sympodial to sympodial",
                    "monopodial and sympodial to absent",
                    "monopodial and sympodial to sympodial",
                    "axis 1 - monopodial axis 2 - sympodial to monopodial",
                    "monopodial to sympodial",
                    "monopodial to absent",
                    "monopodial to monopodial",
                    "monopodial to monopodial and sympodial",
                    "monopodial to sympodial",
                    "monopodial to absent",
                    "nodata",
                    "sympodial to absent",
                    "sympodial to monopodial"
                    )
)

conversion.FloweringAxis <-
    data.frame(
        original = sort(unique(c(levels(as.factor(cultivation[cultivation$Trait == "FloweringAxis", "Direction.WR2.to.Crop"])),
                            levels(as.factor(cultivation[cultivation$Trait == "FloweringAxis", "Direction.WR2.to.CWR1"]))
                            ))),
        grouped = c("absent to terminal",
                    "lateral and terminal to lateral",
                    "lateral and terminal to terminal",
                    "lateral to lateral and terminal",
                    "lateral to terminal",
                    "lateral to lateral and terminal",
                    "nodata",
                    "lateral and terminal to lateral",
                    "lateral and terminal to terminal",
                    "terminal to absent",
                    "terminal to lateral",
                    "terminal to lateral and terminal",
                    "terminal to lateral or terminal",
                    "terminal to lateral and terminal",
                    "terminal to lateral or terminal")
    )


conversion.MeristemFunction <-
    data.frame(
        original = sort(unique(c(levels(as.factor(cultivation[cultivation$Trait == "MeristemFunction", "Direction.WR2.to.Crop"])),
                            levels(as.factor(cultivation[cultivation$Trait == "MeristemFunction", "Direction.WR2.to.CWR1"]))
                            ))),
        grouped = c("determinate to indeterminate",
                    "determinate to indeterminate",
                    "indeterminate to determinate",
                    "indeterminate to determinate",
                    "nodata")
    )


conversion.Phyllotaxis <-
    data.frame(
        original = sort(unique(c(levels(as.factor(cultivation[cultivation$Trait == "Phyllotaxis", "Direction.WR2.to.Crop"])),
                            levels(as.factor(cultivation[cultivation$Trait == "Phyllotaxis", "Direction.WR2.to.CWR1"]))
                            ))),
        grouped = c("axis 1 - spiral axis 2 - distichous to spiral",
                    "spiral to axis 1 - spiral axis 2 - distichous",
                    "distichous to opposite",
                    "distichous to spiral",
                    "distichous to spiral and distichous",
                    "nodata",
                    "opposite to distichous",
                    "opposite to spiral",
                    "opposite to opposite/spiral",
                    "opposite/spiral to spiral",
                    "spiral to axis 1 - spiral axis 2 - distichous",
                    "spiral and distichous to distichous",
                    "spiral and distichous to spiral",
                    "spiral to axis 1 - spiral axis 2 - distichous",
                    "spiral to distichous",
                    "spiral to opposite",
                    "spiral to opposite/spiral",
                    "spiral to spiral and distichous",
                    "spiral to opposite/spiral",
                    "opposite/spiral to spiral")
    )

fullconversion <- rbind(conversion.BranchingMechanism,
                        conversion.BranchingPosition,
                        conversion.BranchingType,
                        conversion.FloweringAxis,
                        conversion.Growthdirection,
                        conversion.MeristemFunction,
                        conversion.Phyllotaxis)
                        

cultivation.new <- cultivation[!(cultivation$Trait %in% c("Branching", "ShortShoots")),]
cultivation.new$Direction.WR2.to.Crop <- fullconversion$grouped[match(cultivation.new$Direction.WR2.to.Crop, fullconversion$original)]
cultivation.new$Direction.WR2.to.CWR1 <- fullconversion$grouped[match(cultivation.new$Direction.WR2.to.CWR1, fullconversion$original)]



## # Group corresponding phenotypes
## group_value_function <- function(x) {
##     ## Branching Mechanism
##     x <- gsub("rhytmic", "rhythmic", x)  # Correct typos
##     x <- gsub("rhty", "rhythmic", x)  # Correct typos

##     ## Branching Position
##     x <- gsub("acrotonic and basitoic", "basitonic and acrotonic", x, fixed = TRUE)
##     x <- gsub("acrotonic & basitonic", "basitonic and acrotonic", x, fixed = TRUE)
##     x <- gsub("acrotonic & basitonic", "basitonic and acrotonic", x, fixed = TRUE)
##     x <- gsub("basitonic & mesotonic", "mesotonic and basitonic", x, fixed = TRUE)
##     x <- gsub("basitonic & mesotonic", "mesotonic and basitonic", x, fixed = TRUE)
##     x <- gsub("acrotonic & mesotonic", "mesotonic and acrotonic", x, fixed = TRUE)
##     x <- gsub("acrotonic & mesotonic", "mesotonic and acrotonic", x, fixed = TRUE)

##     ## Branching Type
##     x <- gsub("axe 1 monopodial axis 2 sympodial", "monopodial and sympodial", x, fixed = TRUE)
##     x <- gsub("axis 1 monopodial axis 2 sympodial", "monopodial and sympodial", x, fixed = TRUE)
##     x <- gsub("axis 1 monopodial", "monopodial", x, fixed = TRUE)
##     x <- gsub("axis1 monopodial", "monopodial", x, fixed = TRUE)
##     x <- gsub("monopodial &  sympodial", "monopodial and sympodial", x, fixed = TRUE)
##     x <- gsub("monopodial & sympodial", "monopodial and sympodial", x, fixed = TRUE)
##     x <- gsub("monopodial  axis 1 and sympodial  axis 2", "monopodial and sympodial", x, fixed = TRUE)
##     x <- gsub("monopodial or sympodial", "monopodial and sympodial", x, fixed = TRUE)
##     ## x <- gsub("Monopodial +Axis 1", "Axis 1 - Monopodial", x, ignore.case = TRUE)
##     ## x <- gsub("Sympodial +Axis 2", "Axis 2 - Sympodial", x, ignore.case = TRUE)
##     x <- gsub("Monopodial +Axis 1", "monopodial", x, ignore.case = TRUE)
##     x <- gsub("Sympodial +Axis 2", "sympodial", x, ignore.case = TRUE)

##     ## Flowering Axis
##     x <- gsub("terminodatal", "terminal", x)  # Correct specific errors
##     x <- gsub("lareral", "lateral", x)  # Correct specific errors
##     x <- gsub("lateral & terminodatal", "lateral and terminal", x)
##     x <- gsub("lateral or terminodatal", "lateral or terminal", x)
##     x <- gsub("terminal  and lateral", "lateral and terminal", x)
##     x <- gsub("terminal and lateral", "lateral and terminal", x)
##     x <- gsub("terminal or lateral", "lateral or terminal", x)
##     x <- gsub("terminodatal and lateral", "lateral and terminal", x)
##     x <- gsub("terminodatal or lateral", "lateral or terminal", x)
##     x <- gsub("terminodatale", "terminal", x)
##     x <- gsub("terminale", "terminal", x)

##     ## Growth Direction
##     x <- gsub("axis 1 orthtropic axis 2 plagiotropic", "axis 1 - orthotropic axis 2 - plagiotropic", x, fixed = TRUE)
##     x <- gsub("axis1 orthtropic & axis 2 plagiotropic", "axis 1 - orthotropic axis 2 - plagiotropic", x, fixed = TRUE)
##     x <- gsub("axis 1 orthotropic", "orthotropic", x)
##     x <- gsub("axis 1- orthotropic", "orthotropic", x)
##     x <- gsub("orthortropic", "orthotropic", x)
##     x <- gsub("orthotropi ", "orthotropic", x, fixed = TRUE)
##     x <- gsub("orthtropic", "orthotropic", x)

##     ## Meristem Function
##     x <- gsub("axis 1 indeterminate", "indeterminate", x)
##     x <- gsub("indeterminodatate", "indeterminate", x)
##     x <- gsub("determinodatate", "determinate", x)
##     ## x <- gsub("determinodatate and indeterminodatate", "determinate", x)
##     ## x <- gsub("determinodatate or indeterminodatate", "determinate", x)
##     ## x <- gsub("indeterminodatate or determinodatate", "determinate", x)

##     ## Phyllotaxis
##     x <- gsub("axis 1 - spiral axis 2 - distichous", "axis 1 - spiral axis 2 - distichous", x, fixed = TRUE)
##     x <- gsub("axis 1 spiral and aes 2 distichous", "axis 1 - spiral axis 2 - distichous", x, fixed = TRUE)
##     x <- gsub("axis 1 spiral and axis 2 distichous", "axis 1 - spiral axis 2 - distichous", x, fixed = TRUE)
##     x <- gsub("opposite & spiral", "opposite/spiral", x, fixed = TRUE)
##     x <- gsub("opposite and spiral", "opposite/spiral", x, fixed = TRUE)
##     x <- gsub("spiral & distichous", "axis 1 - spiral axis 2 - distichous", x, fixed = TRUE)
##     x <- gsub("spiral and distichous", "axis 1 - spiral axis 2 - distichous", x, fixed = TRUE)
##     x <- gsub("spiral/opposite", "opposite/spiral", x, fixed = TRUE)
##     x <- gsub("\\bspira\\b(?!l)", "spiral", x, perl = TRUE)
##     x <- gsub("spipral", "spiral", x)  # Correct typos
##     return(x)
## }

# Apply the grouping function to multiple columns in the 'cultivation' dataframe
## cultivation <- cultivation %>%
##   mutate(VALUE = group_value_function(VALUE)) %>%
##   mutate(Direction.WR2.to.Crop = group_value_function(Direction.WR2.to.Crop)) %>%
##   mutate(Direction.WR2.to.CWR1 = group_value_function(Direction.WR2.to.CWR1))

# Optionally, print the updated dataframe to verify changes
#print(cultivation)

# Get the list of unique traits from the cultivation dataset
unique_traits <- unique(cultivation.new$Trait)

# Initialize a list to store the combined data for each trait
all_traits_data <- list()

# Loop through each trait
for (trait in unique_traits) {
  # Filter data for the current trait where conditions are met
  test <- cultivation.new %>%
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

write.csv(final_data, here::here("Data/TraitProportions.csv"), row.names = FALSE)

ggplot(final_data, aes(x = State, y = Proportion, fill = Before)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Trait, scales = "free_y") + # Faceting by Trait with free y scales if needed
  scale_fill_manual(values = rainbow(length(unique(final_data$Before)))) + # Setting colors
  labs(x = "State", y = "Proportion", fill = "Before") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Adjust text angle for clarity
