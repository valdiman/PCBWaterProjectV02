

# Install packages
install.packages("dplyr")
install.packages("tidyr")
install.packages("readr")

# Load libraries
library(dplyr)
library(tidyr)
library(readr)

# Read the CSV file into a data frame
PS_data <- read_csv("Data/PassaicRiver/2017-2019_OU2_PDI_Water_Column_20210924.csv")

# Arrange the data to ensure it's ordered correctly (if needed)
PS_data <- PS_data %>% arrange(SAMPLE_NAME, SAMPLE_DATE)

# Filter and keep only the columns you need
PS_data <- PS_data %>%
  select(SAMPLE_NAME, SAMPLE_DATE, SAMPLETIME, SAMPLE_TYPE_CODE, MATRIX_CODE,
         ANALYTIC_METHOD, LATITUDE, LONGITUDE, REPORT_RESULT_UNIT, CHEMICAL_NAME,
         RESULT_NUMERIC, TEST_TYPE)

# Filter rows based on the ANALYTIC_METHOD condition & just water samples (i.e., pg/l)
# and initial analysis
PS_data <- PS_data %>%
  filter(ANALYTIC_METHOD == "E1668A", REPORT_RESULT_UNIT == "pg/l",
         TEST_TYPE == "INITIAL", SAMPLE_TYPE_CODE == "N")

# Modify the CHEMICAL_NAME column to extract "PCB X" where X is the number
PS_data <- PS_data %>%
  mutate(CHEMICAL_NAME = sub(".*\\(PCB (\\d+)\\).*", "PCB \\1", CHEMICAL_NAME))

# Filter rows where CHEMICAL_NAME contains "PCB" and not "Total"
# and "Polychlorinated Biphenyl (PCB)"
PS_data <- PS_data %>%
  filter(grepl("PCB", CHEMICAL_NAME, fixed = TRUE) &
           !grepl("Total", CHEMICAL_NAME, fixed = TRUE) &
           !grepl("Polychlorinated Biphenyl", CHEMICAL_NAME))

# Create a new data frame with transposed values
transposed_data <- PS_data %>%
  pivot_wider(
    id_cols = c(SAMPLE_NAME, SAMPLE_DATE, SAMPLETIME, SAMPLE_TYPE_CODE,
                MATRIX_CODE, ANALYTIC_METHOD, LATITUDE, LONGITUDE,
                REPORT_RESULT_UNIT),
    names_from = CHEMICAL_NAME,
    values_from = RESULT_NUMERIC,
    values_fn = list  # Treat duplicate values as lists
  )

# Read the JSON file with new congener list from code NewPCBList.R
pcb_groups <- read_json("pcb_groups.json")

# Create an empty data frame to store the grouped data
grouped_data <- PS_data %>%
  distinct(SAMPLE_NAME, SAMPLE_DATE, SAMPLETIME, SAMPLE_TYPE_CODE, MATRIX_CODE,
           ANALYTIC_METHOD, LATITUDE, LONGITUDE, REPORT_RESULT_UNIT)

# Remove spaces and special characters from column names in transposed_data
colnames(transposed_data) <- gsub("[^[:alnum:]]", "", colnames(transposed_data))

# Function to safely coerce to numeric, replacing non-numeric values with NAs
safe_as_numeric <- function(x) {
  as.numeric(as.character(x))
}

# Iterate through the pcb_groups list and sum columns
for (group_name in names(pcb_groups)) {
  group_columns <- pcb_groups[[group_name]]
  group_columns <- gsub("[^[:alnum:]]", "", group_columns)  # Clean group column names
  
  # Coerce group_columns to numeric
  transposed_data[group_columns] <- lapply(transposed_data[group_columns], safe_as_numeric)
  
  # Calculate row sums, replacing NA with 0
  grouped_data[group_name] <- rowSums(transposed_data[group_columns], na.rm = TRUE)
  
  # Clean up: replace 0 with NA
  grouped_data[group_name][grouped_data[group_name] == 0] <- NA
}

# Remove SAMPLE_TYPE_CODE and MATRIX_CODE columns from grouped_data
grouped_data <- grouped_data %>%
  select(-SAMPLE_TYPE_CODE, -MATRIX_CODE, -ANALYTIC_METHOD)

# Create a new column named "tPCB" that sums columns 6 to 109
grouped_data <- grouped_data %>%
  mutate(tPCB = rowSums(select(., starts_with("PCB")), na.rm = TRUE))

# Replace 0s with NA for columns 6 to 109
grouped_data <- grouped_data %>%
  mutate_at(vars(starts_with("PCB")), ~ ifelse(. == 0, NA, .))

# Change SYS_SAMPLE_CODE
grouped_data <- grouped_data %>%
  mutate(SAMPLE_NAME = substr(SAMPLE_NAME, 1, nchar(SAMPLE_NAME) - 3))

# Remove SAMPLE_NAME with no name (NA)
grouped_data <- grouped_data %>%
  filter(!is.na(SAMPLE_NAME))

# Change the name of columns to be consistent
colnames(grouped_data)[colnames(grouped_data) == "REPORT_RESULT_UNIT"] <- "RESULT_UNIT"

# Export results
write.csv(grouped_data,
          file = "Data/PassaicRiver/2017-2019_OU2_PDI_Water_Column_20210924V0.csv")

# Remove SAMPLETIME columns from grouped_data & create a new data.frame
grouped_dataFV <- grouped_data %>%
  select(-SAMPLETIME)

# Export results
write.csv(grouped_dataFV,
          file = "Data/PassaicRiver/2017-2019_OU2_PDI_Water_Column_20210924FV.csv")

# in CombinePassaicData: Pass.9
