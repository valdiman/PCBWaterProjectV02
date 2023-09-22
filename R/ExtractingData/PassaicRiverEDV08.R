

# Install packages
install.packages("dplyr")
install.packages("tidyr")
install.packages("readr")

# Load libraries
library(dplyr)
library(tidyr)
library(readr)

# Read the CSV file into a data frame
PS_data <- read_csv("Data/PassaicRiver/2010-2016_Diamond-Alkali_OU3_EPA-NBSA-Split-samples.csv")

# Arrange the data to ensure it's ordered correctly (if needed)
PS_data <- PS_data %>% arrange(SYS_LOC_CODE, SAMPLE_DATE)

# Filter and keep only the columns you need
PS_data <- PS_data %>%
  select(SYS_LOC_CODE, SAMPLE_DATE, SAMPLE_TYPE_CODE, MATRIX_CODE, ANALYTIC_METHOD,
         LATITUDE, LONGITUDE, RESULT_UNIT, CHEMICAL_NAME, RESULT_NUMERIC)

# Filter rows based on the ANALYTIC_METHOD condition & just water samples
PS_data <- PS_data %>%
  filter(ANALYTIC_METHOD == "E1668 AXYS MLA010", RESULT_UNIT == "pg/l")


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
    id_cols = c(SYS_LOC_CODE, SAMPLE_DATE, SAMPLE_TYPE_CODE, MATRIX_CODE,
                ANALYTIC_METHOD, LATITUDE, LONGITUDE, RESULT_UNIT),
    names_from = CHEMICAL_NAME,
    values_from = RESULT_NUMERIC,
    values_fn = list  # Treat duplicate values as lists
  )

# Define new PCB list
pcb_groups <- list(
  "PCB1" = c("PCB1"), "PCB2" = c("PCB2"), "PCB3" = c("PCB3"),
  "PCB4+10" = c("PCB4", "PCB10"), "PCB5+8" = c("PCB5", "PCB8"),
  "PCB6" = c("PCB6"), "PCB7+9" = c("PCB7", "PCB9"),
  "PCB11" = c("PCB11"), "PCB12+13" = c("PCB12", "PCB13"),
  "PCB14" = c("PCB14"), "PCB15" = c("PCB15"),
  "PCB16+32" = c("PCB16", "PCB32"), "PCB17" = c("PCB17"),
  "PCB18+30" = c("PCB18", "PCB30"), "PCB19" = c("PCB19"),
  "PCB20+21+28+31+33+50+53" = c("PCB20", "PCB21", "PCB28", "PCB31", "PCB33",
                                "PCB50", "PCB53"), "PCB22" = c("PCB22"),
  "PCB23" = c("PCB23"), "PCB24+27" = c("PCB24", "PCB27"),
  "PCB25" = c("PCB25"), "PCB26+29" = c("PCB26", "PCB29"), "PCB34" = c("PCB34"),
  "PCB35" = c("PCB35"), "PCB36" = c("PCB36"), "PCB37+42" = c("PCB37", "PCB42"),
  "PCB38" = c("PCB38"), "PCB39" = c("PCB39"),
  "PCB40+41+64+71+72" = c("PCB40", "PCB41", "PCB64", "PCB71", "PCB72"),
  "PCB43+49+52+69+73" = c("PCB43", "PCB49", "PCB52", "PCB69", "PCB73"),
  "PCB44+47+65" = c("PCB44", "PCB47", "PCB65"), "PCB45+51" = c("PCB45", "PCB51"),
  "PCB46" = c("PCB46"), "PCB48+59+62+75" = c("PCB48", "PCB59", "PCB62", "PCB75"),
  "PCB54" = c("PCB54"), "PCB55" = c("PCB55"), "PCB56+60" = c("PCB56", "PCB60"),
  "PCB57" = c("PCB57"), "PCB58" = c("PCB58"),
  "PCB61+66+70+74+76+93+95+98+100+102" = c("PCB61", "PCB66", "PCB70", "PCB74",
                                           "PCB76", "PCB93", "PCB95", "PCB98",
                                           "PCB100", "PCB102"), "PCB63" = c("PCB63"),
  "PCB67" = c("PCB67"), "PCB68" = c("PCB68"),
  "PCB77+85+110+111+115+116+117" = c("PCB77", "PCB85", "PCB110", "PCB111",
                                     "PCB115", "PCB116", "PCB117"),
  "PCB78" = c("PCB78"), "PCB79" = c("PCB79"), "PCB80" = c("PCB80"),
  "PCB81+86+87+97+107+108+109+112+119+124+125" = c("PCB81", "PCB86", "PCB87",
                                                   "PCB97", "PCB107", "PCB108",
                                                   "PCB109", "PCB112", "PCB119",
                                                   "PCB124", "PCB125"),
  "PCB82+135+144+151+154" = c("PCB82", "PCB135", "PCB144", "PCB151", "PCB154"),
  "PCB83+99" = c("PCB83", "PCB99"), "PCB84+92" = c("PCB84", "PCB92"),
  "PCB88+91" = c("PCB88", "PCB91"), "PCB89" = c("PCB89"),
  "PCB90+101+113" = c("PCB90", "PCB101", "PCB113"),
  "PCB94" = c("PCB94"), "PCB96" = c("PCB96"), "PCB103" = c("PCB103"),
  "PCB104" = c("PCB104"), "PCB105" = c("PCB105"),
  "PCB106+118" = c("PCB106", "PCB118"),
  "PCB114+122+131+133+142+146+165" = c("PCB114", "PCB122", "PCB131", "PCB133",
                                       "PCB142", "PCB146", "PCB165"),
  "PCB120" = c("PCB120"), "PCB121" = c("PCB121"),
  "PCB123+139+140+147+149" = c("PCB123", "PCB139", "PCB140", "PCB147", "PCB149"),
  "PCB126" = c("PCB126"), "PCB127" = c("PCB127"),
  "PCB128+162+166+167" = c("PCB128", "PCB162", "PCB166", "PCB167"),
  "PCB129+137+138+158+160+163+164+176+178" = c("PCB129", "PCB137", "PCB138",
                                               "PCB158", "PCB160", "PCB163",
                                               "PCB164", "PCB176", "PCB178"),
  "PCB130" = c("PCB130"), "PCB132+153+161+168" = c("PCB132", "PCB153",
                                                   "PCB161", "PCB168"),
  "PCB134+143" = c("PCB134", "PCB143"), "PCB136" = c("PCB136"),
  "PCB141" = c("PCB141"), "PCB145" = c("PCB145"), "PCB148" = c("PCB148"),
  "PCB150" = c("PCB150"), "PCB152" = c("PCB152"), "PCB155" = c("PCB155"),
  "PCB156+157+172+197+200" = c("PCB156", "PCB157", "PCB172", "PCB197", "PCB200"),
  "PCB159" = c("PCB159"), "PCB169" = c("PCB169"),
  "PCB170+190" = c("PCB170", "PCB190"), "PCB171+173+202" = c("PCB171", "PCB173",
                                                             "PCB202"),
  "PCB174" = c("PCB174"), "PCB175" = c("PCB175"), "PCB177" = c("PCB177"),
  "PCB179" = c("PCB179"), "PCB180+193" = c("PCB180", "PCB193"),
  "PCB181" = c("PCB181"), "PCB182+187" = c("PCB182", "PCB187"),
  "PCB183+185" = c("PCB183", "PCB185"), "PCB184" = c("PCB184"),
  "PCB186" = c("PCB186"), "PCB188" = c("PCB188"), "PCB189" = c("PCB189"),
  "PCB191" = c("PCB191"), "PCB192" = c("PCB192"),
  "PCB194+205" = c("PCB194", "PCB205"), "PCB195+208" = c("PCB195", "PCB208"),
  "PCB196+203" = c("PCB196", "PCB203"), "PCB198+199+201" = c("PCB198", "PCB199",
                                                             "PCB201"),
  "PCB204" = c("PCB204"), "PCB206" = c("PCB206"), "PCB207" = c("PCB207"),
  "PCB209" = c("PCB209")
)

# Create an empty data frame to store the grouped data
grouped_data <- PS_data %>%
  distinct(SYS_LOC_CODE, SAMPLE_DATE, SAMPLE_TYPE_CODE, MATRIX_CODE,
           ANALYTIC_METHOD, LATITUDE, LONGITUDE, RESULT_UNIT)

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

# Remove SAMPLE_NAME with no name (NA)
grouped_data <- grouped_data %>%
  filter(!is.na(SYS_LOC_CODE))

# Export results
write.csv(grouped_data, file = "Data/PassaicRiver/2010-2016_Diamond-Alkali_OU3_EPA-NBSA-Split-samplesFV.csv")

