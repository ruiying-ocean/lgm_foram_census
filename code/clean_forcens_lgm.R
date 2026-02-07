## ForCenS_LGM dataset (replaces the 5 regional MARGO files)
## Produces the same output format as clean_margo.R

data <- fread("raw/LGM_foraminifera_assemblages_20240110.csv")

## -----------------------------------
## Remove composite / morphotype columns and identifiers
## -----------------------------------

cols_to_remove <- c(
  "Globorotalia_menardii_Globorotalia_tumida",
  "Globigerinoides_ruber_Globigerinoides_white",
  "Turborotalita_humilis_Berggrenia_pumilio",
  grep("_dextral_coiling$|_sinistral_coiling$|_w_sac_chamber$|_wo_sac_chamber$",
    names(data),
    value = TRUE
  ),
  "unidentified", "siteID", "sampleID"
)

data <- data %>% select(-any_of(cols_to_remove))

## -----------------------------------
## Rename metadata columns
## -----------------------------------

data <- data %>% rename(
  `Coring device` = CoringDevice,
  `Water depth (m)` = WaterDepth_m,
  `Sample depth - upper (m)` = SampleDepthTop_m
)

## -----------------------------------
## Clean species names (underscore -> space, then standardize taxonomy)
## -----------------------------------

## Identify species columns (everything after QCNote up to isCounts)
meta_cols <- c(
  "Core", "Coring device", "Latitude", "Longitude", "Water depth (m)",
  "Ocean", "SampleLabel", "Sample depth - upper (m)",
  "SampleDepthBottom_m", "SampleDepthMid_m",
  "CalendarAge_cal_ky_BP", "ChronozoneLevel", "SedimentationRate_cm_ky",
  "Laboratory", "Publication", "ChronologySource", "AddedBy",
  "DateOfAddition", "Source", "TaxonomicNotes", "Count", "QCNote",
  "isCounts", "inMARGO"
)
sp_cols <- setdiff(names(data), meta_cols)

## Replace underscores with spaces in species column names
new_sp_names <- gsub("_", " ", sp_cols)
setnames(data, sp_cols, new_sp_names)

## Standardize taxonomy using revise_sp_name
data <- data %>%
  revise_sp_name("Globorotalia menardii", "Globorotalia cultrata") %>%
  revise_sp_name("Globorotalia theyeri", "Globorotalia eastropacia") %>%
  revise_sp_name("Tenuitella iota", "Tenuitellita iota") %>%
  revise_sp_name("Globigerinoides ruber", "Globigerinoides ruber ruber") %>%
  revise_sp_name("Globigerinoides white", "Globigerinoides ruber albus") %>%
  revise_sp_name("Globoconella inflata", "Globorotalia inflata")

## Filter out rows with missing coordinates
data <- data %>% dplyr::filter(!is.na(Longitude) & !is.na(Latitude))

## Update species column list after renaming
sp_cols <- setdiff(names(data), meta_cols)

## -----------------------------------
## Helper: remove columns utility
## -----------------------------------

remove_columns <- function(data, column_names) {
  existing_columns <- intersect(column_names, names(data))
  if (length(existing_columns) > 0) {
    data <- data %>% select(-any_of(existing_columns))
  }
  return(data)
}

## -----------------------------------
## Relative abundance
## -----------------------------------

forcens_lgm_relative_abundance <- function(data) {
  perc_data <- data %>% dplyr::filter(isCounts == FALSE)
  count_data <- data %>% dplyr::filter(isCounts == TRUE)

  ## Convert count rows to percentages
  count_data <- count_data %>% mutate(across(
    all_of(sp_cols),
    ~ . / Count * 100
  ))

  ## Combine
  new_data <- rbind(perc_data, count_data)

  ## Remove isCounts and Count columns
  new_data <- new_data %>% remove_columns(c("isCounts", "Count", "inMARGO"))

  ## Convert percentages to decimals
  new_data <- new_data %>% mutate(across(
    all_of(sp_cols),
    ~ . / 100
  ))

  return(new_data)
}

## -----------------------------------
## Absolute abundance
## -----------------------------------

forcens_lgm_abs_abundance <- function(data) {
  count_data <- data %>% dplyr::filter(isCounts == TRUE)
  perc_data <- data %>% dplyr::filter(isCounts == FALSE)

  ## Convert percentage rows to counts (need Count column)
  perc_data <- perc_data %>%
    drop_na(Count) %>%
    mutate(across(
      all_of(sp_cols),
      ~ ceiling(. * Count / 100)
    ))

  ## Combine
  new_data <- rbind(count_data, perc_data)

  ## Remove isCounts and Count columns
  new_data <- new_data %>% remove_columns(c("isCounts", "Count", "inMARGO"))

  return(new_data)
}

## -----------------------------------
## Convert to long format and abbreviate species
## -----------------------------------

forcens_lgm_to_long <- function(data) {
  data <- data %>% remove_columns(c(
    "SampleLabel", "SampleDepthBottom_m", "SampleDepthMid_m",
    "CalendarAge_cal_ky_BP", "ChronozoneLevel", "SedimentationRate_cm_ky",
    "Laboratory", "Publication", "ChronologySource", "AddedBy",
    "DateOfAddition", "Source", "TaxonomicNotes", "QCNote"
  ))

  pivot_longer(data,
    cols = -c(
      Core, `Coring device`, Latitude, Longitude,
      `Water depth (m)`, Ocean, `Sample depth - upper (m)`
    ),
    names_to = "Species", values_to = "Abundance"
  )
}

## -----------------------------------
## Build outputs
## -----------------------------------

data_r <- forcens_lgm_relative_abundance(data)
data_a <- forcens_lgm_abs_abundance(data)

lgm_r <- forcens_lgm_to_long(data_r) %>%
  rename("Relative Abundance" = "Abundance") %>%
  rowwise() %>%
  mutate_at(.vars = "Species", .funs = species_abbrev)

lgm_a <- forcens_lgm_to_long(data_a) %>%
  rename("Absolute Abundance" = "Abundance") %>%
  rowwise() %>%
  mutate_at(.vars = "Species", .funs = species_abbrev)

fwrite(lgm_r, "sp/lgm_forcens_lgm_sp_r.csv")
fwrite(lgm_a, "sp/lgm_forcens_lgm_sp_a.csv")

## -----------------------------------
## Functional groups
## -----------------------------------

forcens_lgm_group_and_aggregate <- function(data) {
  symbiosis_short_tbl <- symbiosis_tbl %>%
    select(!c(Species)) %>%
    distinct()
  data_merged <- merge(data, symbiosis_short_tbl, by.x = "Species", by.y = "short_name") %>% select(!"Species")

  data_merged <- data_merged %>% remove_columns(c(
    "SampleLabel", "SampleDepthBottom_m", "SampleDepthMid_m",
    "CalendarAge_cal_ky_BP", "ChronozoneLevel", "SedimentationRate_cm_ky",
    "Laboratory", "Publication", "ChronologySource", "AddedBy",
    "DateOfAddition", "Source", "TaxonomicNotes", "QCNote"
  ))

  data_merged <- data_merged %>%
    group_by(
      Core, `Coring device`, Latitude, Longitude,
      `Water depth (m)`, Ocean, `Sample depth - upper (m)`,
      Symbiosis, Spine
    ) %>%
    summarise_all(.funs = sum, na.rm = T) %>%
    ungroup()

  return(data_merged)
}

lgm_r %>%
  forcens_lgm_group_and_aggregate() %>%
  dplyr::filter(`Relative Abundance` < 1.1) %>%
  fwrite("fg/lgm_forcens_lgm_fg_r.csv")

lgm_a %>%
  forcens_lgm_group_and_aggregate() %>%
  fwrite("fg/lgm_forcens_lgm_fg_a.csv")
