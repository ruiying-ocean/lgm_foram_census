#---------------------------------------------------
# merge all relative abundance (LGM Species)
#---------------------------------------------------
process_file <- function(file, ecolvl, abundance_col, file_pattern) {
  print(file)
  df <- read_csv(file)

  # If contains Core column, rename it to "Event"
  if ("Core" %in% colnames(df)) {
    df <- df %>% rename(Event = Core)
  }

  # Candidate columns for Depth
  candidate_cols <- c(
    "Depth sed [m]", "Depth top [m]", "depth",
    "Depth sed", "Sample depth - upper (m)"
  )
  common_cols <- intersect(candidate_cols, colnames(df))

  # If contains any of the candidate Depth columns, rename it to "Depth [m]"
  if (length(common_cols) > 0) {
    df <- df %>% rename(`Depth [m]` = common_cols[1])
  }

  # Select specific columns
  df <- df %>% select(c("Event", "Latitude", "Longitude", "Depth [m]", all_of(ecolvl), all_of(abundance_col)))

  # Extract additional information from the file name
  df <- df %>% mutate(Data_Source = str_extract(basename(file), file_pattern), .after = Event)

  return(df)
}

filenames <- list.files("sp", pattern = "^lgm.*_r\\.csv$", full.names = TRUE)
ldf <- lapply(filenames, process_file,
  ecolvl = "Species",
  abundance_col = "Relative Abundance",
  file_pattern = "(?<=lgm_).*?(?=_sp_r\\.csv)"
)

wdf <- ldf %>%
  rbindlist() %>%
  pivot_wider(
    id_cols = c("Event", "Data_Source", "Depth [m]", "Latitude", "Longitude"),
    values_from = "Relative Abundance",
    names_from = "Species",
    values_fill = NA,
    values_fn = list
  )

wdf <- wdf %>% unnest(where(is.list))
## remove duplicate rows
wdf <- wdf[!duplicated(wdf), ]

wdf %>% fwrite("tidy/lgm_sp_r_tidy.csv")

#---------------------------------------------------
# merge all absolute abundance (LGM Species)
#---------------------------------------------------

filenames <- list.files("sp", pattern = "^lgm.*_a\\.csv$", full.names = TRUE)
ldf <- lapply(filenames, process_file,
  ecolvl = "Species",
  abundance_col = "Absolute Abundance",
  file_pattern = "(?<=lgm_).*?(?=_sp_a\\.csv)"
)

wdf <- ldf %>%
  rbindlist() %>%
  pivot_wider(
    id_cols = c("Event", "Data_Source", "Depth [m]", "Latitude", "Longitude"),
    values_from = "Absolute Abundance",
    names_from = "Species",
    values_fill = NA,
    values_fn = list
  )

wdf <- wdf %>% unnest(where(is.list))

## remove duplicate rows
wdf <- wdf[!duplicated(wdf), ]

wdf %>% fwrite("tidy/lgm_sp_a_tidy.csv")

#---------------------------------------------------
# merge all absolute abundance (LGM Functional Group)
#---------------------------------------------------

filenames <- list.files("fg", pattern = "^lgm.*_a\\.csv$", full.names = TRUE)

ldf <- ldf <- lapply(filenames, process_file,
  ecolvl = c("Symbiosis", "Spine"),
  abundance_col = "Absolute Abundance",
  file_pattern = "(?<=lgm_).*?(?=_fg_a\\.csv)"
)

ldf <- ldf %>% rbindlist()

ldf <- ldf %>% mutate(`Functional Group` = paste(Symbiosis, Spine, sep = " "))

wdf <- ldf %>%
  pivot_wider(
    id_cols = c("Event", "Data_Source", "Depth [m]", "Latitude", "Longitude"),
    values_from = "Absolute Abundance",
    names_from = "Functional Group",
    values_fill = NA,
    values_fn = list
  )

wdf <- wdf %>% unnest(where(is.list))

## remove duplicate rows
wdf <- wdf[!duplicated(wdf), ]

wdf %>% fwrite("tidy/lgm_fg_a_tidy.csv")

#---------------------------------------------------
# merge all relative abundance (LGM Functional Group)
#---------------------------------------------------

filenames <- list.files("fg", pattern = "^lgm.*_r\\.csv$", full.names = TRUE)
ldf <- lapply(filenames, process_file,
  ecolvl = c("Symbiosis", "Spine"),
  abundance_col = "Relative Abundance",
  file_pattern = "(?<=lgm_).*?(?=_fg_r\\.csv)"
)

ldf <- ldf %>%
  rbindlist() %>%
  mutate(`Functional Group` = paste(Symbiosis, Spine, sep = " "))

wdf <- ldf %>% pivot_wider(
  id_cols = c("Event", "Data_Source", "Depth [m]", "Latitude", "Longitude"),
  values_from = "Relative Abundance",
  names_from = "Functional Group",
  values_fill = NA,
  values_fn = list
)

wdf <- wdf %>% unnest(where(is.list))

## remove duplicate rows
wdf <- wdf[!duplicated(wdf), ]

wdf %>% fwrite("tidy/lgm_fg_r_tidy.csv")

#--------------------
# Modern fg relative/absolute abundance
#--------------------
pi_a <- read_csv("fg/forcens_fg_a.csv") %>% mutate(`Functional Group` = paste(Symbiosis, Spine, sep = " "))

pi_aw <- pi_a %>%
  pivot_wider(
    id_cols = c("Latitude", "Longitude"),
    values_from = "Absolute Abundance",
    names_from = "Functional Group",
    values_fill = NA,
    values_fn = list
  ) %>%
  distinct() %>%
  unnest(cols = `symbiont-barren spinose`:`symbiont-facultative spinose`)

pi_aw %>% fwrite("tidy/forcens_fg_a_tidy.csv")

pi_r <- read_csv("fg/forcens_fg_r.csv") %>% mutate(`Functional Group` = paste(Symbiosis, Spine, sep = " "))

pi_rw <- pi_r %>%
  pivot_wider(
    id_cols = c("Latitude", "Longitude"),
    values_from = "Relative Abundance",
    names_from = "Functional Group",
    values_fill = NA,
    values_fn = list
  ) %>%
  distinct() %>%
  unnest(cols = `symbiont-barren spinose`:`symbiont-facultative spinose`)

pi_rw %>% fwrite("tidy/forcens_fg_r_tidy.csv")

## forcens species data already exported in clean_forcens.R
## so ignore here
