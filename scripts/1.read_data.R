# 1. read_data.R ---------------------------------------------------------------

library(readxl)


# ------------------------------------------------------------------------------
# 0. File and dataset settings
# ------------------------------------------------------------------------------

data_file <- "data (raw).xlsx"

expected_raw_columns <- 49L
expected_final_columns <- 46L
expected_rows <- 866L

expected_start_date <- as.Date("2020-01-30")
expected_end_date <- as.Date("2022-06-13")


# ------------------------------------------------------------------------------
# 1. Check and import the Excel file
# ------------------------------------------------------------------------------

if (!file.exists(data_file)) {
  stop(
    "Cannot find the input file: ",
    data_file
  )
}

df <- readxl::read_excel(
  path = data_file,
  na = c(
    "",
    "NA",
    "N/A",
    "NULL",
    "null"
  ),
  trim_ws = TRUE
)

df <- as.data.frame(
  df,
  check.names = FALSE
)

if (ncol(df) != expected_raw_columns) {
  stop(
    "Expected ",
    expected_raw_columns,
    " columns in the raw Excel file, but found ",
    ncol(df),
    ". Check whether the workbook has changed."
  )
}

if (nrow(df) != expected_rows) {
  stop(
    "Expected ",
    expected_rows,
    " observations in the raw Excel file, but found ",
    nrow(df),
    ". Check whether the workbook has changed."
  )
}


# ------------------------------------------------------------------------------
# 2. Assign consistent analysis variable names
# ------------------------------------------------------------------------------

original_column_names <- names(df)

analysis_column_names <- c(
  "Date",
  "Excess.mortality",
  "Schools",
  "Face.masks",
  "Lockdown.severity",
  "Majority.COVID.19.variant",
  "Flights.7.day.moving.average",
  "OpenTable.restaurant.bookings.London.index",
  "Google.homeworking.Greater.London.mobility.index",
  "Google.workplace.Greater.London.mobility.index",
  "Apple.walking.London.mobility.index",
  "Google.parks.Greater.London.mobility.index",
  "Google.retail.recreation.Greater.London.mobility.index",
  "Google.grocery.pharmacy.Greater.London.mobility.index",
  "Google.transit.stations.mobility.index",
  "TfL.Tube.mobility.index",
  "TfL.Bus.mobility.index",
  "Citymapper.journeys.mobility.index",
  "Season",
  "PCR.tests",
  "PCR.tests.capacity",
  "Antibody.tests",
  "Antibody.tests.capacity",
  "Pillar.1.NHS.and.UKHSA.capacity",
  "Pillar.2.UK.Government.capacity",
  "Pillar.3.Antibody.capacity",
  "Pillar.4.Surveillance.capacity",
  "Pillar.1.NHS.and.UKHSA.tests",
  "Pillar.2.UK.Government.tests",
  "Pillar.3.Antibody.tests",
  "Pillar.4.Surveillance.tests",
  "Tests.across.all.4.Pillars",
  "New.cases",
  "New.infections",
  "Reinfections",
  "New.cases.specimen.date.7.day.change.",
  "New.cases.specimen.date.7.day...change.",
  "Hospital.admissions",
  "Patients.in.hospital",
  "Patients.in.MVBs",
  "Vaccinations.total",
  "Vaccinations.1st.dose",
  "Vaccinations.2nd.dose",
  "Vaccinations.3rd.dose",
  "First.dose.uptake",
  "Second.dose.uptake",
  "Third.dose.uptake",
  "Deaths.within.28.days.of.COVID.19",
  "COVID-19.deaths.on.certificate"
)

if (length(analysis_column_names) != ncol(df)) {
  stop(
    "The number of analysis variable names does not ",
    "match the number of columns in the Excel file."
  )
}

analysis_column_names <- make.names(
  analysis_column_names,
  unique = TRUE
)

names(df) <- analysis_column_names

# Save a record linking the original Excel headings to the analysis names.
column_name_mapping <- data.frame(
  original_name = original_column_names,
  analysis_name = analysis_column_names,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 3. Remove variables that are not used in the analysis
# ------------------------------------------------------------------------------

columns_to_remove <- c(
  "New.cases.specimen.date.7.day.change.",
  "New.cases.specimen.date.7.day...change.",
  "Deaths.within.28.days.of.COVID.19"
)

missing_removal_columns <- setdiff(
  columns_to_remove,
  names(df)
)

if (length(missing_removal_columns) > 0L) {
  stop(
    "The following columns intended for removal were not found: ",
    paste(
      missing_removal_columns,
      collapse = ", "
    )
  )
}

df <- df[
  ,
  !names(df) %in% columns_to_remove,
  drop = FALSE
]

if (ncol(df) != expected_final_columns) {
  stop(
    "Expected ",
    expected_final_columns,
    " columns after removing unused variables, but found ",
    ncol(df),
    "."
  )
}


# ------------------------------------------------------------------------------
# 4. Parse and validate the date variable
# ------------------------------------------------------------------------------

parse_analysis_date <- function(x) {
  
  if (inherits(x, "Date")) {
    return(
      as.Date(x)
    )
  }
  
  if (inherits(x, "POSIXt")) {
    return(
      as.Date(x)
    )
  }
  
  if (is.numeric(x)) {
    return(
      as.Date(
        x,
        origin = "1899-12-30"
      )
    )
  }
  
  date_text <- trimws(
    as.character(x)
  )
  
  parsed_date <- as.Date(
    date_text,
    format = "%Y-%m-%d"
  )
  
  unresolved <- is.na(parsed_date) &
    !is.na(date_text)
  
  if (any(unresolved)) {
    
    parsed_date[unresolved] <- as.Date(
      date_text[unresolved],
      format = "%d/%m/%Y"
    )
  }
  
  unresolved <- is.na(parsed_date) &
    !is.na(date_text)
  
  if (any(unresolved)) {
    
    parsed_date[unresolved] <- as.Date(
      date_text[unresolved],
      format = "%d-%m-%Y"
    )
  }
  
  return(parsed_date)
}


df$Date <- parse_analysis_date(
  df$Date
)

if (anyNA(df$Date)) {
  
  invalid_date_rows <- which(
    is.na(df$Date)
  )
  
  stop(
    "Some dates could not be parsed. Affected row(s): ",
    paste(
      invalid_date_rows,
      collapse = ", "
    )
  )
}

df <- df[
  order(df$Date),
  ,
  drop = FALSE
]

rownames(df) <- NULL


# ------------------------------------------------------------------------------
# 5. Check the date range, ordering and daily continuity
# ------------------------------------------------------------------------------

if (anyDuplicated(df$Date)) {
  
  duplicated_dates <- unique(
    df$Date[
      duplicated(df$Date)
    ]
  )
  
  stop(
    "Duplicated dates were found: ",
    paste(
      duplicated_dates,
      collapse = ", "
    )
  )
}

date_differences <- as.integer(
  diff(df$Date)
)

if (any(date_differences != 1L)) {
  
  gap_positions <- which(
    date_differences != 1L
  )
  
  first_gap <- gap_positions[1]
  
  stop(
    "The observations are not consecutive daily observations. ",
    "The first problem occurs between ",
    as.character(df$Date[first_gap]),
    " and ",
    as.character(df$Date[first_gap + 1L]),
    "."
  )
}

if (min(df$Date) != expected_start_date ||
    max(df$Date) != expected_end_date) {
  
  stop(
    "Expected the date range ",
    as.character(expected_start_date),
    " to ",
    as.character(expected_end_date),
    ", but found ",
    as.character(min(df$Date)),
    " to ",
    as.character(max(df$Date)),
    "."
  )
}


# ------------------------------------------------------------------------------
# 6. Define the categorical and numeric variables
# ------------------------------------------------------------------------------

categorical_columns <- c(
  "Schools",
  "Face.masks",
  "Lockdown.severity",
  "Majority.COVID.19.variant",
  "Season"
)

missing_categorical_columns <- setdiff(
  categorical_columns,
  names(df)
)

if (length(missing_categorical_columns) > 0L) {
  stop(
    "The following expected categorical variables were not found: ",
    paste(
      missing_categorical_columns,
      collapse = ", "
    )
  )
}

numeric_columns <- setdiff(
  names(df),
  c(
    "Date",
    categorical_columns
  )
)

if (length(numeric_columns) != 40L) {
  stop(
    "Expected 40 numeric variables after preprocessing, but found ",
    length(numeric_columns),
    "."
  )
}


# ------------------------------------------------------------------------------
# 7. Convert and validate all numeric variables
# ------------------------------------------------------------------------------

convert_numeric_column <- function(
    x,
    column_name
) {
  
  if (is.numeric(x)) {
    
    converted_values <- as.numeric(x)
    
  } else {
    
    cleaned_values <- trimws(
      as.character(x)
    )
    
    cleaned_values[
      cleaned_values %in% c(
        "",
        "NA",
        "N/A",
        "NULL",
        "null"
      )
    ] <- NA_character_
    
    # Remove thousands separators if any values were stored as text.
    cleaned_values <- gsub(
      ",",
      "",
      cleaned_values,
      fixed = TRUE
    )
    
    converted_values <- suppressWarnings(
      as.numeric(cleaned_values)
    )
    
    invalid_values <- !is.na(cleaned_values) &
      is.na(converted_values)
    
    if (any(invalid_values)) {
      
      invalid_examples <- unique(
        cleaned_values[
          invalid_values
        ]
      )
      
      stop(
        "Non-numeric values were found in ",
        column_name,
        ": ",
        paste(
          head(invalid_examples, 10L),
          collapse = ", "
        )
      )
    }
  }
  
  nonfinite_values <- !is.na(converted_values) &
    !is.finite(converted_values)
  
  if (any(nonfinite_values)) {
    stop(
      "Non-finite values were found in ",
      column_name,
      "."
    )
  }
  
  return(converted_values)
}


for (column_name in numeric_columns) {
  
  df[[column_name]] <- convert_numeric_column(
    df[[column_name]],
    column_name = column_name
  )
}


# ------------------------------------------------------------------------------
# 8. Convert the five categorical variables to factors
# ------------------------------------------------------------------------------

convert_factor_column <- function(x) {
  
  factor_values <- trimws(
    as.character(x)
  )
  
  factor_values[
    factor_values %in% c(
      "",
      "NA",
      "N/A",
      "NULL",
      "null"
    )
  ] <- NA_character_
  
  factor(
    factor_values
  )
}


df[categorical_columns] <- lapply(
  df[categorical_columns],
  convert_factor_column
)

df <- droplevels(df)


# ------------------------------------------------------------------------------
# 9. Final dataset validation
# ------------------------------------------------------------------------------

if (nrow(df) != expected_rows) {
  stop(
    "Expected ",
    expected_rows,
    " final observations, but found ",
    nrow(df),
    "."
  )
}

if (ncol(df) != expected_final_columns) {
  stop(
    "Expected ",
    expected_final_columns,
    " final variables, but found ",
    ncol(df),
    "."
  )
}

numeric_type_check <- vapply(
  df[numeric_columns],
  is.numeric,
  logical(1)
)

if (!all(numeric_type_check)) {
  stop(
    "Some expected numeric variables are not numeric: ",
    paste(
      names(numeric_type_check)[
        !numeric_type_check
      ],
      collapse = ", "
    )
  )
}

factor_type_check <- vapply(
  df[categorical_columns],
  is.factor,
  logical(1)
)

if (!all(factor_type_check)) {
  stop(
    "Some expected categorical variables are not factors: ",
    paste(
      names(factor_type_check)[
        !factor_type_check
      ],
      collapse = ", "
    )
  )
}

if (!inherits(df$Date, "Date")) {
  stop("The Date variable was not stored as class Date.")
}

date_column <- df$Date
