# ==============================================================================
# Demographic and Clinical Data Preprocessing
# ==============================================================================
#
# Description:
#   This script loads and preprocesses the clinical dataset used for the
#   TMT analysis. Missing TMT measurements are supplemented using measurements
#   from the opposite side or, where necessary, from the corresponding
#   measurements of the other raters.
#
#   The script additionally:
#     - calculates patient age at the time of MRI examination,
#     - prepares diagnosis categories,
#     - generates descriptive demographic visualizations.
#
# Input:
#   Excel file containing the clinical and TMT measurements.
#
# Required variables:
#   - RC
#   - Gender
#   - Datum.MRI
#   - TMT.SIN...mm....MEDIK
#   - TMT.DX...mm....MEDIK
#   - TMT.SIN...mm....RADIOLOG
#   - TMT.DX...mm....RADIOLOG
#   - L.sin.Kaplanova
#   - L.dx.Kaplanova
#   - Diag (optional)
#
# Output:
#   - Preprocessed dataset `data`
#   - Patient age at MRI
#   - Imputed TMT measurements
#   - Age-by-diagnosis boxplot
#   - Age distribution by diagnosis and gender
#
# Requirements:
#   R >= 4.0
#
# ==============================================================================

#author: Jana Schwarzerova

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

required_packages <- c(
    "readxl",
    "dplyr",
    "ggplot2"
)

missing_packages <- required_packages[
    !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
    stop(
        "The following packages are not installed: ",
        paste(missing_packages, collapse = ", ")
    )
}

invisible(
    lapply(
        required_packages,
        library,
        character.only = TRUE
    )
)


# ------------------------------------------------------------------------------
# 2. Load input data
# ------------------------------------------------------------------------------

# The input file should be placed in the project data directory.
#
# Example:
# data <- readxl::read_excel("data/Reguli_et_al_input.xlsx")

input_file <- "data/Reguli_et_al_input.xlsx"

if (!file.exists(input_file)) {
    stop(
        "Input file not found: ",
        input_file,
        "\nPlease update `input_file` to the location of the input dataset."
    )
}

data <- readxl::read_excel(input_file)

# Standardize column names for reproducible access in R.
names(data) <- make.names(names(data))


# ------------------------------------------------------------------------------
# 3. Data type conversion
# ------------------------------------------------------------------------------

data <- data %>%
    mutate(
        Gender = as.factor(Gender),
        RC = as.character(RC),
        Datum.MRI = as.Date(Datum.MRI)
    )


# ------------------------------------------------------------------------------
# 4. TMT preprocessing
# ------------------------------------------------------------------------------

# Missing measurements for the left/right side are supplemented using the
# corresponding measurement from the opposite side.
#
# For the third rater, missing measurements are first supplemented from the
# opposite side and subsequently estimated from the corresponding measurements
# of raters 1 and 2 when both are available.

data <- data %>%
    mutate(
        
        # --------------------------------------------------------------------------
        # Rater 1: Medical measurement
        # --------------------------------------------------------------------------
        
        TMT.SIN...mm....MEDIK = if_else(
            is.na(TMT.SIN...mm....MEDIK),
            TMT.DX...mm....MEDIK,
            TMT.SIN...mm....MEDIK
        ),
        
        TMT.DX...mm....MEDIK = if_else(
            is.na(TMT.DX...mm....MEDIK),
            TMT.SIN...mm....MEDIK,
            TMT.DX...mm....MEDIK
        ),
        
        
        # --------------------------------------------------------------------------
        # Rater 2: Radiological measurement
        # --------------------------------------------------------------------------
        
        TMT.SIN...mm....RADIOLOG = if_else(
            is.na(TMT.SIN...mm....RADIOLOG),
            TMT.DX...mm....RADIOLOG,
            TMT.SIN...mm....RADIOLOG
        ),
        
        TMT.DX...mm....RADIOLOG = if_else(
            is.na(TMT.DX...mm....RADIOLOG),
            TMT.SIN...mm....RADIOLOG,
            TMT.DX...mm....RADIOLOG
        ),
        
        
        # --------------------------------------------------------------------------
        # Rater 3: Kaplanova
        # --------------------------------------------------------------------------
        
        L.dx.Kaplanova = if_else(
            is.na(L.dx.Kaplanova),
            L.sin.Kaplanova,
            L.dx.Kaplanova
        ),
        
        L.sin.Kaplanova = if_else(
            is.na(L.sin.Kaplanova),
            L.dx.Kaplanova,
            L.sin.Kaplanova
        ),
        
        # If both measurements from rater 3 are unavailable, estimate the
        # corresponding side from raters 1 and 2.
        L.sin.Kaplanova = if_else(
            is.na(L.sin.Kaplanova),
            rowMeans(
                cbind(
                    TMT.SIN...mm....MEDIK,
                    TMT.SIN...mm....RADIOLOG
                ),
                na.rm = TRUE
            ),
            L.sin.Kaplanova
        ),
        
        L.dx.Kaplanova = if_else(
            is.na(L.dx.Kaplanova),
            rowMeans(
                cbind(
                    TMT.DX...mm....MEDIK,
                    TMT.DX...mm....RADIOLOG
                ),
                na.rm = TRUE
            ),
            L.dx.Kaplanova
        )
    )


# ------------------------------------------------------------------------------
# 5. Calculate patient age at MRI examination
# ------------------------------------------------------------------------------

# The patient identifier (`RC`) encodes the date of birth in YYMMDD format.
# For female patients, 50 is added to the month component according to the
# coding convention used in the source dataset.
#
# Age is calculated at the date of MRI examination.


# Extract birth year.

birth_year <- as.numeric(
    paste0(
        "19",
        substr(data$RC, 1, 2)
    )
)


# Extract birth month.

birth_month <- as.numeric(
    substr(data$RC, 3, 4)
)

birth_month[data$Gender == "F"] <-
    birth_month[data$Gender == "F"] - 50


# Extract birth day.

birth_day <- as.numeric(
    substr(data$RC, 5, 6)
)


# Construct dates of birth.

date_of_birth <- as.Date(
    paste(
        birth_year,
        birth_month,
        birth_day,
        sep = "-"
    )
)


# Calculate age at MRI.

data <- data %>%
    mutate(
        Age = as.integer(
            floor(
                as.numeric(
                    difftime(
                        Datum.MRI,
                        date_of_birth,
                        units = "days"
                    )
                ) / 365.25
            )
        )
    )


# ------------------------------------------------------------------------------
# 6. Order diagnosis categories
# ------------------------------------------------------------------------------

# If a diagnosis variable is available, place "Other" as the final category.

if ("Diag" %in% names(data)) {
    
    diagnosis_levels <- unique(
        as.character(data$Diag)
    )
    
    diagnosis_levels <- c(
        setdiff(diagnosis_levels, "Other"),
        "Other"
    )
    
    data$Diag <- factor(
        data$Diag,
        levels = diagnosis_levels
    )
}


# ------------------------------------------------------------------------------
# 7. Demographic visualization
# ------------------------------------------------------------------------------

if ("Diag" %in% names(data)) {
    
    # --------------------------------------------------------------------------
    # Age distribution by diagnosis and gender
    # --------------------------------------------------------------------------
    
    age_by_diagnosis_plot <- ggplot(
        data,
        aes(
            x = Diag,
            y = Age,
            fill = Gender
        )
    ) +
        geom_boxplot(
            position = position_dodge(
                width = 0.75
            )
        ) +
        labs(
            x = "Diagnosis",
            y = "Age",
            fill = "Gender"
        ) +
        theme_minimal()
    
    print(age_by_diagnosis_plot)
    
    
    # --------------------------------------------------------------------------
    # Age histogram by diagnosis and gender
    # --------------------------------------------------------------------------
    
    age_distribution_plot <- ggplot(
        data,
        aes(
            x = Age,
            fill = Gender
        )
    ) +
        geom_histogram(
            position = "dodge",
            binwidth = 5,
            color = "black"
        ) +
        facet_wrap(
            ~ Diag
        ) +
        labs(
            x = "Age",
            y = "Number of Patients",
            fill = "Gender"
        ) +
        theme_minimal()
    
    print(age_distribution_plot)
}


# ------------------------------------------------------------------------------
# 8. End of analysis
# ------------------------------------------------------------------------------