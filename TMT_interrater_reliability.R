# ==============================================================================
# TMT Inter-Rater Reliability and Agreement Analysis
# ==============================================================================
#
# Description:
#   This script performs preprocessing and inter-rater reliability analysis
#   of TMT measurements obtained from three independent raters.
#
#   The analysis includes:
#     1. TMT preprocessing and outlier removal
#     2. Calculation of TMT summary measurements
#     3. Overall inter-rater reliability using the intraclass correlation
#        coefficient (ICC)
#     4. Pairwise ICC analysis between individual raters
#     5. Bland-Altman agreement analysis
#
# Input:
#   A data.frame/tibble containing TMT measurements for three raters.
#
#   Required variables:
#     - MEAN1
#     - SIN
#     - DX
#     - PatientID
#     - TMT.SIN...mm....MEDIK
#     - TMT.DX...mm....MEDIK
#     - TMT.SIN...mm....RADIOLOG
#     - TMT.DX...mm....RADIOLOG
#     - L.sin.Kaplanova
#     - L.dx.Kaplanova
#
# Expected input objects:
#   - `df_min`: preprocessed dataset containing the variable `MEAN1`
#   - `data`: dataset containing individual TMT measurements
#
# Output:
#   - Cleaned TMT dataset
#   - Overall ICC across all three raters
#   - Pairwise ICC estimates
#   - Bland-Altman bias and limits of agreement
#   - Bland-Altman plots
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
    "dplyr",
    "tidyr",
    "irr",
    "psych",
    "BlandAltmanLeh"
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
# 2. Input data
# ------------------------------------------------------------------------------

# The analysis assumes that the required datasets have already been loaded.
#
# Example:
# data <- readxl::read_excel("path/to/input_data.xlsx")
#
# If required, standardize column names:
# names(data) <- make.names(names(data))

if (!exists("data")) {
    stop(
        "Input dataset `data` was not found. ",
        "Please load the dataset before running the analysis."
    )
}

if (!exists("df_min")) {
    stop(
        "Input dataset `df_min` was not found. ",
        "Please load the preprocessed dataset before running the analysis."
    )
}


# ------------------------------------------------------------------------------
# 3. TMT preprocessing
# ------------------------------------------------------------------------------

# The primary TMT summary measure is based on the variable `MEAN1`.

df_min <- df_min %>%
    mutate(
        TMT_avg = MEAN1
    )


# ------------------------------------------------------------------------------
# 4. Outlier detection using the IQR method
# ------------------------------------------------------------------------------

# Outliers are identified using the conventional 1.5 × IQR rule.

Q1 <- quantile(
    df_min$TMT_avg,
    probs = 0.25,
    na.rm = TRUE
)

Q3 <- quantile(
    df_min$TMT_avg,
    probs = 0.75,
    na.rm = TRUE
)

IQR_value <- Q3 - Q1

cat("TMT interquartile range (IQR):", IQR_value, "\n")


# Remove observations outside the 1.5 × IQR range.

df_clean <- df_min %>%
    filter(
        TMT_avg >= Q1 - 1.5 * IQR_value,
        TMT_avg <= Q3 + 1.5 * IQR_value
    )


# Retain observations with complete left- and right-side measurements.

df_clean <- df_clean %>%
    filter(
        !is.na(SIN),
        !is.na(DX),
        !is.na(TMT_avg)
    )


# ------------------------------------------------------------------------------
# 5. Helper function for PCA input preparation
# ------------------------------------------------------------------------------

# This function creates a standardized long-format dataset containing
# left/right TMT measurements for a specified rater.

build_pca_input <- function(
    sin_name,
    dx_name,
    rater_label,
    data
) {
    
    sin <- as.numeric(data[[sin_name]])
    dx <- as.numeric(data[[dx_name]])
    id <- data$PatientID
    
    pca_input <- data.frame(
        ID = id,
        SIN = sin,
        DX = dx,
        Rater = rater_label
    ) %>%
        tidyr::drop_na()
    
    pca_input
}


# ------------------------------------------------------------------------------
# 6. Prepare data for inter-rater reliability analysis
# ------------------------------------------------------------------------------

# For each rater, the mean TMT measurement is calculated from the left
# and right measurements.

df_icc <- data.frame(
    rater1 = rowMeans(
        cbind(
            as.numeric(data$TMT.SIN...mm....MEDIK),
            as.numeric(data$TMT.DX...mm....MEDIK)
        ),
        na.rm = TRUE
    ),
    
    rater2 = rowMeans(
        cbind(
            as.numeric(data$TMT.SIN...mm....RADIOLOG),
            as.numeric(data$TMT.DX...mm....RADIOLOG)
        ),
        na.rm = TRUE
    ),
    
    rater3 = rowMeans(
        cbind(
            as.numeric(data$L.sin.Kaplanova),
            as.numeric(data$L.dx.Kaplanova)
        ),
        na.rm = TRUE
    )
)


# Retain complete cases for the reliability analysis.

df_icc <- df_icc %>%
    tidyr::drop_na()


# ------------------------------------------------------------------------------
# 7. Overall inter-rater reliability
# ------------------------------------------------------------------------------

# Overall ICC is calculated across all three raters using psych::ICC.

cat("\n")
cat("============================================================\n")
cat("Overall ICC - All Raters\n")
cat("============================================================\n")

overall_icc <- psych::ICC(
    df_icc
)

print(overall_icc)


# ------------------------------------------------------------------------------
# 8. Pairwise inter-rater reliability
# ------------------------------------------------------------------------------

# Pairwise ICC is calculated using a two-way agreement model
# with single-rater units.

rater_pairs <- list(
    c("rater1", "rater2"),
    c("rater1", "rater3"),
    c("rater2", "rater3")
)


pairwise_icc_results <- list()

for (pair in rater_pairs) {
    
    df_pair <- df_icc[, pair]
    
    icc_result <- irr::icc(
        df_pair,
        model = "twoway",
        type = "agreement",
        unit = "single"
    )
    
    pair_name <- paste(
        pair[1],
        "vs",
        pair[2]
    )
    
    pairwise_icc_results[[pair_name]] <- icc_result
    
    cat("\n")
    cat("============================================================\n")
    cat("Pairwise ICC:", pair_name, "\n")
    cat("============================================================\n")
    
    print(icc_result)
}


# ------------------------------------------------------------------------------
# 9. Bland-Altman analysis
# ------------------------------------------------------------------------------

# The Bland-Altman analysis estimates:
#   - mean difference (bias)
#   - lower limit of agreement
#   - upper limit of agreement
#
# Limits of agreement are calculated as:
#
#   bias ± 1.96 × SD of the paired differences

bland_altman_results <- function(x, y) {
    
    difference <- x - y
    
    bias <- mean(
        difference,
        na.rm = TRUE
    )
    
    sd_difference <- sd(
        difference,
        na.rm = TRUE
    )
    
    loa_lower <- bias - 1.96 * sd_difference
    loa_upper <- bias + 1.96 * sd_difference
    
    list(
        bias = bias,
        loa_lower = loa_lower,
        loa_upper = loa_upper
    )
}


# ------------------------------------------------------------------------------
# 10. Run Bland-Altman analysis for all rater pairs
# ------------------------------------------------------------------------------

bland_altman_summary <- list()

for (pair in rater_pairs) {
    
    x <- df_icc[[pair[1]]]
    y <- df_icc[[pair[2]]]
    
    result <- bland_altman_results(
        x,
        y
    )
    
    pair_name <- paste(
        pair[1],
        "vs",
        pair[2]
    )
    
    bland_altman_summary[[pair_name]] <- result
    
    cat("\n")
    cat("============================================================\n")
    cat("Bland-Altman Analysis:", pair_name, "\n")
    cat("============================================================\n")
    
    cat(
        "Bias:",
        round(result$bias, 2),
        "mm\n"
    )
    
    cat(
        "Lower limit of agreement:",
        round(result$loa_lower, 2),
        "mm\n"
    )
    
    cat(
        "Upper limit of agreement:",
        round(result$loa_upper, 2),
        "mm\n"
    )
    
    # Bland-Altman plot
    bland.altman.plot(
        x,
        y,
        main = paste(
            "Bland-Altman:",
            pair_name
        )
    )
}


# ------------------------------------------------------------------------------
# 11. End of analysis
# ------------------------------------------------------------------------------

