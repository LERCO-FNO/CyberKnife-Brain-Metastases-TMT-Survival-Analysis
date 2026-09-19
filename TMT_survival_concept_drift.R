# ==============================================================================
# Survival Analysis and Concept Drift Detection Based on TMT Measurements
# ==============================================================================
#
# Description:
#   This script performs survival analysis and longitudinal concept drift
#   detection using TMT measurements from multiple raters. Drift is evaluated
#   across predefined temporal segments using Cox proportional hazards,
#   random survival forest (RSF), and survival support vector machine (SVM)
#   models.
#
#   Age and gender are retained as patient-level characteristics for subsequent
#   characterization of patients with detected drift.
#
# Input:
#   A data.frame/tibble named `data` containing at least the following variables:
#     - Gender
#     - Age
#     - Datum.MRI
#     - Umrti
#     - TMT.SIN...mm....MEDIK
#     - TMT.DX...mm....MEDIK
#     - TMT.SIN...mm....RADIOLOG
#     - TMT.DX...mm....RADIOLOG
#     - L.sin.Kaplanova
#     - L.dx.Kaplanova
#
# Output:
#   - Patient-level concept drift results for Cox, RSF, and SVM models
#   - Descriptive summaries of age and gender among patients with detected drift
#   - Optional visualization of the age distribution
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
    "lubridate",
    "ggplot2",
    "survival",
    "randomForestSRC",
    "survivalsvm"
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

invisible(lapply(required_packages, library, character.only = TRUE))


# ------------------------------------------------------------------------------
# 2. Input data
# ------------------------------------------------------------------------------

# The analysis assumes that the input dataset is already loaded as `data`.
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

data <- as.data.frame(data)


# ------------------------------------------------------------------------------
# 3. Data type conversion
# ------------------------------------------------------------------------------

data <- data %>%
    mutate(
        Gender = as.factor(Gender),
        Datum.MRI = as.Date(Datum.MRI),
        RC = as.character(RC)
    )


# ------------------------------------------------------------------------------
# 4. Preprocess TMT measurements
# ------------------------------------------------------------------------------

# Missing left/right measurements are imputed using the measurement from
# the opposite side when available.
#
# For the third rater, missing left/right measurements are subsequently
# estimated from the corresponding radiologist and medical measurements.

data <- data %>%
    mutate(
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
# 5. Calculate mean TMT measurements for each rater
# ------------------------------------------------------------------------------

df_surv <- data %>%
    mutate(
        TMT_avg_rater1 = rowMeans(
            cbind(
                as.numeric(TMT.SIN...mm....MEDIK),
                as.numeric(TMT.DX...mm....MEDIK)
            ),
            na.rm = TRUE
        ),
        TMT_avg_rater2 = rowMeans(
            cbind(
                as.numeric(TMT.SIN...mm....RADIOLOG),
                as.numeric(TMT.DX...mm....RADIOLOG)
            ),
            na.rm = TRUE
        ),
        TMT_avg_rater3 = rowMeans(
            cbind(
                as.numeric(L.sin.Kaplanova),
                as.numeric(L.dx.Kaplanova)
            ),
            na.rm = TRUE
        )
    )


# ------------------------------------------------------------------------------
# 6. Calculate overall survival
# ------------------------------------------------------------------------------

# Overall survival time is calculated as the number of days between
# MRI examination and death.
#
# IMPORTANT:
# The current implementation assumes that all observations represent events
# (OS_event = 1). If censored observations are present, OS_event must be
# defined accordingly.

data <- data %>%
    mutate(
        OS_time = as.numeric(Umrti - Datum.MRI),
        OS_event = 1
    )

df_surv <- df_surv %>%
    mutate(
        time = data$OS_time,
        status = data$OS_event,
        Gender = data$Gender,
        Age = data$Age,
        Datum.MRI = data$Datum.MRI
    )


# ------------------------------------------------------------------------------
# 7. Define temporal segments
# ------------------------------------------------------------------------------

# Temporal segments are used to evaluate changes in model predictions
# between consecutive periods.

df_surv <- df_surv %>%
    mutate(
        year_segment = cut(
            Datum.MRI,
            breaks = as.Date(
                c(
                    "2010-01-01",
                    "2014-12-31",
                    "2018-12-31",
                    "2022-12-31",
                    "2024-08-31"
                )
            ),
            labels = c(
                "2010-2014",
                "2015-2018",
                "2019-2022",
                "2023-2024"
            ),
            include.lowest = TRUE
        )
    )


# ------------------------------------------------------------------------------
# 8. Concept drift detection
# ------------------------------------------------------------------------------

detect_drift_per_patient <- function(
    df,
    model_type = c("Cox", "RSF", "SVM"),
    min_patients = 10,
    threshold = 0.1,
    ntree = 500
) {
    
    model_type <- match.arg(model_type)
    
    drift_patients <- data.frame(
        PatientRow = integer(),
        segment = character(),
        Gender = character(),
        Age = numeric(),
        stringsAsFactors = FALSE
    )
    
    segments <- sort(unique(na.omit(df$year_segment)))
    
    if (length(segments) < 2) {
        warning("At least two temporal segments are required.")
        return(drift_patients)
    }
    
    for (i in seq_along(segments)[-1]) {
        
        train_seg <- segments[i - 1]
        test_seg <- segments[i]
        
        df_train <- df %>%
            filter(year_segment == train_seg)
        
        df_test <- df %>%
            filter(year_segment == test_seg)
        
        # Skip comparisons when either segment contains too few observations.
        if (
            nrow(df_train) < min_patients ||
            nrow(df_test) < min_patients
        ) {
            next
        }
        
        # --------------------------------------------------------------------------
        # Fit survival model
        # --------------------------------------------------------------------------
        
        if (model_type == "Cox") {
            
            model <- coxph(
                Surv(time, status) ~
                    TMT_avg_rater1 +
                    TMT_avg_rater2 +
                    TMT_avg_rater3,
                data = df_train
            )
            
            pred_test <- predict(
                model,
                newdata = df_test
            )
            
        } else if (model_type == "RSF") {
            
            model <- rfsrc(
                Surv(time, status) ~
                    TMT_avg_rater1 +
                    TMT_avg_rater2 +
                    TMT_avg_rater3,
                data = df_train,
                ntree = ntree
            )
            
            pred_test <- predict(
                model,
                newdata = df_test
            )$predicted
            
        } else if (model_type == "SVM") {
            
            model <- survivalsvm(
                Surv(time, status) ~
                    TMT_avg_rater1 +
                    TMT_avg_rater2 +
                    TMT_avg_rater3,
                data = df_train,
                kernel = "lin_kernel",
                gamma.mu = 1,
                opt.meth = "quadprog"
            )
            
            pred_test <- as.numeric(
                predict(
                    model,
                    newdata = df_test
                )$predicted
            )
        }
        
        # --------------------------------------------------------------------------
        # Standardize predictions within the test segment
        # --------------------------------------------------------------------------
        
        pred_scaled <- as.numeric(scale(pred_test))
        
        # --------------------------------------------------------------------------
        # Identify patients exceeding the drift threshold
        # --------------------------------------------------------------------------
        
        drift_idx <- which(abs(pred_scaled) > threshold)
        
        if (length(drift_idx) > 0) {
            
            drift_patients <- rbind(
                drift_patients,
                data.frame(
                    PatientRow = as.integer(rownames(df_test)[drift_idx]),
                    segment = as.character(test_seg),
                    Gender = as.character(df_test$Gender[drift_idx]),
                    Age = df_test$Age[drift_idx],
                    stringsAsFactors = FALSE
                )
            )
        }
    }
    
    drift_patients
}


# ------------------------------------------------------------------------------
# 9. Run concept drift detection
# ------------------------------------------------------------------------------

drift_cox_patients <- detect_drift_per_patient(
    df_surv,
    model_type = "Cox"
)

drift_rsf_patients <- detect_drift_per_patient(
    df_surv,
    model_type = "RSF"
)

drift_svm_patients <- detect_drift_per_patient(
    df_surv,
    model_type = "SVM"
)


# ------------------------------------------------------------------------------
# 10. Summarize patients with detected drift
# ------------------------------------------------------------------------------

# Gender distribution
gender_summary_cox <- table(
    drift_cox_patients$Gender
)

print(gender_summary_cox)


# Age distribution
age_summary_cox <- summary(
    drift_cox_patients$Age
)

print(age_summary_cox)


# ------------------------------------------------------------------------------
# 11. Visualize age distribution
# ------------------------------------------------------------------------------

age_plot_cox <- ggplot(
    drift_cox_patients,
    aes(
        x = Age,
        fill = Gender
    )
) +
    geom_histogram(
        binwidth = 5,
        alpha = 0.7,
        position = "dodge"
    ) +
    theme_minimal() +
    labs(
        title = "Age Distribution of Patients with Detected Concept Drift",
        x = "Age",
        y = "Number of Patients",
        fill = "Gender"
    )

print(age_plot_cox)


# ------------------------------------------------------------------------------
# End of analysis
# ------------------------------------------------------------------------------

