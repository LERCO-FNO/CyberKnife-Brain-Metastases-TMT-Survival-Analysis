# ==============================================================================
# TMT Survival Analysis
# ==============================================================================
#
# Description:
#   This script performs Kaplan-Meier survival analysis based on TMT
#   measurements obtained from three independent raters.
#
#   For each rater, patients are divided into two groups according to the
#   median TMT measurement (Low vs High). Kaplan-Meier survival curves are
#   generated for:
#
#     1. The complete study population
#     2. Female patients
#     3. Male patients
#
#   Survival time is calculated from the date of MRI examination to the date
#   of death. The event indicator is defined as 1 for observed death events.
#
# Input:
#   A clinical dataset containing:
#     - Datum.MRI
#     - Umrti
#     - Gender
#     - TMT.SIN...mm....MEDIK
#     - TMT.DX...mm....MEDIK
#     - TMT.SIN...mm....RADIOLOG
#     - TMT.DX...mm....RADIOLOG
#     - L.sin.Kaplanova
#     - L.dx.Kaplanova
#
# Output:
#   - Kaplan-Meier survival curves for each rater
#   - Risk tables
#   - Log-rank p-values
#   - 95% confidence intervals
#   - Sex-stratified Kaplan-Meier curves
#
# Requirements:
#   R >= 4.0
#
# ==============================================================================

# author: Jana Schwarzerova

# ------------------------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------------------------

required_packages <- c(
    "dplyr",
    "survival",
    "survminer",
    "patchwork"
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

# The analysis assumes that the clinical dataset has already been loaded
# as `data`.
#
# Example:
# data <- readxl::read_excel("data/Reguli_et_al_input.xlsx")
#
# If column names require standardization:
# names(data) <- make.names(names(data))

if (!exists("data")) {
    stop(
        "Input dataset `data` was not found. ",
        "Please load the dataset before running the analysis."
    )
}


# ------------------------------------------------------------------------------
# 3. Prepare survival variables
# ------------------------------------------------------------------------------

data <- data %>%
    mutate(
        Datum.MRI = as.Date(Datum.MRI),
        Umrti = as.Date(Umrti)
    )


# ------------------------------------------------------------------------------
# 4. Calculate survival time
# ------------------------------------------------------------------------------

# Overall survival time is defined as the number of days between
# MRI examination and death.
#
# IMPORTANT:
# This implementation assumes that all patients experienced the event
# (death). If censored observations are present, `OS_event` must be defined
# according to the study-specific censoring information.

data <- data %>%
    mutate(
        OS_time = as.numeric(
            difftime(
                Umrti,
                Datum.MRI,
                units = "days"
            )
        ),
        OS_event = 1
    )


# ------------------------------------------------------------------------------
# 5. Calculate mean TMT measurement for each rater
# ------------------------------------------------------------------------------

# For each rater, the left- and right-side measurements are averaged to obtain
# a single patient-level TMT value.

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
    ) %>%
    select(
        OS_time,
        OS_event,
        Gender,
        TMT_avg_rater1,
        TMT_avg_rater2,
        TMT_avg_rater3
    )


# ------------------------------------------------------------------------------
# 6. Kaplan-Meier analysis by TMT median
# ------------------------------------------------------------------------------

plot_survival_by_rater <- function(
    df,
    tmt_col,
    rater_label
) {
    
    # Calculate the median TMT value.
    median_value <- median(
        df[[tmt_col]],
        na.rm = TRUE
    )
    
    # Divide patients into Low and High TMT groups.
    df <- df %>%
        mutate(
            TMT_group = if_else(
                .data[[tmt_col]] >= median_value,
                "High",
                "Low"
            ),
            TMT_group = factor(
                TMT_group,
                levels = c("Low", "High")
            )
        ) %>%
        filter(
            !is.na(OS_time),
            !is.na(OS_event),
            !is.na(TMT_group)
        )
    
    # Kaplan-Meier model.
    survival_fit <- survfit(
        Surv(OS_time, OS_event) ~ TMT_group,
        data = df
    )
    
    # Plot.
    survival_plot <- ggsurvplot(
        survival_fit,
        data = df,
        risk.table = TRUE,
        pval = TRUE,
        conf.int = TRUE,
        legend.title = paste0(
            "TMT ",
            rater_label
        ),
        legend.labs = c(
            "Low",
            "High"
        ),
        title = paste0(
            "Kaplan-Meier Survival Analysis: ",
            rater_label
        ),
        xlab = "Time (days)",
        ylab = "Survival probability"
    )
    
    return(survival_plot)
}


# ------------------------------------------------------------------------------
# 7. Kaplan-Meier curves for all raters
# ------------------------------------------------------------------------------

km_rater1 <- plot_survival_by_rater(
    df_surv,
    "TMT_avg_rater1",
    "Rater 1"
)

km_rater2 <- plot_survival_by_rater(
    df_surv,
    "TMT_avg_rater2",
    "Rater 2"
)

km_rater3 <- plot_survival_by_rater(
    df_surv,
    "TMT_avg_rater3",
    "Rater 3"
)


# ------------------------------------------------------------------------------
# 8. Display Kaplan-Meier curves
# ------------------------------------------------------------------------------

all_raters_plot <- (
    km_rater1$plot +
        km_rater2$plot +
        km_rater3$plot
) +
    plot_layout(
        ncol = 3
    )

print(all_raters_plot)


# ------------------------------------------------------------------------------
# 9. Sex-stratified Kaplan-Meier analysis
# ------------------------------------------------------------------------------

plot_survival_by_gender <- function(
    df,
    tmt_col,
    rater_label
) {
    
    median_value <- median(
        df[[tmt_col]],
        na.rm = TRUE
    )
    
    df <- df %>%
        mutate(
            TMT_group = if_else(
                .data[[tmt_col]] >= median_value,
                "High",
                "Low"
            ),
            TMT_group = factor(
                TMT_group,
                levels = c("Low", "High")
            )
        ) %>%
        filter(
            !is.na(OS_time),
            !is.na(OS_event),
            !is.na(TMT_group),
            !is.na(Gender)
        )
    
    gender_levels <- unique(
        as.character(df$Gender)
    )
    
    plots <- list()
    
    for (gender in gender_levels) {
        
        df_gender <- df %>%
            filter(
                as.character(Gender) == gender
            )
        
        # Skip groups with insufficient observations or without both TMT groups.
        if (
            nrow(df_gender) < 2 ||
            length(unique(df_gender$TMT_group)) < 2
        ) {
            next
        }
        
        survival_fit <- survfit(
            Surv(OS_time, OS_event) ~ TMT_group,
            data = df_gender
        )
        
        plots[[gender]] <- ggsurvplot(
            survival_fit,
            data = df_gender,
            risk.table = TRUE,
            pval = TRUE,
            conf.int = TRUE,
            legend.title = paste0(
                "TMT ",
                rater_label
            ),
            legend.labs = c(
                "Low",
                "High"
            ),
            title = paste0(
                "Kaplan-Meier Survival: ",
                rater_label,
                " (Gender: ",
                gender,
                ")"
            ),
            xlab = "Time (days)",
            ylab = "Survival probability"
        )
    }
    
    return(plots)
}


# ------------------------------------------------------------------------------
# 10. Generate sex-stratified curves
# ------------------------------------------------------------------------------

km_rater1_gender <- plot_survival_by_gender(
    df_surv,
    "TMT_avg_rater1",
    "Rater 1"
)

km_rater2_gender <- plot_survival_by_gender(
    df_surv,
    "TMT_avg_rater2",
    "Rater 2"
)

km_rater3_gender <- plot_survival_by_gender(
    df_surv,
    "TMT_avg_rater3",
    "Rater 3"
)


# ------------------------------------------------------------------------------
# 11. Display sex-stratified curves
# ------------------------------------------------------------------------------

# Helper function for displaying available gender-specific plots.

combine_gender_plots <- function(
    plot_list,
    ncol = 2
) {
    
    if (length(plot_list) == 0) {
        return(NULL)
    }
    
    plot_objects <- lapply(
        plot_list,
        function(x) x$plot
    )
    
    wrap_plots(
        plot_objects,
        ncol = ncol
    )
}


# Rater 1
if (length(km_rater1_gender) > 0) {
    print(
        combine_gender_plots(
            km_rater1_gender
        )
    )
}


# Rater 2
if (length(km_rater2_gender) > 0) {
    print(
        combine_gender_plots(
            km_rater2_gender
        )
    )
}


# Rater 3
if (length(km_rater3_gender) > 0) {
    print(
        combine_gender_plots(
            km_rater3_gender
        )
    )
}


# ------------------------------------------------------------------------------
# 12. Combined sex-stratified visualization
# ------------------------------------------------------------------------------

gender_plot_rater1 <- combine_gender_plots(
    km_rater1_gender
)

gender_plot_rater2 <- combine_gender_plots(
    km_rater2_gender
)

gender_plot_rater3 <- combine_gender_plots(
    km_rater3_gender
)

available_plots <- Filter(
    Negate(is.null),
    list(
        gender_plot_rater1,
        gender_plot_rater2,
        gender_plot_rater3
    )
)

if (length(available_plots) > 0) {
    
    combined_gender_plots <- wrap_plots(
        available_plots,
        ncol = 1
    )
    
    print(combined_gender_plots)
}


# ------------------------------------------------------------------------------
# End of analysis
# ------------------------------------------------------------------------------
