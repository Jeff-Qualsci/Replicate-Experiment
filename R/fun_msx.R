# Functionalized version of msx.R
library(tidyverse)
library(patchwork)
library(ggrepel)

# Test data generation functions ---------------------------------------
# set.seed(12252023) #For reproducible examples 

# Generate Potency test data for Replicate-Experiment
msr_data <- function(SmplNum, TstMSR, Shift = 1) {

  TstSD <- log10(TstMSR) / 2 # convert MSR to the SD on the log10 scale
  Shift <- log10(Shift)  # will be added to the 2nd replicate to simulate a fold-shift
  Sample <- sample(c(1000000 : 9999999), SmplNum)
  TrueMeas <- runif(SmplNum, -3, 3) # True log10(Potency) values
  ExpData <- map(TrueMeas, \(m) rnorm(2, m, TstSD)) # Generate 2 random values for each TrueMeas

  TstData <- tibble(Sample, ExpData) %>% 
    unnest_wider(col = ExpData, names_sep = '_') %>%
    rename(Exp1 = ExpData_1, Exp2 = ExpData_2) %>%
    mutate(Exp2 = Exp2 + Shift,
            across(-Sample, ~ 10 ^ .x))
}

# Generate Efficacy test data for Replicate-Experiment- consistent SD
msd_data <- function(SmplNum, TstMSD, Shift = 0) {

  TstSD <- TstMSD /2
  Sample <- sample(c(1000000 : 9999999), SmplNum)
  TrueMeas <- runif(SmplNum, -20, 120)
  ExpData <- map(TrueMeas, \(m) rnorm(2, m, TstSD))

  TstData <- tibble(Sample, ExpData) %>%
    unnest_wider(col = ExpData, names_sep = '_') %>%
    rename(Exp1 = ExpData_1, Exp2 = ExpData_2) %>%
    mutate(Exp2 = Exp2 + Shift)
}

# Generate Efficacy test data for Replicate-Experiment - consistent cv
msd_cv_data <- function(SmplNum, cv, Shift = 0) {

  Sample <- sample(c(1000000 : 9999999), SmplNum)
  TrueMeas <- runif(SmplNum, -10, 110)
  ExpData <- map(TrueMeas, \(m) rnorm(2, m, abs(cv*m)))

  TstData <- tibble(Sample, ExpData) %>% 
    unnest_wider(col = ExpData, names_sep = '_') %>%
    rename(Exp1 = ExpData_1, Exp2 = ExpData_2) %>%
    mutate(Exp2 = Exp2 + Shift)
}

# Generate inactive efficacy test data for Replicate-Experiment ----------------
msd_inact <- function(SmplNum, TstMSD, Shift = 0) {

  TstSD <- TstMSD /(2 * sqrt(2))
  Sample <- sample(c(1000000 : 9999999), SmplNum)
  TrueMeas <- runif(SmplNum, -20, 20)
  ExpData <- map(TrueMeas, \(m) rnorm(2, m, TstSD))

  TstData <- tibble(Sample, ExpData) %>% 
    unnest_wider(col = ExpData, names_sep = '_') %>%
    rename(Exp1 = ExpData_1, Exp2 = ExpData_2) %>%
    mutate(Exp2 = Exp2 + Shift)
}

# Generate active efficacy test data for Replicate-Experiment ------------
msd_act <- function(SmplNum, TstMSD, Shift = 0) {

  TstSD <- TstMSD /(2 * sqrt(2))
  Sample <- sample(c(1000000 : 9999999), SmplNum)
  TrueMeas <- runif(SmplNum, 50, 120)
  ExpData <- map(TrueMeas, \(m) rnorm(2, m, TstSD))

  TstData <- tibble(Sample, ExpData) %>%
    unnest_wider(col = ExpData, names_sep = '_') %>%
    rename(Exp1 = ExpData_1, Exp2 = ExpData_2) %>%
    mutate(Exp2 = Exp2 + Shift)
}

# Replicate-Experiment Analysis Functions -------------------------------

# Summary Stats
repexp.stats <- function(df) {
  rSpearman <- cor(x = df[["Exp1"]], y = df[["Exp2"]], method = 'spearman')

  summarise(df,
            n = n(),
            t = qt(0.975, n - 1),
            MeanDiff = mean(Difference),
            StdDev = sd(Difference),
            MSD = 2 * t * StdDev,
            UDL = MeanDiff + (t * StdDev/sqrt(n)),
            LDL = MeanDiff - (t * StdDev/sqrt(n)),
            ULSA = MeanDiff + (3 * StdDev),
            LLSA = MeanDiff- (3 * StdDev)) %>%
    mutate(r = rSpearman,
           r2 = r ^ 2) %>%
    select(-t, -StdDev)
}

# Mean Difference Plot
mdplot <- function(Data, Stats) {
  MSDLabel <- paste("MSD =", format(Stats[["MSD"]], digits = 2))

  ggplot(Data, aes(x = Mean, y = Difference)) +
    geom_point(shape = Data[["Outlier"]], size = 3) +
    geom_text_repel(aes(label = Label), na.rm = TRUE, hjust = -.5) +
    geom_hline(yintercept = Stats[["MeanDiff"]], color = "mediumblue") +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = Stats[["UDL"]], color = "mediumblue", linetype = "dashed") +
    geom_text(x = 1, y = Stats[["UDL"]], label = "Upper Difference Limit", color = "mediumblue", vjust = -0.3) +
    geom_hline(yintercept = Stats[["LDL"]], color = "mediumblue", linetype = "dashed") +
    geom_text(x = 1, y = Stats[["LDL"]], label = "Lower Difference Limit", color = "mediumblue", vjust = 1) +
    geom_hline(yintercept = Stats[["ULSA"]], color = "#D55E00", linetype = "dotdash") +
    geom_text(x = 1, y = Stats[["ULSA"]], label = "Upper Limit of Agreement", color = "#D55E00", vjust = -0.3) +
    geom_hline(yintercept = Stats[["LLSA"]], color = "#D55E00", linetype = "dotdash") + 
    geom_text(x = 1, y = Stats[["LLSA"]], label = "Lower Limit of Agreement", color = "#D55E00", vjust = 1) +
    theme_minimal() +
    labs(title = "Difference vs Mean Efficacy",
         subtitle = MSDLabel,
         x = "Mean Efficacy", y = "Difference")
}

# Mean Ratio Plot
mrplot <- function(Data, Stats) {
  MSRLabel <- paste("MSR =", Stats[["MSR"]])

  ggplot(Data, aes(x =  GeometricMean, y = Ratio)) +
    geom_point(shape = Data[["Outlier"]], size = 3) +
    geom_text_repel(aes(label = Label), na.rm = TRUE) +
    geom_hline(yintercept = Stats[["MeanRatio"]], color = "mediumblue") +
    geom_text(x = log10(max(Data[["GeometricMean"]]))-0.5, y = log10(Stats[["MeanRatio"]]), label = 'Mean Ratio', color = 'mediumblue', vjust = 1) +
    geom_hline(yintercept = 1, color = "black") +
    geom_hline(yintercept = Stats[["URL"]], color = "mediumblue", linetype = "dashed") +
    geom_text(x = log10(max(Data[["GeometricMean"]]))-0.5, y = log10(Stats[["URL"]]), label = "Upper Ratio Limit", color = "mediumblue", vjust = -0.3) +
    geom_hline(yintercept = Stats[["LRL"]], color = "mediumblue", linetype = "dashed") +
    geom_text(x = log10(max(Data[["GeometricMean"]]))-0.5, y = log10(Stats[["LRL"]]), label = "Lower Ratio Limit", color = "mediumblue", vjust = 1) +
    geom_hline(yintercept = Stats[["ULSA"]], color = "#D55E00", linetype = "dotdash") +
    geom_text(x = log10(max(Data[["GeometricMean"]]))-0.5, y = log10(Stats[["ULSA"]]), label = "Upper Limit of Agreement", color = "#D55E00", vjust = -0.3) +
    geom_hline(yintercept = Stats[["LLSA"]], color = "#D55E00", linetype = "dotdash") +
    geom_text(x = log10(max(Data[["GeometricMean"]]))-0.5, y = log10(Stats[["LLSA"]]), label = "Lower Limit of Agreement", color = "#D55E00", vjust = 1) +
    theme_minimal() +
    scale_x_continuous(trans='log10') +
    scale_y_continuous(trans='log10') +
    labs(title = "Potency Ratio vs Geometric Mean Potency",
         subtitle = MSRLabel,
         x = "Geometric Mean Potency", y = "Ratio")
}

# R1R2 Plot - Run1/Run2 with correlation.
r1r2plot <- function(Data, Stats) {
  ggplot(Data,aes(x = Exp1, y = Exp2)) +
      geom_point() +
      geom_text_repel(aes(label = Label), na.rm = TRUE, hjust = -.5) +
      geom_smooth(method = "lm", se = TRUE, linetype = "dashed", linewidth = 2) +
      geom_abline(slope = 1) +
      labs(title = "Correlation Exp1 vs Exp2",
           subtitle = paste("Concordance Correlation r =", Stats[["r"]]),
           x = "Exp1",
           y = "Exp2") +
      theme_minimal()

}

# Replicate-Experiment Efficacy --------------------------------------------
repexp.efficacy<- function(df) {

  RepExp_Data <- df %>%
    mutate(Mean = (Exp1 + Exp2) / 2,
            Difference = Exp1 - Exp2)

  RepExp_Stats <- repexp.stats(RepExp_Data) %>%
    mutate(across(-n, \(x) signif(x, digits = 3)))

  # Flag data pairs with potential outliers (MeasDiff outside of limits of agreement)

  RepExp_Data <- RepExp_Data %>% 
    mutate(across(-Sample,  \(x) signif(x, digits = 3)),
            Outlier = Difference > RepExp_Stats[["ULSA"]] | Difference < RepExp_Stats[["LLSA"]],
            Label = if_else(Outlier, Sample, NA))

  MeanDifferencePlot <- mdplot(RepExp_Data, RepExp_Stats)

  R1R2CorrelationPlot <- r1r2plot(RepExp_Data, RepExp_Stats)

  list(Data = RepExp_Data, Stats = RepExp_Stats, MDPlot = MeanDifferencePlot, CorrPlot = R1R2CorrelationPlot)
}

# Replicate-Experiment Potency -------------------------------
repexp.potency <- function(df) {

  RepExp_Data <- df %>%
    mutate(across(-Sample, log10),
            Mean = (Exp1 + Exp2) / 2,
            Difference = Exp1 - Exp2)

  RepExp_Stats <- repexp.stats(RepExp_Data)

  # Flag data pairs with potential outliers (MeasDiff outside of limits of agreement)

  RepExp_Data <- RepExp_Data %>%
    mutate(Outlier = Difference > RepExp_Stats[["ULSA"]]| Difference < RepExp_Stats[["LLSA"]],
            Label = if_else(Outlier, Sample, NA),
            across(-c(1, 6, 7), ~10 ^ .x),
            across(-c(1, 6, 7),  \(x) signif(x, digits = 3))) %>% 
    rename(GeometricMean = Mean, Ratio = Difference)

  RepExp_Stats <- RepExp_Stats %>%
    rename(MSR = MSD,
           MeanRatio = MeanDiff,
           URL = UDL,
           LRL = LDL) %>%
    mutate(across(-c(1, 8, 9), ~10 ^ .x),
           across(-c(1), \(x) signif(x, digits = 3)))

  # MSRn Table
  n <- c(1:6)
  s <- log10(RepExp_Stats[["MSR"]]) / 2
  MSRn <- 10^((2 * s)/(sqrt(n)))

  MSRnTbl <- tibble(n, MSRn) %>%
    pivot_wider(names_from = n, names_prefix = 'n = ', values_from = MSRn)

  MeanRatioPlot <- mrplot(RepExp_Data, RepExp_Stats)

  R1R2CorrelationPlot <- r1r2plot(RepExp_Data, RepExp_Stats) +
    scale_x_continuous(trans='log10') +
    scale_y_continuous(trans='log10')

  list(Data = RepExp_Data, Stats = RepExp_Stats, MSRnTbl = MSRnTbl, MDPlot = MeanRatioPlot, CorrPlot = R1R2CorrelationPlot)
}

# Write report files ----------------------
repexp.save <- function(report, path) {

  # Top-level all reports directory
  all_reports_dir = file.path("Reports")

  # If the directory does not exist, then create directory
  if (!dir.exists(all_reports_dir)) {
    dir.create(all_reports_dir)
  }

  # Individual report directory
  report_dir <- file.path(all_reports_dir, path)

  # If the directory does not exist, then create directory
  if (!dir.exists(report_dir)) {
    dir.create(report_dir)
  }

  # Save the dataframe with calculated data
  write_csv(report[["Data"]], file = file.path(report_dir, "CalcData.csv"))

  # Save the summary statistics
  write_csv(report[["Stats"]], file = file.path(report_dir, "Stats.csv"))

  # Save the Bland-Altman Figure
  ggsave(filename = file.path(report_dir, "MDPlot.png"), plot = report[["MDPlot"]], height = 4, width = 6, units = 'in')

  # Save the Correlation Figure
  ggsave(filename = file.path(report_dir, "CorrPlot.png"), plot = report[["CorrPlot"]], height = 4, width = 4, units = 'in')
}

# Replicate-Experiment Example Analysis ------------------------------

# UsrData <- read_csv('Data/RepExpPotencyShift.csv')
# 
# MSR_3_Shift_Report <- repexp.potency(UsrData)
