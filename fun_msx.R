# Functionalized version of msx.R

library(tidyverse)

# Test data generation functions ---------------------------------------
set.seed(12252023) #For reproducible examples 

# Generate Potency test data for Replicate Experiment
msr_data <- function(SmplNum, TstMSR, Shift = 0) {
  
  TstSD <- log10(TstMSR) / (2 * sqrt(2))
  Sample <- sample(c(1000000 : 9999999), SmplNum)
  TrueMeas <- runif(SmplNum, -3, 3)
  ExpData <- map(TrueMeas, \(m) rnorm(2, m, TstSD))
  
  TstData <- tibble(Sample, ExpData) %>% 
    unnest_wider(col = ExpData, names_sep = '_') %>% 
    rename(Exp1 = ExpData_1, Exp2 = ExpData_2) %>% 
    mutate(Exp2 = Exp2 + Shift,
           across(-Sample, ~ 10 ^ .x))
}

# Generate Efficacy test data for Replicate Experiment
msd_data <- function(SmplNum, TstMSD, Shift = 0) {
  
  TstSD <- TstMSD /(2 * sqrt(2))
  Sample <- sample(c(1000000 : 9999999), SmplNum)
  TrueMeas <- runif(SmplNum, -20, 120)
  ExpData <- map(TrueMeas, \(m) rnorm(2, m, TstSD))
  
  TstData <- tibble(Sample, ExpData) %>% 
    unnest_wider(col = ExpData, names_sep = '_') %>% 
    rename(Exp1 = ExpData_1, Exp2 = ExpData_2) %>% 
    mutate(Exp2 = Exp2 + Shift)
}

# Generate inactive efficacy test data for Replicate Experiment
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

# Generate active efficacy test data for Replicate Experiment
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

# Concordance Correlation V. Devanarayan --------------------------------
conc.corr <- function(x, y, alpha = 0.05)
{
  qval <- qnorm(1. - alpha/2.)
  N <- length(x)
  lm.xy <- lm(y ~ x)
  Intercept <- c(lm.xy$coeff[1.])
  Slope <- c(lm.xy$coeff[2.])
  A <- 2. * Slope * var(x)
  B <- var(y) + var(x)
  C <- (Intercept + (Slope - 1.) * (mean(x)))^2.
  rho <- A/(B + C)
  U <- (mean(x) - mean(y))^2./(sqrt(var(x)) * sqrt(var(y)))
  # cat("Zhat = ",Zhat,"\n")
  # cat("D=",D,"\n")
  # cat("E=",E,"\n")
  Zhat <- 0.5 * logb((1. + rho)/(1. - rho))
  D <- (4. * rho^3. * (1. - rho) * U)/(rho * (1. - rho^2.)^2.)
  E <- (2. * rho^4. * U^2.)/(rho^2. * (1. - rho^2.)^2.)
  # cat("var.Zhat=",var.Zhat,"\n")
  var.Zhat <- 1./(N - 2.) * (D - E)
  lowlim <- Zhat - qval * sqrt(var.Zhat)
  # cat(lowlim,upplim,"\n")
  # Now transform back:
  upplim <- Zhat + qval * sqrt(var.Zhat)
  lowlim <- (exp(2. * lowlim) - 1.)/(1. + exp(2. * lowlim))
  upplim <- (exp(2. * upplim) - 1.)/(1. + exp(2. * upplim))
  lowlim <- max(lowlim, -1.)
  upplim <- min(upplim, 1.)
  list(rho = rho, lowlim = lowlim, upplim = upplim)
}


# Replicate-Experiment Analysis Functions -------------------------------

# Summary Stats
repexp.stats <- function(df) {
  ConcCorr <- conc.corr(df$Exp1, df$Exp2)
  rho <- as.double(ConcCorr[1])  
  
  summarise(df,
            n = n(), 
            t = qt(0.975, n - 1),
            MeanDiff = mean(Difference),
            StdDev = sd(Difference),
            MSD = 2 * StdDev,
            UDL = MeanDiff + (t * StdDev/sqrt(n)),
            LDL = MeanDiff - (t * StdDev/sqrt(n)),
            ULSA = MeanDiff + (t * StdDev),
            LLSA = MeanDiff- (t * StdDev)) %>% 
    mutate(r = rho) %>% 
    select(-t, -StdDev)
}

#Mean Difference Plot
mdplot <- function(Data, Stats) {
  MSDLabel <- paste("MSD =", format(Stats$MSD, digits = 2))
  
  ggplot(Data, aes(x = Mean, y = Difference)) +
    geom_point(shape = Data$Outlier, size = 3) +
    geom_text(aes(label = Label), na.rm = TRUE, hjust = -.5) +
    geom_hline(yintercept = Stats$MeanDiff, color = "mediumblue") +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = Stats$UDL, color = "mediumblue", linetype = "dashed") +
    geom_text(x = 1, y = Stats$UDL, label = "Upper Difference Limit", color = "mediumblue", vjust = -0.3) +
    geom_hline(yintercept = Stats$LDL, color = "mediumblue", linetype = "dashed") +
    geom_text(x = 1, y = Stats$LDL, label = "Lower Difference Limit", color = "mediumblue", vjust = 1) +
    geom_hline(yintercept = Stats$ULSA, color = "#D55E00", linetype = "dotdash") +
    geom_text(x = 1, y = Stats$ULSA, label = "Upper Limit of Aggreement", color = "#D55E00", vjust = -0.3) +
    geom_hline(yintercept = Stats$LLSA, color = "#D55E00", linetype = "dotdash") + 
    geom_text(x = 1, y = Stats$LLSA, label = "Lower Limit of Aggreement", color = "#D55E00", vjust = 1) +
    theme_minimal() +
    labs(title="Difference vs Mean Efficacy",
         subtitle = MSDLabel,
         x ="Mean Efficacy", y = "Difference")
}

# Mean Ratio Plot
mrplot <- function(Data, Stats) {
  MSRLabel <- paste("MSR =", Stats$MSR)
  
  ggplot(Data, aes(x =  GeometricMean, y = Ratio)) +
    geom_point(shape = Data$Outlier, size = 3) +
    geom_text(aes(label = Label), na.rm = TRUE, hjust = -.5) +
    geom_hline(yintercept = Stats$MeanRatio, color = "mediumblue") +
    geom_hline(yintercept = 1, color = "black") +
    geom_hline(yintercept = Stats$UDL, color = "mediumblue", linetype = "dashed") +
    geom_text(x = 1, y = Stats$UDL, label = "Upper Ratio Limit", color = "mediumblue", vjust = -0.3) +
    geom_hline(yintercept = Stats$LDL, color = "mediumblue", linetype = "dashed") +
    geom_text(x = 1, y = Stats$LDL, label = "Lower Ratio Limit", color = "mediumblue", vjust = 1) +
    geom_hline(yintercept = Stats$ULSA, color = "#D55E00", linetype = "dotdash") +
    geom_text(x = 1, y = Stats$ULSA, label = "Upper Limit of Aggreement", color = "#D55E00", vjust = -0.3) +
    geom_hline(yintercept = Stats$LLSA, color = "#D55E00", linetype = "dotdash") + 
    geom_text(x = 1, y = Stats$LLSA, label = "Lower Limit of Aggreement", color = "#D55E00", vjust = 1) +
    theme_minimal() +
    scale_x_continuous(trans='log10') + 
    scale_y_continuous(trans='log10') +
    labs(title="Potency Ratio vs Geometric Mean Potency",
         subtitle = MSRLabel,
         x ="Geometric Mean Potency", y = "Ratio")
}

# R1R2 Plot - Run1/Run2 with correlation.

r1r2plot <- function(Data, Stats) {
  ggplot(Data,aes(x = Exp1, y = Exp2)) +
    geom_point() +
    geom_text(aes(label = Label), na.rm = TRUE, hjust = -.5) +
    geom_smooth(method = "lm", se = TRUE, linetype = "dashed", linewidth = 2) +
    geom_abline(slope = 1) +
    labs(title = "Correlation Run 1 vs Run 2",
         subtitle = paste("Concordance Correlation r =", Stats$r),
         x = "Run 1",
         y = "Run2") +
    theme_minimal()
  
}

#Replicate-Experiment Efficacy --------------------------------------------

repexp.efficacy<- function(df) {
  
RepExp_Data <- df %>% 
    mutate(Mean = (Exp1 + Exp2) / 2,
           Difference = Exp1 - Exp2)
 
RepExp_Stats <- repexp.stats(RepExp_Data) %>% 
  mutate(across(-n, \(x) signif(x, digits = 3)))

# Flag data pairs with potential outliers (MeasDiff outside of limits of agreement)

RepExp_Data <- RepExp_Data %>% 
  mutate(Outlier = Difference > RepExp_Stats$ULSA | Difference < RepExp_Stats$LLSA,
         Label = if_else(Outlier, Sample, NA))

MeanDifferencePlot <- mdplot(RepExp_Data, RepExp_Stats)

R1R2CorrelationPlot <- r1r2plot(RepExp_Data,RepExp_Stats)

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
    mutate(Outlier = Difference > RepExp_Stats$ULSA | Difference < RepExp_Stats$LLSA,
           Label = if_else(Outlier, Sample, NA),
           across(-c(1, 6, 7), ~10 ^ .x),
           across(-c(1, 6, 7),  \(x) signif(x, digits = 3))) %>% 
    rename(GeometricMean = Mean, Ratio = Difference)
  
  RepExp_Stats <- RepExp_Stats %>% 
    rename(MSR = MSD,
           MeanRatio = MeanDiff) %>% 
    mutate(across(-c(1, 8), ~10 ^ .x),
           across(-c(1), \(x) signif(x, digits = 3)))
  
  MeanRatioPlot <- mrplot(RepExp_Data, RepExp_Stats)
  
  R1R2CorrelationPlot <- r1r2plot(RepExp_Data,RepExp_Stats) +
    scale_x_continuous(trans='log10') + 
    scale_y_continuous(trans='log10')
  
  list(Data = RepExp_Data, Stats = RepExp_Stats, MRPlot = MeanRatioPlot, CorrPlot = R1R2CorrelationPlot)
}
  

# Replicate-Experiment Example Analysis ------------------------------

UsrData <- msd_data(320, 20)

MSD_320_20 <- repexp.efficacy(UsrData)

MSD_320_20[[1]]
MSD_320_20[[2]]
MSD_320_20[[3]]
MSD_320_20[[4]]

UsrData <- msr_data(40, 3.5)

MSR_40_3.5 <- repexp.potency(UsrData)
MSR_40_3.5[[1]]
MSR_40_3.5[[2]]
print(MSR_40_3.5[[3]])
MSR_40_3.5[[4]]