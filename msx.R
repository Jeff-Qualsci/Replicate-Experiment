# Base functionality for MSX application for proposed Shiny Application

# Set up  server environment

library(tidyverse)

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


Pot <- TRUE # flag TRUE for potency data and FALSE for efficacy data

# Input Data and prepare for analysis

UsrData <- read_csv('Data/my_RepExpData.csv', col_names = TRUE, col_types = 'cnn') %>% 
  mutate(Potency = Pot,
         Meas1 = if_else(Potency, log10(Exp1), Exp1),
         Meas2 = if_else(Potency, log10(Exp2), Exp2), 
         MeasMean = (Meas1 + Meas2) / 2,
         MeasDiff = Meas1 - Meas2,
         GeoMean = 10 ^ MeasMean,
         Ratio = 10 ^ MeasDiff)

# Calculate MSX Statistics

MSxStats <- UsrData %>%
  summarise(n = n(), 
            t = qt(0.975, n - 1),
            MeanDiff = mean(MeasDiff),
            StdDev = sd(MeasDiff),
            MSx = 2 * StdDev,
            UDL = MeanDiff + (t * StdDev/sqrt(n)),
            LDL = MeanDiff - (t * StdDev/sqrt(n)),
            ULSA = MeanDiff + (t * StdDev),
            LLSA = MeanDiff- (t * StdDev))

ConCorr <- conc.corr(UsrData$Meas1, UsrData$Meas2)
ConCorrLabel <- paste("Concordance Correlation r =", format(ConCorr, digits = 3))

# Flag data pairs with potential outliers (MeasDiff outside of limits of agreement)

UsrData <- UsrData %>% 
  mutate(Outlier = MeasDiff > MSxStats$ULSA | MeasDiff < MSxStats$LLSA,
         Label = if_else(Outlier, Sample, NA))

# Generate MD Plot (Mean/Difference for MSD or Geometric Mean/Ratio for MSR) for report.

MDPlot <- if (Pot) {
  
  MSRLabel <- paste("MSR =", format(10 ^ MSxStats$MSx, digits = 2))
  
  ggplot(UsrData, aes(x =  GeoMean, y = Ratio)) +
    geom_point(shape = UsrData$Outlier, size = 3) +
    geom_text(aes(label = Label), na.rm = TRUE, hjust = -.5) +
    geom_hline(yintercept = 10 ^ MSxStats$MeanDiff, color = "mediumblue") +
    geom_hline(yintercept = 1, color = "black") +
    geom_hline(yintercept = 10 ^ MSxStats$UDL, color = "mediumblue", linetype = "dashed") +
    geom_text(x = 1, y = MSxStats$UDL, label = "Upper Ratio Limit", color = "mediumblue", vjust = -0.3) +
    geom_hline(yintercept = 10 ^MSxStats$LDL, color = "mediumblue", linetype = "dashed") +
    geom_text(x = 1, y = MSxStats$LDL, label = "Lower Ratio Limit", color = "mediumblue", vjust = 1) +
    geom_hline(yintercept = 10^ MSxStats$ULSA, color = "#D55E00", linetype = "dotdash") +
    geom_text(x = 1, y = MSxStats$ULSA, label = "Upper Limit of Aggreement", color = "#D55E00", vjust = -0.3) +
    geom_hline(yintercept = 10 ^ MSxStats$LLSA, color = "#D55E00", linetype = "dotdash") + 
    geom_text(x = 1, y = MSxStats$LLSA, label = "Lower Limit of Aggreement", color = "#D55E00", vjust = 1) +
    theme_minimal() +
    scale_x_continuous(trans='log10') + 
    scale_y_continuous(trans='log10') +
    labs(title="Potency Ratio vs Geometric Mean Potency",
         subtitle = MSRLabel,
         x ="Geometric Mean Potency", y = "Ratio")
  
} else {
  MSDLabel <- paste("MSD =", format(MSxStats$MSx, digits = 2))
  
  ggplot(UsrData, aes(x = MeasMean, y = MeasDiff)) +
    geom_point(shape = UsrData$Outlier, size = 3) +
    geom_text(aes(label = Label), na.rm = TRUE, hjust = -.5) +
    geom_hline(yintercept = MSxStats$MeanDiff, color = "mediumblue") +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = MSxStats$UDL, color = "mediumblue", linetype = "dashed") +
    geom_text(x = 1, y = MSxStats$UDL, label = "Upper Difference Limit", color = "mediumblue", vjust = -0.3) +
    geom_hline(yintercept = MSxStats$LDL, color = "mediumblue", linetype = "dashed") +
    geom_text(x = 1, y = MSxStats$LDL, label = "Lower Difference Limit", color = "mediumblue", vjust = 1) +
    geom_hline(yintercept = MSxStats$ULSA, color = "#D55E00", linetype = "dotdash") +
    geom_text(x = 1, y = MSxStats$ULSA, label = "Upper Limit of Aggreement", color = "#D55E00", vjust = -0.3) +
    geom_hline(yintercept = MSxStats$LLSA, color = "#D55E00", linetype = "dotdash") + 
    geom_text(x = 1, y = MSxStats$LLSA, label = "Lower Limit of Aggreement", color = "#D55E00", vjust = 1) +
    theme_minimal() +
    labs(title="Difference vs Mean Efficacy",
         subtitle = MSDLabel,
         x ="Mean Efficacy", y = "Difference")
  
}

MDPlot

ggsave(filename = 'Output/MDPlot.png', plot = MDPlot, height = 4, width = 6, units = "in")

# R1R2 Plot - Run1/Run2 with correlation.


R1R2Plot <- ggplot(UsrData,aes(x = Meas1, y = Meas2)) +
  geom_point() +
  geom_text(aes(label = Label), na.rm = TRUE, hjust = -.5) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", linewidth = 2) +
  geom_abline(slope = 1) +
  labs(title = "Correlation Run 1 vs Run 2",
       subtitle = ConCorrLabel,
       x = "Run 1",
       y = "Run2") +
  theme_minimal()

R1R2Plot

ggsave(filename = 'Output/CorrelationPlot.png', plot = R1R2Plot, height = 4, width = 4, units = "in")

# Generate output data files for user

  RepExpData <- if(Pot) { UsrData %>% 
      select(-Potency, -Meas1, -Meas2, -MeasDiff, -MeasMean, -Label) %>% 
      rename(`Geometric Mean` = GeoMean, `Ratio (Run 1/Run 2` = Ratio)
    
} else { UsrData %>% 
      select(-Potency, -Meas1, -Meas2, -GeoMean, -Ratio) %>% 
      rename(Mean = MeasMean, Difference = MeasDiff, -Label)
  
}
  
 RepExpStats <- if(Pot) {
   MSxStats %>% 
     select(-t, -StdDev) %>% 
     mutate(MeanDiff = 10 ^ MeanDiff,
            MSx = 10 ^ MSx,
            UDL = 10 ^ UDL,
            LDL = 10 ^ LDL,
            ULSA = 10 ^ ULSA,
            LLSA = 10 ^ LLSA) %>% 
     rename(MSR = MSx,
            `Upper Ratio Limit` = UDL,
            `Lower Ratio Limit` = LDL,
            `Upper Agreement Limit` = ULSA,
             `Lower Agreement Limit` = LLSA)
   
   
 } else {
   MSxStats %>% 
     select(-t, -StdDev) %>%
     rename(MSD = MSx,
            `Upper Difference Limit` = UDL,
            `Lower Difference Limit` = LDL,
            `Upper Agreement Limit` = ULSA,
            `Lower Agreement Limit` = LLSA)
 }
 
 write_csv(RepExpData, file = 'Output/RepExpData.csv')
 write_csv(RepExpStats, file = 'Output/RepExpStats.csv')
 