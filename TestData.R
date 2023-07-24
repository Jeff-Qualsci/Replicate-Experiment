# Create test data for replicate experiment AGM webtool

library(tidyverse)
setSessionTimeLimit(12252023) #For reproducible examples 

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

# Test Data Sets
# Standard MSD with 320 samples MSD = 20

MSD20Data320 <- msd_data(SmplNum = 320, TstMSD = 20, Shift = 0)
write_csv(MSD20Data320, file = 'Data/MSD20Data320.csv')

# MSD with 320 samples MSD = 20 Shift = 20

MSD20Data320Shift20 <- msd_data(SmplNum = 320, TstMSD = 20, Shift = 20)
write_csv(MSD20Data320Shift20, file = 'Data/MSD20Shift20Data3820.csv')

# Standard MSD with 10000 samples MSD = 30

MSD30Data10000 <- msd_data(SmplNum = 10000, TstMSD = 30, Shift = 0)
write_csv(MSD30Data10000, file = 'Data/MSD30Data10000.csv')

# Screening MSD set with 95% inactive and 5% active

ScrnInactMSDdata <- msd_inact(SmplNum = 9500, TstMSD = 10, Shift = 0)
write_csv(ScrnInactMSDdata, file = 'Data/ScrnInact.csv')

ScrnActMSDdata <- msd_act(SmplNum = 500, TstMSD = 30, Shift = 0)
write_csv(ScrnActMSDdata, file = 'Data/ScrnAct.csv')

ScrnMSDdata <- bind_rows(ScrnInactMSDdata, ScrnActMSDdata)
write.csv(ScrnMSDdata, file = 'Data/ScrnMSD.csv')

# MSR 32 samples MSR = 3

MSR3data32 <- msr_data(SmplNum = 32, TstMSR = 3)
write_csv(MSR3data32, file = 'Data/MSR3data32.csv')

# MSR 32 samples MSR = 3 Shift = 0.3

MSR3data32shift <- msr_data(SmplNum = 32, TstMSR = 3, Shift = 0.3)
write_csv(MSR3data32shift, file = 'Data/MSR3data32shift.csv')
