###################################################################################
#analysis of Nitrogen use (concentration and isotopic data) across historical Suillus and plant host specimens 
#from the Bell Museum in MN
#all associated code can be found on in the git repo https://github.com/MycoPunk/Suillus_N_use
#this script was run in R version 4.4.1
#data collection and preliminary data analysis preformed by T. Michaud
#code and final analysis written and preformed by L. Lofgren
#last updated 10.May.2025
###################################################################################


##### INITIATION AND DATA CLEANING #####

#set packages 
library(data.table)
library(dplyr)
library(tidyverse)
library(ggtree)
library(ape)
library(phytools)
library(ggplot2)
library(ggpubr)
library(devtools)
library(gridExtra)
library(lmerTest)
library(scales)
library(patchwork)
library(performance)
library(emmeans)
library(svglite)
sessionInfo()


#set seed for reproducibility
set.seed(666)

#set wd
setwd("~/Desktop/Project_N_use/")

#read in the data
data<- as.data.frame(fread("data/N_use_isotopes.csv", header = TRUE, sep = ",")) #note, this data sheet supplied by Talia includes data collected from WA, which is not being included in this analysis. 

#set color scheme 
fungus_colors <- c(
  "S. americanus" = "#5E3023",  # Rich brown
  "S. ampliporus" = "#D4AC6E",  # Warm sand/wheat
  "S. clintonianus" = "#355E2B",  # Richer forest green
  "S. luteus" = "#9B9B7A",  # Sage/olive
  "S. spraguei" = "#B07156",  # Terracotta/clay
  "S. weaverae" = "#6D5A72"   # Muted purple/aubergine
)


host_colors <- c(
  "L. laricina" = "#3A5A7E",  # Deeper blue/indigo
  "P. resinosa"  = "#954535",  # Burnt sienna
  "P. strobus" = "#D17A22"    # dusty rose
)

##clean the data by removing species other than the 6 fungal species and 3 host species used in this analysis 
#set fungal species to keep
keep_fungi<- c("Suillus weaverae", "Suillus spraguei", "Suillus americanus", "Suillus clintonianus", "Suillus ampliporus", "Suillus luteus")

#retain only the 6 fungal species used in this analysis 
data_filtered_fungi <- data %>%
  filter(species %in% keep_fungi)

#check that it worked
unique(data_filtered_fungi$species) #looks good
dim(data_filtered_fungi) #112 fungal samples in total
table(data_filtered_fungi$species)

#set host species to keep
keep_host<- c("Pinus strobus", "Pinus resinosa", "Larix laricina")
#retain only the 3 host species used in the analysis
data_filtered_host <- data %>%
  filter(species %in% keep_host)

#check that it worked
unique(data_filtered_host$species) #looks good
dim(data_filtered_host) #131 in total
table(data_filtered_host$species)

#due to isotope fractionation and differential distribution of protein in different tissues, the expectation is that cap will have higher N15 than stipe. 
#samples that do not adhere to this expectation may indicate contamination, abnormal physiology, or issues with preservation.
#here we use d15Ndif (which was pre-calculated by Talia by subtracting each d15N_stipe value from d15N_cap value), as well as d13Cdif (obtained by subtracting d13C_stipe from d13C_cap)
#T0 implement QC: If the difference is positive, the sample is retained. Otherwise, it is excluded 
data_filtered_fungi_clean <- subset(data_filtered_fungi, sym == "myco" & d15Ndif > 0 & d13Cdif > 0);dim(data_filtered_fungi_clean)

#quick look at the data that was excluded
data_excluded <- subset(data_filtered_fungi, !(sym == "myco" & d15Ndif > 0 & d13Cdif > 0 | sym == "host"))
dim(data_excluded) #we eliminated 5 samples, reducing the total from 112 to 107 fungal samples
table(data_excluded$species) #no species signal in the samples removed: we excluded 1 S. americanus, 2 S. clintonianus, and 2 S. weaverae

#how many of each species are in our final data set? 
table(data_filtered_fungi_clean$species)
#Suillus americanus   Suillus ampliporus Suillus clintonianus       Suillus luteus     Suillus spraguei     Suillus weaverae 
#20                   11                   24                    9                   15                   28 



##fix species names to match how they should appear in the manuscript
#create df for name map for fungi
names_df<- data.frame(cbind(raw_name = c("Suillus clintonianus", "Suillus spraguei", "Suillus americanus", "Suillus weaverae", "Suillus luteus", "Suillus ampliporus"), 
                            manuscript_name = c("S. clintonianus", "S. spraguei", "S. americanus", "S. weaverae", "S. luteus", "S. ampliporus")))
#turn the df into a named vector
rename_vector <- setNames(names_df$manuscript_name, names_df$raw_name)
#rename the names column 
data_filtered_fungi_clean <- data_filtered_fungi_clean %>%
  mutate(species = dplyr::recode(species, !!!as.list(rename_vector)))
#check that is worked
data_filtered_fungi_clean$species #looks good. 

#create df for name map for hosts
names_df2<- data.frame(cbind(raw_name = c("Larix laricina", "Pinus resinosa", "Pinus strobus"), 
                            manuscript_name = c("L. laricina", "P. resinosa", "P. strobus")))
#turn the df into a named vector
rename_vector2 <- setNames(names_df2$manuscript_name, names_df2$raw_name)
#rename the names column 
data_filtered_host <- data_filtered_host %>%
  mutate(species = dplyr::recode(species, !!!as.list(rename_vector2)))
#check
data_filtered_host$species #looks good. 


##over what years were the samples collected?
#for fungi
unique(sort(data_filtered_fungi_clean$year)) #fungi were collected between  1961 and 2019
#for hosts
unique(sort(data_filtered_host$year)) #hosts were collected between  1961 and 2021





##### MODELING #####

#we will treat year, species, and interactions between year and species as fixed effects.
#and county and month as random effects. 
#however, because there is some missing data for month, we have to remove that fist or the model won't run.
#we will add these values back into the model if they are not significant after model finding

#identify  missing values
#for fungi
missing_summary_myco <- data.frame(
  Column = names(data_filtered_fungi_clean),
  Missing_Values = sapply(data_filtered_fungi_clean, function(x) sum(is.na(x)))
)
missing_summary_myco #there are 4 fungal samples with no month, and none missing county, year or species
#remove the rows with missing month values
data_filtered_fungi_clean_minus_missing_month <- data_filtered_fungi_clean[!is.na(data_filtered_fungi_clean$month), ]
#check that it worked
dim(data_filtered_fungi_clean)
dim(data_filtered_fungi_clean_minus_missing_month) #looks good

#for hosts
missing_summary_host <- data.frame(
  Column = names(data_filtered_host),
  Missing_Values = sapply(data_filtered_host, function(x) sum(is.na(x)))
)
missing_summary_host #no missing data for month, county, year, or species



## Model total N ##

#for Suillus, model total N over time. Here "n" is the average total n between cap and stipe, pre-calculated by Talia. 
#initiate the full model for selection using the data set with missing data (month) removed
N_MycoInit <- lmer(n ~ year + species + year:species + (1|county)+ (1|month), 
                     data = data_filtered_fungi_clean_minus_missing_month);summary(N_MycoInit)
step(N_MycoInit) #model selection using AIC = starts with full model and sequentially removes terms that don't improve the model fit according to AIC
#Backward reduced fixed-effect table:
#  Eliminated Df Sum of Sq    RSS     AIC F value   Pr(>F)   
#year:species          1  5    1.1923 22.712 -141.72  1.0084 0.417446   
#year                  0  1    2.5183 25.230 -132.89 10.6445 0.001530 **
#  species               0  5    5.1415 27.854 -130.70  4.3464 0.001311 **
#Model found:
#  n ~ year + species

#the above indicates that the simpler model with just year and species is a better model

#run the model identified above, using the full data set since month was not significant and not included in the model
N_MycoFin<- lm(n ~ year + species , data = data_filtered_fungi_clean);summary(N_MycoFin)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.95045 -0.31627 -0.07775  0.26083  1.48803 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)   
#(Intercept)            -15.782250   5.647316  -2.795  0.00623 **
#  year                     0.009364   0.002845   3.291  0.00138 **
#  speciesS. ampliporus     0.058962   0.187441   0.315  0.75375   
#speciesS. clintonianus  -0.193921   0.150562  -1.288  0.20073   
#speciesS. luteus        -0.317045   0.201016  -1.577  0.11791   
#speciesS. spraguei       0.468384   0.169299   2.767  0.00675 **
#  speciesS. weaverae       0.080874   0.144978   0.558  0.57820   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.4951 on 100 degrees of freedom
#Multiple R-squared:  0.2311,	Adjusted R-squared:  0.1849 
#F-statistic: 5.008 on 6 and 100 DF,  p-value: 0.0001588

#conclusion = Suillus fruitbody N is significantly increasing overtime


### How much N accumulation is that per year? ###
#get all unique species in the data set
unique_species <- unique(data_filtered_fungi_clean$species)

#create prediction dfs for each species at both years
predictions_1961 <- list()
predictions_2019 <- list()

#to get the CIs, using the predict() with interval="confidence"
for (i in 1:length(unique_species)) {
  current_species <- unique_species[i]
  new_data_1961 <- data.frame(year = 1961, species = as.character(current_species))
  new_data_2019 <- data.frame(year = 2019, species = as.character(current_species))
  
  # Get predictions with confidence intervals
  predictions_1961[[i]] <- predict(N_MycoFin, newdata = new_data_1961, interval = "confidence")
  predictions_2019[[i]] <- predict(N_MycoFin, newdata = new_data_2019, interval = "confidence")
}

#calculate means across species for each component (fit, lower, upper).
mean_1961_fit <- mean(sapply(predictions_1961, function(x) x[,"fit"]))
mean_1961_lwr <- mean(sapply(predictions_1961, function(x) x[,"lwr"]))
mean_1961_upr <- mean(sapply(predictions_1961, function(x) x[,"upr"]))

mean_2019_fit <- mean(sapply(predictions_2019, function(x) x[,"fit"]))
mean_2019_lwr <- mean(sapply(predictions_2019, function(x) x[,"lwr"]))
mean_2019_upr <- mean(sapply(predictions_2019, function(x) x[,"upr"]))

#calculate absolute and percent change
absolute_change <- mean_2019_fit - mean_1961_fit
percent_change <- (absolute_change / mean_1961_fit) * 100

#calculate percent change CIs based on the ratio of the extremes
percent_change_lower <- ((mean_2019_lwr - mean_1961_upr) / mean_1961_upr) * 100
percent_change_upper <- ((mean_2019_upr - mean_1961_lwr) / mean_1961_lwr) * 100

#print results
cat(sprintf("1961 estimated mean N concentration: %.2f%% (95%% CI: %.2f%%, %.2f%%)\n", 
            mean_1961_fit, mean_1961_lwr, mean_1961_upr))
cat(sprintf("2019 estimated mean N concentration: %.2f%% (95%% CI: %.2f%%, %.2f%%)\n", 
            mean_2019_fit, mean_2019_lwr, mean_2019_upr))
cat(sprintf("Absolute change in N concentration: %.2f%% (95%% CI: %.2f%%, %.2f%%)\n", 
            absolute_change, mean_2019_lwr - mean_1961_upr, mean_2019_upr - mean_1961_lwr))
cat(sprintf("Percent change from 1961 to 2019: %.1f%% (95%% CI: %.1f%%, %.1f%%)\n", 
            percent_change, percent_change_lower, percent_change_upper))




#For host, model total N
N_HostInit <- lmer(n ~ year + species + year:species + (1|county)+ (1|month), 
                     data = data_filtered_host);summary(N_HostInit)
step(N_HostInit)


#Backward reduced random-effect table:
#  Eliminated npar  logLik    AIC      LRT Df Pr(>Chisq)
#<none>                     9 -87.154 192.31                       
#(1 | county)          1    8 -87.154 190.31 0.000000  1     1.0000
#(1 | month)           2    7 -87.164 188.33 0.020863  1     0.8852

#Backward reduced fixed-effect table:
#  Eliminated Df Sum of Sq    RSS     AIC F value Pr(>F)    
#year:species          1  2    0.1152 22.070 -225.31  0.3281 0.7209    
#year                  2  1    0.0047 22.074 -227.28  0.0268 0.8703    
#species               0  2   18.5839 40.658 -151.27 53.8800 <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Model found:
#  n ~ species

#the above indicates that the model with n concentration and species is a better model, eliminating time.

#run the model identified above.
#note, were adding 0 to the model to get the exact values for all comparisons. 
N_HostFin<- lm(n ~ 0 + species , data = data_filtered_host);summary(N_HostFin)


#Residuals:
#  Min       1Q   Median       3Q      Max 
#-1.04309 -0.24476 -0.02676  0.17691  1.60691 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#speciesL. laricina  2.09309    0.05600   37.38   <2e-16 ***
#  speciesP. resinosa  1.15676    0.07122   16.24   <2e-16 ***
#  speciesP. strobus   1.65976    0.06408   25.90   <2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.4153 on 128 degrees of freedom
#Multiple R-squared:  0.948,	Adjusted R-squared:  0.9467 
#F-statistic: 777.3 on 3 and 128 DF,  p-value: < 2.2e-16


#becasue we're using ln here instead of lmer, re-run without forcing  0, or it will mess up your R2
N_HostFin<- lm(n ~ species , data = data_filtered_host);summary(N_HostFin)
#R2 0.4486
#Residual standard error: 0.4153 on 128 degrees of freedom
#Multiple R-squared:  0.4571,	Adjusted R-squared:  0.4486 
#F-statistic: 53.88 on 2 and 128 DF,  p-value: < 2.2e-16



### make figure panel of total N by year for host and fungi ###
#get range to make scales equal 
y_limits <- range(c(data_filtered_host$n, data_filtered_fungi_clean$n), na.rm = TRUE)


#plot total host N
p_host <- ggplot(data_filtered_host, aes(x = year, y = n)) +
  geom_point(aes(color = species), pch=16, size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, aes(color = species), linewidth = 1) +
  scale_color_manual(
    values = host_colors,
    labels = c(
      "L. laricina" = expression(italic("L. laricina")),
      "P. resinosa" = expression(italic("P. resinosa")),
      "P. strobus" = expression(italic("P. strobus"))
    )
  ) +
  scale_y_continuous(limits = y_limits) +
  labs(y = "[N]", x = "Year", color = "Host") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )



#plot total fungus N
#since the overall trend is significant, extract statistical information from the model to add to the graph.
model_stats <- summary(N_MycoFin)
adj_r2 <- round(model_stats$adj.r.squared, 2)
overall_r2 <- round(model_stats$r.squared, 2)
f_stat <- model_stats$fstatistic
overall_p <- pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
if (as.numeric(overall_p) < 0.001) overall_p <- "p < 0.001"

p_fungus <- ggplot(data_filtered_fungi_clean, aes(x = year, y = n)) +
  #add species-specific points and trend lines
  geom_point(aes(color = species), pch=16, size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, aes(color = species), linewidth = 1) +
  #add overall trend line with confidence intervals
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray50", linetype = "solid", linewidth = 1.2) +
  scale_color_manual(
    values = fungus_colors,
    breaks = c("S. americanus", "S. ampliporus", "S. clintonianus", "S. luteus", "S. spraguei", "S. weaverae"),
    labels = c(
      "S. americanus" = expression(italic("S. americanus")),
      "S. ampliporus" = expression(italic("S. ampliporus")),
      "S. clintonianus" = expression(italic("S. clintonianus")),
      "S. luteus" = expression(italic("S. luteus")),
      "S. spraguei" = expression(italic("S. spraguei")),
      "S. weaverae" = expression(italic("S. weaverae"))
    )
  ) +
  scale_y_continuous(limits = y_limits) +
  labs(y = "[N]", x = "Year", color = "Fungus") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  #add statistical annotation to the lower right corner
  annotate("text", 
           x = max(data_filtered_fungi_clean$year, na.rm = TRUE) - 2, 
           y = min(data_filtered_fungi_clean$n, na.rm = TRUE) + 0.1,
           label = paste0("R² = ", round(overall_r2, 2),
                          "\nR² (adjusted) = ", round(adj_r2, 2), 
                          "\np = 0.003"),
           hjust = 1, vjust = 1, size = 3.5)

p_fungus


annotate("text", 
         x = max(data_filtered_host$year, na.rm = TRUE) - 5, 
         y = min(data_filtered_host$d13C, na.rm = TRUE) + 0.3,
         label = paste0("R² (fixed) = ", round(r2_values_host_C$R2_marginal, 2), 
                        "\nR² (fixed+random) = ", round(r2_values_host_C$R2_conditional, 2),
                        "\np = 0.003"),
         hjust = 0, 
         size = 3.5)



## Model δ15N ##

#for fungus, model δ15N over time. Here "d15N" is the average total δ15N between cap and stipe.
N15_MycoInit <- lmer(d15N ~ year + species + year:species + (1|county)+ (1|month), 
                   data = data_filtered_fungi_clean_minus_missing_month);summary(N_MycoInit)
step(N15_MycoInit) 

#Backward reduced random-effect table:
  
#  Eliminated npar  logLik    AIC    LRT Df Pr(>Chisq)  
#<none>                    15 -247.71 525.42                       
#(1 | month)           1   14 -248.12 524.24 0.8161  1    0.36632  
#(1 | county)          0   13 -249.84 525.68 3.4475  1    0.06335 .
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Backward reduced fixed-effect table:
#  Degrees of freedom method: Satterthwaite 

#Eliminated Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#year:species          1  47.39   9.479     5 84.419  1.7209    0.1386    
#year                  2   1.77   1.767     1 87.444  0.2944    0.5888    
#species               0 631.90 126.380     5 94.243 21.0983 4.248e-14 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Model found:
#  d15N ~ species + (1 | county)

#best model preserves county as a random effect 

#run on final model - note we're using lmer here rather than lm because it has a random effect (county), which lm can't handle.
#we're running this on the full data set because month is not identified as a significant random effect. 
N15_MycoFin <- lmer(d15N ~ 0 + species  + (1|county), data = data_filtered_fungi_clean);summary(N15_MycoFin)


# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: d15N ~ 0 + species + (1 | county)
#    Data: data_filtered_fungi_clean
# 
# REML criterion at convergence: 493.8
# 
# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -2.61351 -0.67942  0.00788  0.65966  2.33795 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  county   (Intercept) 1.353    1.163   
#  Residual             5.810    2.410   
# Number of obs: 107, groups:  county, 22
# 
# Fixed effects:
#                        Estimate Std. Error      df t value Pr(>|t|)    
# speciesS. americanus     8.1105     0.6566 47.0159  12.353 2.28e-16 ***
# speciesS. ampliporus     3.7076     0.8034 83.4903   4.615 1.41e-05 ***
# speciesS. clintonianus   6.6561     0.6219 37.5442  10.704 5.80e-13 ***
# speciesS. luteus         8.1372     0.8859 88.7260   9.186 1.59e-14 ***
# speciesS. spraguei      10.1229     0.7208 63.1707  14.044  < 2e-16 ***
# speciesS. weaverae      12.2622     0.5633 31.8647  21.768  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             spcsS.amr spcsS.amp spcS.c spcS.l spcS.s
# spcsS.amplp 0.110                                   
# spcsS.clntn 0.195     0.151                         
# specisS.lts 0.147     0.085     0.139               
# specsS.sprg 0.220     0.094     0.208  0.100        
# specisS.wvr 0.283     0.145     0.206  0.130  0.220 


#get marginal and conditional R2 values using the performance package
#rerun model without the 0 intercept (we don't need to do this here because we're uisng lmer and not lm, but a good policy anyway to reset things)
N15_MycoFin <- lmer(d15N ~ species  + (1|county), data = data_filtered_fungi_clean);summary(N15_MycoFin)
r2_values <- r2(N15_MycoFin)
r2_values


#for host, model δ15N concentration over time. Here "d15N" is the total δ15N
N15_HostInit <- lmer(d15N ~ year + species + year:species + (1|county)+ (1|month), 
                     data = data_filtered_host);summary(N15_HostInit)

step(N15_HostInit) 
# Backward reduced random-effect table:
#   
#   Eliminated npar  logLik    AIC    LRT Df Pr(>Chisq)  
# <none>                     9 -318.56 655.13                       
# (1 | month)           1    8 -318.56 653.13 0.0000  1    0.99999  
# (1 | county)          0    7 -320.42 654.85 3.7228  1    0.05368 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Backward reduced fixed-effect table:
#   Degrees of freedom method: Satterthwaite 
# 
# Eliminated  Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)    
# year:species          1    1.18    0.59     2 115.34  0.0954 0.9091    
# year                  2   10.18   10.18     1 115.41  1.6714 0.1987    
# species               0 1140.89  570.45     2 121.82 92.0187 <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Model found:
#   d15N ~ species + (1 | county)

#again, county is retained. 

#run final model
N15_HostFin <- lmer(d15N ~ 0 + species  + (1|county), data = data_filtered_host);summary(N15_HostFin)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: d15N ~ 0 + species + (1 | county)
#    Data: data_filtered_host
# 
# REML criterion at convergence: 624.6
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -4.2422 -0.4440  0.0647  0.4568  2.3540 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  county   (Intercept) 1.157    1.075   
#  Residual             6.199    2.490   
# Number of obs: 131, groups:  county, 39
# 
# Fixed effects:
#                    Estimate Std. Error      df t value Pr(>|t|)    
# speciesL. laricina  -0.2383     0.4063 56.7125  -0.587     0.56    
# speciesP. resinosa   7.5665     0.4890 91.0625  15.473  < 2e-16 ***
# speciesP. strobus    3.7697     0.4489 76.8902   8.397 1.76e-12 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             spcL.l spcP.r
# specsP.rsns 0.175        
# spcsP.strbs 0.154  0.177 


#get marginal and conditional R2 values using the performance package
#first re-run model without the 0 intercept, again don't need to but good policy 
N15_HostFin <- lmer(d15N ~ species + (1|county), data = data_filtered_host);summary(N15_HostFin)

r2_values_host <- r2(N15_HostFin)
r2_values_host



### make figure panel of total δN15 by year for host and fungi ###
#get range to make scales equal 
y_limits_delta <- range(c(data_filtered_host$d15N, data_filtered_fungi_clean$d15N), na.rm = TRUE)


#host δN15 plot
p_host_delta_15 <- ggplot(data_filtered_host, aes(x = year, y = d15N)) +
  geom_point(aes(color = species), pch=16, size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, aes(color = species), linewidth = 1) +
  scale_color_manual(
    values = host_colors,
    labels = c(
      "L. laricina" = expression(italic("L. laricina")),
      "P. resinosa" = expression(italic("P. resinosa")),
      "P. strobus" = expression(italic("P. strobus"))
    )
  ) +
  scale_y_continuous(limits = y_limits_delta) +
  labs(y = "[δ15N]", x = "Year", color = "Host") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

#fungus δN15 plot
p_fungus_delta_15 <- ggplot(data_filtered_fungi_clean, aes(x = year, y = d15N)) +
  geom_point(aes(color = species), pch=16, size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, aes(color = species), linewidth = 1) +
  scale_color_manual(
    values = fungus_colors,
    breaks = c("S. spraguei", "S. americanus", "S. weaverae", "S. luteus", "S. clintonianus", "S. ampliporus"),
    labels = c(
      "S. americanus" = expression(italic("S. americanus")),
      "S. ampliporus" = expression(italic("S. ampliporus")),
      "S. clintonianus" = expression(italic("S. clintonianus")),
      "S. luteus" = expression(italic("S. luteus")),
      "S. spraguei" = expression(italic("S. spraguei")),
      "S. weaverae" = expression(italic("S. weaverae"))
    )
  ) +
  scale_y_continuous(limits = y_limits_delta) +
  labs(y = "[δ15N]", x = "Year", color = "Fungus") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )




### graph each species separately ###

#convert species to factor to reverse the order in the legend
data_filtered_fungi_clean$species <- factor(
  data_filtered_fungi_clean$species,
  levels = rev(c(species_order <-c("S. spraguei", "S. americanus", "S. weaverae", "S. luteus", "S. clintonianus", "S. ampliporus")
))
)


#set colors
d15N_boxplot <- ggplot(data_filtered_fungi_clean, 
                       aes(x = species, y = d15N, fill = species)) +
  #box plot layer
  geom_boxplot(alpha = 0.5, width = 0.5, outlier.shape = NA) +
  #data points layer
  geom_jitter(aes(color = species), width = 0.25, height = 0, 
              pch=16, alpha = 0.7, size = 2) +
  #apply colors
  scale_fill_manual(values = fungus_colors) +
  scale_color_manual(values = fungus_colors) +
  #flip coordinates
  coord_flip() +
  #set labels
  labs(
    x = "",  #remove x-axis label
    y = "δ15N"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none", 
    axis.text.y = element_text(face = "italic")
  )


d15N_boxplot

#make a second version wihtout species labels for saving

#convert species to factor to correrct order in the legend
data_filtered_fungi_clean$species <- factor(
  data_filtered_fungi_clean$species,
  levels = rev(c("S. spraguei", "S. americanus", "S. weaverae", "S. luteus", "S. clintonianus", "S. ampliporus"))
)


#set colors
d15N_boxplot_nolabs <- ggplot(data_filtered_fungi_clean, 
                       aes(x = species, y = d15N, fill = species)) +
  #box plot layer
  geom_boxplot(alpha = 0.5, width = 0.5, outlier.shape = NA) +
  #data points layer
  geom_jitter(aes(color = species), width = 0.25, height = 0, 
              pch=16, alpha = 0.7, size = 2) +
  #apply colors
  scale_fill_manual(values = fungus_colors) +
  scale_color_manual(values = fungus_colors) +
  #flip coordinates
  coord_flip() +
  #set labels
  labs(
    x = "",  #remove x-axis label
    y = "δ15N"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none", 
    axis.text.y = element_blank()
  )


d15N_boxplot_nolabs



## Model δ13C ##

#For Suillus, model δ13C concentration over time. Here "d13C" is the average total δ13C between cap and stipe.
C13_MycoInit <- lmer(d13C ~ year + species + year:species + (1|county)+ (1|month), 
                     data = data_filtered_fungi_clean_minus_missing_month);summary(C13_MycoInit)
step(C13_MycoInit) 


# boundary (singular) fit: see help('isSingular')
# Backward reduced random-effect table:
#   
#   Eliminated npar  logLik    AIC     LRT Df Pr(>Chisq)
# <none>                    15 -154.47 338.93                      
# (1 | county)          1   14 -154.47 336.93 0.00000  1     1.0000
# (1 | month)           2   13 -154.56 335.12 0.18677  1     0.6656
# 
# Backward reduced fixed-effect table:
#   Eliminated Df Sum of Sq    RSS     AIC F value Pr(>F)  
# year:species          0  5    11.823 88.531 -1.5912  2.8052 0.0211 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Model found:
#   d13C ~ year + species + year:species

#best described by, year, species and interaction between year and species (no month tho to use full data set)

#run final model, but do not set the 0 intercept separately, use emtrends instead to keep interpretation consistent with the rest of the paper
C13_MycoFin <- lm(d13C ~ year + species + year:species, 
                     data = data_filtered_fungi_clean);summary(C13_MycoFin)


#because we used lm instead of lmer as for the N15 values (because there are no mixed effects) to get the slopes and CI's you need to run an additional step
#get species-specific year slopes directly
slopes <- emtrends(C13_MycoFin, ~species, var="year")
summary(slopes, infer=TRUE)


## Make figure for delta13C for each species of Suillus ##

#convert species to a factor to set order
data_filtered_fungi_clean$species <- factor(
  data_filtered_fungi_clean$species,
  levels = c("S. spraguei", "S. americanus", "S. weaverae", "S. luteus", "S. clintonianus", "S. ampliporus")
)

#run separate linear models for each species to get R² values
r2_values <- data.frame(species = unique(data_filtered_fungi_clean$species),
                        r2 = NA)

for (sp in levels(data_filtered_fungi_clean$species)) {
  sp_data <- subset(data_filtered_fungi_clean, species == sp)
  sp_model <- lm(d13C ~ year, data = sp_data)
  r2_values$r2[r2_values$species == sp] <- summary(sp_model)$adj.r.squared
}

#create the basic plot with ordered facets
d13C_over_time_p <- ggplot(data_filtered_fungi_clean, aes(x = year, y = d13C, color = species)) +
  geom_point(pch=16, alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = species), alpha = 0.6) +
  facet_wrap(~ species, ncol = 3) +
  scale_color_manual(values = fungus_colors) +
  scale_fill_manual(values = fungus_colors) +
  labs(
    x = "Year",
    y = "δ13C"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = NA),
    strip.text = element_text(face = "italic"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

#get the significant species from the slopes summary
slopes_summary <- summary(slopes, infer=TRUE)

#create the annotation data frame
annotation_data <- data.frame(
  species = slopes_summary$species,
  p_value = slopes_summary$p.value,
  significant = slopes_summary$p.value < 0.05
)

#merge annotations with R² values
annotation_data <- merge(annotation_data, r2_values)

#create formatted text with both p-value and R²
annotation_data$label <- ifelse(
  annotation_data$significant,
  paste0("R² = ", round(annotation_data$r2, 2), "\n", 
         ifelse(annotation_data$p_value < 0.001, "p < 0.001", 
                paste0("p = ", round(annotation_data$p_value, 3)))),
  NA
)

#filter to keep only significant species
annotation_data <- annotation_data[annotation_data$significant, ]

#for each significant species, find a good position for the annotation
annotation_positions <- data.frame()
for (sp in annotation_data$species) {
  sp_data <- subset(data_filtered_fungi_clean, species == sp)
  # Calculate position in the lower right of each panel
  x_range <- range(sp_data$year, na.rm = TRUE)
  y_range <- range(sp_data$d13C, na.rm = TRUE)
  
  #position at 80% of x-range and 5% above bottom of y-range
  x_pos <- x_range[1] + 0.8 * diff(x_range)
  y_pos <- y_range[1] + 0.05 * diff(y_range)
  
  annotation_positions <- rbind(annotation_positions, 
                                data.frame(species = sp, x = x_pos, y = y_pos))
}

#merge positions with annotation text
annotation_data <- merge(annotation_data, annotation_positions)

#add the annotations using geom_text
d13C_over_time_p <- d13C_over_time_p + 
  geom_text(data = annotation_data, 
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            hjust = 1, 
            color = "black",
            size = 3,
            fontface = "plain")

#plot
d13C_over_time_p




#For Host, model δC concentration over time. Here "C13N" is the average total C13
C13_HostInit <- lmer(d13C ~ year + species + year:species + (1|county)+ (1|month), 
                     data = data_filtered_host);summary(C13_HostInit)

step(C13_HostInit) 

# Backward reduced random-effect table:
#   
#   Eliminated npar  logLik    AIC   LRT Df Pr(>Chisq)  
# <none>                     9 -214.97 447.93                      
# (1 | month)           1    8 -214.97 445.93 0.000  1    1.00000  
# (1 | county)          0    7 -216.44 446.88 2.955  1    0.08561 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Backward reduced fixed-effect table:
#   Degrees of freedom method: Satterthwaite 
# 
# Eliminated  Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# year:species          1  0.2943  0.1471     2 114.26  0.1239 0.883630   
# species               2  5.8079  2.9039     2 121.61  2.4747 0.088417 . 
# year                  0 11.0382 11.0382     1 117.20  9.3879 0.002711 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Model found:
#   d13C ~ year + (1 | county)

#run final model, - do not need to use 0 to set intercept because we're using lmer
C13_HostFin <- lmer(d13C ~ year + (1 | county), 
                  data = data_filtered_host);summary(C13_HostFin)
step(C13_HostFin) 


#get species-specific p values
#first add species back in
C13_HostFin_with_species <- lmer(d13C ~ year*species + (1 | county), 
                                 data = data_filtered_host)
host_slopes <- emtrends(C13_HostFin_with_species, ~species, var="year")
summary(host_slopes, infer=TRUE)


#to get marginal and conditional R2
r2_values_host_C <- r2(C13_HostFin)
r2_values_host_C



## Make figure for delta13C for each species of Suillus ##

#note, we don't plot each species separately as in Suillus, because 'species' was not significant in the model.
d13C_host_combined <- ggplot(data_filtered_host, aes(x = year, y = d13C, color = species)) +
  #points
  geom_point(pch=16, alpha = 0.6, size = 3) +
  #species-specific regression lines
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  #overall regression line
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray50", 
              linewidth = 1.5, aes(group = 1)) +
  #apply color palette
  scale_color_manual(
    values = host_colors,
    labels = c(
      "L. laricina" = expression(italic("L. laricina")),
      "P. resinosa" = expression(italic("P. resinosa")),
      "P. strobus" = expression(italic("P. strobus"))
    )
  ) +
  #add labels
  labs(
    x = "",
    y = "",
    color = ""
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "italic"),
    legend.position = "none", 
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +

annotate("text", 
         x = max(data_filtered_host$year, na.rm = TRUE) - 5, 
         y = min(data_filtered_host$d13C, na.rm = TRUE) + 0.3,
         label = paste0("R² (fixed) = ", round(r2_values_host_C$R2_marginal, 2), 
                        "\nR² (fixed+random) = ", round(r2_values_host_C$R2_conditional, 2),
                        "\np = 0.003"),
         hjust = 0, 
         size = 3.5)

d13C_host_combined




### How much δ13C depletion is that per year for host trees? ###

#get all unique tree species in the dataset
unique_species <- unique(data_filtered_host$species)

#get the fixed effects from the model
fixed_effects <- fixef(C13_HostFin)
intercept <- fixed_effects["(Intercept)"]
year_coef <- fixed_effects["year"]

#calculate predicted values for 1961 and 2019
pred_1961 <- intercept + year_coef * 1961
pred_2019 <- intercept + year_coef * 2019

#get standard errors for confidence intervals
model_summary <- summary(C13_HostFin)
year_se <- model_summary$coefficients["year", "Std. Error"]
intercept_se <- model_summary$coefficients["(Intercept)", "Std. Error"]

#calculate 95% CI
t_value <- qt(0.975, df = model_summary$coefficients["year", "df"])
margin_1961 <- t_value * sqrt(intercept_se^2 + (1961^2 * year_se^2))
margin_2019 <- t_value * sqrt(intercept_se^2 + (2019^2 * year_se^2))

#calculate values with CI
mean_1961_fit <- pred_1961
mean_1961_lwr <- pred_1961 - margin_1961
mean_1961_upr <- pred_1961 + margin_1961

mean_2019_fit <- pred_2019
mean_2019_lwr <- pred_2019 - margin_2019
mean_2019_upr <- pred_2019 + margin_2019

#get absolute change
absolute_change <- mean_2019_fit - mean_1961_fit
absolute_change_per_year <- year_coef # This is already the change per year from the model

#get CI for the difference
diff_se <- year_se * (2019 - 1961)
diff_lwr <- absolute_change - t_value * diff_se
diff_upr <- absolute_change + t_value * diff_se

#get percent change
percent_change <- (absolute_change / abs(mean_1961_fit)) * 100
percent_change_lower <- (diff_lwr / abs(mean_1961_fit)) * 100
percent_change_upper <- (diff_upr / abs(mean_1961_fit)) * 100

#in words:
cat(sprintf("1961 estimated mean δ13C: %.2f‰ (95%% CI: %.2f‰, %.2f‰)\n", 
            mean_1961_fit, mean_1961_lwr, mean_1961_upr))
cat(sprintf("2019 estimated mean δ13C: %.2f‰ (95%% CI: %.2f‰, %.2f‰)\n", 
            mean_2019_fit, mean_2019_lwr, mean_2019_upr))
cat(sprintf("Absolute change in δ13C: %.2f‰ (95%% CI: %.2f‰, %.2f‰)\n", 
            absolute_change, diff_lwr, diff_upr))
cat(sprintf("Change per year: %.4f‰ yr⁻¹\n", absolute_change_per_year))
cat(sprintf("Percent change from 1961 to 2019: %.1f%% (95%% CI: %.1f%%, %.1f%%)\n", 
            percent_change, percent_change_lower, percent_change_upper))




### model difference in deltaN15 between cap and stipe ###

#initial model with all potential fixed and random effects- here also including average N concentration (n), and the difference in n concentration between cap and stipe (xNdif)
N15_diff_init <- lmer(d15Ndif ~ year + species + n + xNdif + year:species + (1|county) + (1|month), 
                      data = data_filtered_fungi_clean_minus_missing_month);summary(N15_diff_init)


#run stepwise model selection to identify the optimal model structure
step(N15_diff_init)

#Model found:
#d15Ndif ~ year + species + n + xNdif + year:species

#based on the results of step(), fit your final model
# (assuming the step result preserves species as a fixed effect and county as a random effect)
N15_diff_final <- lm(d15Ndif ~ 0 + year + species + n + xNdif + year:species, 
                       data = data_filtered_fungi_clean)
summary(N15_diff_final)


#get R² values
#run without forcing 0 intercpt
N15_diff_final <- lm(d15Ndif ~ year + species + n + xNdif + year:species, 
                     data = data_filtered_fungi_clean)
summary(N15_diff_final)
r2_values_diff <- r2(N15_diff_final)
r2_values_diff


### make supp. plots of cap vs. stipe ###

#reformat to long data for the connected line plot
d15N_long <- data_filtered_fungi_clean %>%
  select(tubeLabel, species, year, d15N_stipe, d15N_cap) %>%
  pivot_longer(
    cols = c(d15N_stipe, d15N_cap),
    names_to = "tissue",
    values_to = "d15N"
  ) %>%
  mutate(
    tissue = factor(tissue, levels = c("d15N_stipe", "d15N_cap"), 
                    labels = c("Stipe", "Cap"))
  )

#plot cap vs. stipe delta N values
CS_Ndelta <- ggplot(d15N_long, aes(x = tissue, y = d15N, group = tubeLabel)) +
  #add lines connecting cap and stipe for each sample with transparency
  geom_line(aes(color = species), size = 0.7, alpha = 0.5) +
  #add points for each tissue with transparency
  geom_point(aes(color = species), pch=16, size = 3.5, alpha = 0.6, stroke = 0.5) +
  scale_color_manual(
    values = fungus_colors,
    labels = c(
      "S. americanus" = expression(italic("S. americanus")),
      "S. ampliporus" = expression(italic("S. ampliporus")),
      "S. clintonianus" = expression(italic("S. clintonianus")),
      "S. luteus" = expression(italic("S. luteus")),
      "S. spraguei" = expression(italic("S. spraguei")),
      "S. weaverae" = expression(italic("S. weaverae"))
    )
  ) +
  labs(
    x = NULL,
    y = expression(delta^{15}*"N (‰)")
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold", size = 12)
  )



#reformat to long data for the connected line plot
N_long <- data_filtered_fungi_clean %>%
  select(tubeLabel, species, year, xN_stipe, xN_cap) %>%
  pivot_longer(
    cols = c(xN_stipe, xN_cap),
    names_to = "tissue",
    values_to = "N_concentration"
  ) %>%
  mutate(
    tissue = factor(tissue, levels = c("xN_stipe", "xN_cap"), 
                    labels = c("Stipe", "Cap"))
  )

#plot cap vs. stipe total N concentration values
CS_N<- ggplot(N_long, aes(x = tissue, y = N_concentration, group = tubeLabel)) +
  #add colored lines connecting cap and stipe for each sample
  geom_line(aes(color = species), size = 0.7, alpha = 0.5) +
  #add points for each tissue with transparency
  geom_point(aes(color = species), pch=16, size = 3.5, alpha = 0.6, stroke = 0.5) +
  scale_color_manual(
    values = fungus_colors,
    labels = c(
      "S. americanus" = expression(italic("S. americanus")),
      "S. ampliporus" = expression(italic("S. ampliporus")),
      "S. clintonianus" = expression(italic("S. clintonianus")),
      "S. luteus" = expression(italic("S. luteus")),
      "S. spraguei" = expression(italic("S. spraguei")),
      "S. weaverae" = expression(italic("S. weaverae"))
    )
  ) +
  # Labels and theme
  labs(
    x = NULL,
    y = "[N]"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold", size = 12)
  )

CS_N


## plot relationship between cap and stipe n ##
#where'd15Ndif' is the difference between cap and stipe δ15N and 'n' is the average sporocarp nitrogen concentration
p_relationship <- ggplot(data_filtered_fungi_clean, aes(x = n, y = d15Ndif)) +
  #colored points by species
  geom_point(aes(color = species), pch=16, size = 3, alpha = 0.6) +
  #add an overall regression line
  geom_smooth(method = "lm", color = "black", fill = "gray50", alpha = 0.3, linewidth = 1) +
  scale_color_manual(
    values = fungus_colors,
    labels = c(
      "S. americanus" = expression(italic("S. americanus")),
      "S. ampliporus" = expression(italic("S. ampliporus")),
      "S. clintonianus" = expression(italic("S. clintonianus")),
      "S. luteus" = expression(italic("S. luteus")),
      "S. spraguei" = expression(italic("S. spraguei")),
      "S. weaverae" = expression(italic("S. weaverae"))
    )
  ) +
  labs(
    x = "[N]",
    y = expression("Cap-stipe "*delta^{15}*"N difference"),
    color = "Fungus species"
  ) +
  #add annotation showing the statistical relationship
  annotate(
    "text", 
    x = max(data_filtered_fungi_clean$n, na.rm = TRUE) * 0.8, 
    y = min(data_filtered_fungi_clean$d15Ndif, na.rm = TRUE) * 0.8,
    label = expression(atop(
      paste(italic("t"), " = -2.69"),
      paste(italic("p"), " = 0.008")
    )),
    hjust = 1,
    size = 4
  )+
  # Apply theme
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(face = "italic"),
    plot.title = element_text(size = 12)
  )

p_relationship


#plot relationship between cap-stipe δ15N difference and cap-stipe N concentration difference
p_relationship_xNdif <- ggplot(data_filtered_fungi_clean, aes(x = xNdif, y = d15Ndif)) +
  # Colored points by species
  geom_point(aes(color = species), pch=16, size = 3, alpha = 0.6) +
  # Add an overall regression line
  geom_smooth(method = "lm", color = "black", fill = "gray50", alpha = 0.3, linewidth = 1) +
  scale_color_manual(
    values = fungus_colors,
    labels = c(
      "S. americanus" = expression(italic("S. americanus")),
      "S. ampliporus" = expression(italic("S. ampliporus")),
      "S. clintonianus" = expression(italic("S. clintonianus")),
      "S. luteus" = expression(italic("S. luteus")),
      "S. spraguei" = expression(italic("S. spraguei")),
      "S. weaverae" = expression(italic("S. weaverae"))
    )
  ) +
  labs(
    x = "Cap-stipe [N] difference",
    y = expression("Cap-stipe "*delta^{15}*"N difference"),
    color = "Fungus species"
  ) +
  #add annotation showing the statistical relationship
  annotate(
    "text", 
    x = max(data_filtered_fungi_clean$xNdif, na.rm = TRUE) * 0.8, 
    y = min(data_filtered_fungi_clean$d15Ndif, na.rm = TRUE) * 0.8,
    label = expression(italic("t")~"= 1.93"),
    hjust = 1,
    size = 4
  ) +
  # p-value annotation on next line
  annotate(
    "text", 
    x = max(data_filtered_fungi_clean$xNdif, na.rm = TRUE) * 0.8, 
    y = min(data_filtered_fungi_clean$d15Ndif, na.rm = TRUE) * 0.8 - 0.2,
    label = expression(italic("p")~"= 0.056"),
    hjust = 1,
    size = 4
  )+
  
  # Apply theme
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.subtitle = element_text(face = "italic"),
    plot.title = element_text(size = 12)
  )

p_relationship_xNdif




##### RENDER PLOTS FOR MANUSCRIPT #####

#first plot side by side
row1 <- CS_N + CS_Ndelta + 
  plot_layout(widths = c(1, 1))

row2 <- p_relationship + p_relationship_xNdif + 
  plot_layout(widths = c(1, 1))

#dd the p_relationship plot below, spanning the full width
combined_plot <- row1 / row2 +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = 'A')

combined_plot

#save pdf
ggsave("isotopes_Sup.pdf", combined_plot, width = 8, height = 8)
#save SVG
ggsave("isotopes_Sup.svg", plot = combined_plot, width = 8, height = 6)


### Arrange plots for the manuscript ###

#first panel (A)
panel_A <- p_host + p_fungus + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A', tag_prefix = 'Panel ')

#second panel (B)
panel_B <- p_host_delta_15 + p_fungus_delta_15 + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A', tag_prefix = 'Panel ')

#note- were not just d15N_boxplot but d15N_boxplot_nolabs here to fix the spacing issue
#panel C 
panel_C <- d15N_boxplot_nolabs +
  plot_annotation(title = "C")

#panel D 
panel_D <- d13C_host_combined +
  plot_annotation(title = "D")

#panel E 
panel_E <- d13C_over_time_p +
  plot_annotation(title = "E")

#combine D and E side by side
panel_CD <- panel_C + panel_D + plot_layout(ncol = 2, widths = c(1, 1))

#combine all
all_panels <- (panel_A / panel_B) | (panel_CD / panel_E) + 
  plot_layout(widths = c(2, 2))

all_panels

#save pdf
ggsave("isotopes_figure_raw1.pdf", all_panels, width = 14, height = 8)
#save SVG
ggsave("isotopes_figure_raw1.svg", plot = all_panels, width = 14, height = 8)



