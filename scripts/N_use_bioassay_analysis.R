###################################################################################
#Analysis of Nitrogen use (biomass accumulation) across 6 genome-sequenced Suillus strains 
#all associated code can be found on in the git repo https://github.com/MycoPunk/Suillus_N_use
#this script was run in R version 4.4.1
#bioassays and data collection preformed by F. Maillard
#code and analysis written and preformed by L. Lofgren
#last updated 16.May.2025
###################################################################################

##### INITIATION AND DATA CLEANING #####

#set packages 
library(data.table)
library(tidyverse)
library(ggtree)
library(phytools)
library(ggplot2)
library(ggpubr)
library(devtools)
library(gridExtra)
library(dplyr)
library(tidyr)
library(broom)
library(purrr)
library(car)
library(emmeans)
library(patchwork)
library(cowplot)
sessionInfo()

#set seed for reproducibility
set.seed(666)

#Set dr 
setwd("~/Desktop/Project_N_use")

#set color scheme for species
fungus_colors <- c(
  "S. americanus" = "#5E3023",  # Rich brown
  "S. ampliporus" = "#D4AC6E",  # Warm sand/wheat
  "S. clintonianus" = "#355E2B",  # Richer forest green
  "S. luteus" = "#9B9B7A",  # Sage/olive
  "S. spraguei" = "#B07156",  # Terracotta/clay
  "S. weaverae" = "#6D5A72"   # Muted purple/aubergine
)

#read in data
data<- as.data.frame(fread("data/N_use_bioassay_clean_data.csv", header = TRUE, sep = ",")) #MEROPS BLASTP results

#get current raw names
unique(data$Species)

##remove species other than the 6 target species 
#set species to remove 
remove<- c("S_cothurnatus", "S_hirtellus", "S_lakei")
#remove rows from the data frame where "Species" matches any value in "remove"
data_filtered <- data %>%
  filter(!Species %in% remove)
#check
unique(data_filtered$Species) #looks good

#rename the species for figure generation
df_filtered <- data_filtered %>%
  mutate(Species = case_when(
    Species == "S_clinton" ~ "S. clintonianus",
    Species == "S_spraguei" ~ "S. spraguei",
    Species == "S_americanus" ~ "S. americanus",
    Species == "S_weaverae" ~ "S. weaverae",
    Species == "S_luteus" ~ "S. luteus",
    Species == "S_ampliporus" ~ "S. ampliporus",
    TRUE ~ as.character(Species)
  ))

unique(df_filtered$Species) #looks good

#subtract time point 0 zero values from each treatment - this is the weight of the plug
normalized_df <- df_filtered %>%
  #group by Species and Treatment
  group_by(Species, Treatment) %>%
  #get the initial biomass (at time 0) for each group
  mutate(initial_biomass = first(`Biomass (mg dry mass)`[`Time (days post inoculation)` == 0])) %>%
  #calculate the normalized mass
  mutate(mass_minus_plug = `Biomass (mg dry mass)` - initial_biomass) %>%
  #remove grouping
  ungroup() %>%
  #reorder columns
  select(ID, Species, `Biomass (mg dry mass)`, `Time (days post inoculation)`, Treatment, mass_minus_plug)

#for graphing, we want to show the difference between treatments relative to the control, so we have to split out the main data and the low ammonium control
ammonium_control <- normalized_df %>%
  filter(Treatment %in% 'Low ammonium control')

no_control <- normalized_df %>%
  filter(!Treatment %in% 'Low ammonium control')



### Make Plot - I'm calling this the lollybox ###

#for graphing, we want to look at the difference in growth rate relative to the low-ammonium control, so first we have to average the growth rate of the control,
#and then we're going to subtract that average for each time point from each replicate for each species.
#create ammonium_control_ave dataframe (for each treatment, for each species, averaged over each time point)
ammonium_control_ave <- ammonium_control %>%
  group_by(Species, `Time (days post inoculation)`) %>%
  summarize(
    average_mass_minus_plug = mean(mass_minus_plug, na.rm = TRUE),
    .groups = "drop"
  )

#note, that a couple of these values are slightly negative, which is just measurement error - to graph this, set anything less than 0 to 0. 
ammonium_control_ave <- ammonium_control_ave %>%
  mutate(average_mass_minus_plug = pmax(average_mass_minus_plug, 0))


#now create a new data frame by joining no_control with ammonium_control_ave, and subtracting the correct ave control value for that time point
no_control_normalized <- no_control %>%
  #join with ammonium_control_ave based on Species and Time
  left_join(
    #select only the columns needed from ammonium_control_ave
    ammonium_control_ave %>%
      select(Species, `Time (days post inoculation)`, average_mass_minus_plug),
    #specify the joining columns
    by = c("Species", "Time (days post inoculation)")
  ) %>%
  #create the new column by subtracting the average ammonium control value
  mutate(
    mass_minus_plug_normalized = mass_minus_plug - average_mass_minus_plug
  )

#filter out time point 0, because we don't want it in the graph
no_control_normalized_filtered <- no_control_normalized %>%
  filter(`Time (days post inoculation)` > 0)

#determine the global min and max for consistent scales
#global_min <- min(no_control_normalized_filtered$mass_minus_plug_normalized, na.rm = TRUE)
#global_max <- max(no_control_normalized_filtered$mass_minus_plug_normalized, na.rm = TRUE)

#get the filtered time points for shapes
time_points <- unique(no_control_normalized_filtered$`Time (days post inoculation)`)
time_points <- sort(time_points)  #ensure they're in order

#map time points to shapes 
time_shapes <- rep(15:25, length.out = length(time_points))
names(time_shapes) <- as.character(time_points)

#set species order
species_order <- sort(unique(no_control_normalized_filtered$Species))

#filter data for each treatment
ammonium_data <- no_control_normalized_filtered %>% filter(Treatment == "Ammonium")
bsa_data <- no_control_normalized_filtered %>% filter(Treatment == "BSA")
bsa_tannin_data <- no_control_normalized_filtered %>% filter(Treatment == "BSA-Tannin")
chitin_data <- no_control_normalized_filtered %>% filter(Treatment == "Chitin")

#filter segment for each treatment
ammonium_segments <- segment_data %>% filter(Treatment == "Ammonium")
bsa_segments <- segment_data %>% filter(Treatment == "BSA")
bsa_tannin_segments <- segment_data %>% filter(Treatment == "BSA-Tannin")
chitin_segments <- segment_data %>% filter(Treatment == "Chitin")




### Plot that shit ###

#define explicit breaks for each x-axis with the new limits
x_min <- -2.5
ammonium_max <- 10
ammonium_breaks <- seq(x_min, ammonium_max, by = 2)
bsa_max <- 5
bsa_breaks <- seq(x_min, bsa_max, by = 1)
bsa_tannin_max <- 2.5
bsa_tannin_breaks <- seq(x_min, bsa_tannin_max, by = 1)
chitin_max <- 2.5
chitin_breaks <- seq(x_min, chitin_max, by = 1)

#create base plot function
create_panel <- function(data, segments_data, title, x_max, breaks, show_y_labels = FALSE) {
  p <- ggplot(data, aes(y = factor(Species, levels = rev(species_order)))) +
    geom_segment(data = segments_data,
                 aes(x = 0, xend = mean_x, y = y, yend = y),
                 color = "black", linewidth = 0.7, alpha = 0.5) +
    geom_point(aes(x = mass_minus_plug_normalized, 
                   color = Species, 
                   shape = factor(`Time (days post inoculation)`)),
               position = position_jitter(height = 0.15, width = 0),
               size = 3.2, alpha = 0.5) +
    stat_summary(aes(x = mass_minus_plug_normalized), 
                 fun = mean, geom = "point", 
                 shape = 18, size = 4, color = "black") +
    stat_summary(aes(x = mass_minus_plug_normalized), 
                 fun.data = "mean_se", geom = "errorbar",
                 width = 0.3, color = "black", linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_x_continuous(limits = c(x_min, x_max), 
                       breaks = breaks,
                       expand = c(0, 0)) +
    scale_color_manual(values = fungus_colors) +
    scale_shape_manual(values = time_shapes) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
      legend.position = "none",
      #note minimal consistent margins
      plot.margin = margin(5, 5, 5, 5)
    ) +
    ggtitle(title)
  
  #toggle y-axis text visibility
  if (!show_y_labels) {
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  } else {
    p <- p + theme(
      axis.text.y = element_text(size = 10, face = "italic")
    )
  }
  
  return(p)
}



#render the plots
p1 <- create_panel(ammonium_data, ammonium_segments, "Ammonium", ammonium_max, ammonium_breaks, TRUE)
p2 <- create_panel(bsa_data, bsa_segments, "BSA", bsa_max, bsa_breaks, FALSE)
p3 <- create_panel(bsa_tannin_data, bsa_tannin_segments, "BSA-Tannin", bsa_tannin_max, bsa_tannin_breaks, FALSE)
p4 <- create_panel(chitin_data, chitin_segments, "Chitin", chitin_max, chitin_breaks, FALSE)


#add the legend
#first combine the 4 plots using patchwork
top_row <- p1 + p2 + p3 + p4 + plot_layout(ncol = 4)

#create the legend manually 
direct_legend <- ggplot() +
  geom_point(aes(x = 1:6, y = rep(1, 6), 
                 color = factor(c("S. americanus", "S. ampliporus", "S. clintonianus", 
                                  "S. luteus", "S. spraguei", "S. weaverae"))),
             shape = 16, size = 3) +
  geom_point(aes(x = 8:10, y = rep(1, 3),
                 shape = factor(c(12, 24, 35))),
             color = "black", size = 3) +
  scale_color_manual(values = fungus_colors, name = "Species") +
  scale_shape_manual(values = time_shapes, name = "Time (days)") +
  theme_void() +
  theme(legend.position = "bottom")

#use cowplot to stack them
final_plot <- plot_grid(
  top_row,                    # The combined row of 4 plots
  direct_legend,              # The legend plot
  ncol = 1,                   # Stack vertically
  rel_heights = c(10, 2),     # Relative heights (legend smaller)
  align = 'v',                # Align vertically
  axis = 'lr'                 # Align left and right edges
)


final_plot

#save plot
#ggsave("bioassay_figure_raw.pdf", final_plot, width = 10, height = 7)
ggsave("bioassay_figure_raw.pdf.svg", plot = final_plot, width = 10, height = 7)




### STATS ###

#convert Species and Treatment to factors
no_control_normalized_filtered$Species <- as.factor(no_control_normalized_filtered$Species)
no_control_normalized_filtered$Treatment <- as.factor(no_control_normalized_filtered$Treatment)


#Set up Two-way ANOVA
#using Type III sum of squares to handle unbalanced design
model <- aov(mass_minus_plug_normalized ~ Species * Treatment, 
             data = no_control_normalized_filtered)

#treatment has a much larger effect than species or species:treatment
#Terms:
#  Species Treatment Species:Treatment Residuals
#Sum of Squares   27.2361  630.9746           46.6174  239.1635
#Deg. of Freedom        5         3                15       192


#check ANOVA assumptions
#check residuals for normality
qqnorm(residuals(model))
qqline(residuals(model))
#looks pretty normal 

#check homogeneity of variances
plot(model, 1)  #residuals vs Fitted plot indicates no serious heteroscedasticity


#ANOVA summary with Type III SS to include interaction effects
Anova(model, type = "III")
summary(model)

# Df Sum Sq Mean Sq F value   Pr(>F)    
# Species             5   27.2    5.45   4.373 0.000858 ***
#   Treatment           3  631.0  210.32 168.848  < 2e-16 ***
#   Species:Treatment  15   46.6    3.11   2.495 0.002211 ** 
#   Residuals         192  239.2    1.25                     
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#now calculate effect sizes (eta squared) to compare the relative importance of species vs treatment
eta_squared <- function(aov_model) {
  ss_effect <- unlist(summary(aov_model)[[1]]["Sum Sq"])
  ss_total <- sum(ss_effect)
  return(ss_effect / ss_total)
}

eta_sq <- eta_squared(model)
eta_sq_df <- data.frame(
  Effect = c("Species", "Treatment", "Species:Treatment", "Residuals"),
  Eta_Squared = eta_sq
)
eta_sq_df

# Effect Eta_Squared
# Sum Sq1           Species  0.02885210
# Sum Sq2         Treatment  0.66841123
# Sum Sq3 Species:Treatment  0.04938324
# Sum Sq4         Residuals  0.25335343


#preform post-hoc tests
# 1. For Species effect
species_emmeans <- emmeans(model, ~ Species)
species_pairs <- pairs(species_emmeans, adjust = "tukey")
species_pairs

#get treatment effects
treatment_emmeans <- emmeans(model, ~ Treatment)
treatment_pairs <- pairs(treatment_emmeans, adjust = "tukey")
treatment_pairs

#get interaction effecgts (effect of Treatment within each Species)
interaction_emmeans <- emmeans(model, ~ Treatment | Species)
interaction_pairs <- pairs(interaction_emmeans, adjust = "tukey")
interaction_pairs

#summary statistics
group_means <- no_control_normalized_filtered %>%
  group_by(Species, Treatment) %>%
  summarise(
    mean_biomass = mean(mass_minus_plug_normalized, na.rm = TRUE),
    sd_biomass = sd(mass_minus_plug_normalized, na.rm = TRUE),
    n = n(),
    se_biomass = sd_biomass / sqrt(n),
    .groups = "drop"
  )

print(n= 100, group_means)


#look specifically at Ammonium performance: 
#create a subset of data containing only Ammonium treatment
Ammonium_subset <- no_control_normalized_filtered %>%
  filter(Treatment == "Ammonium")

#run a one-way ANOVA for species effect within Ammonium
Ammonium_model <- aov(mass_minus_plug_normalized ~ Species, data = Ammonium_subset)
summary(Ammonium_model)

#run post-hoc pairwise comparisons with Tukey adjustment
Ammonium_tukey <- TukeyHSD(Ammonium_model)
print(Ammonium_tukey)

#look specifically at BSA-Tannin performance: 
#create a subset of data containing only BSA-Tannin treatment
BSA_Tannin_subset <- no_control_normalized_filtered %>%
  filter(Treatment == "BSA-Tannin")

#run a one-way ANOVA for species effect within BSA-Tannin
BSA_Tannin_model <- aov(mass_minus_plug_normalized ~ Species, data = BSA_Tannin_subset)
summary(BSA_Tannin_model)

# Run post-hoc pairwise comparisons with Tukey adjustment
BSA_Tannin_tukey <- TukeyHSD(BSA_Tannin_model)
print(BSA_Tannin_tukey)


##look specifically at BSA performance: 
#create a subset of data containing only BSA treatment
BSA_subset <- no_control_normalized_filtered %>%
  filter(Treatment == "BSA")

#run a one-way ANOVA for species effect within BSA
BSA_model <- aov(mass_minus_plug_normalized ~ Species, data = BSA_subset)
summary(BSA_model)

#run post-hoc pairwise comparisons with Tukey adjustment
BSA_tukey <- TukeyHSD(BSA_model)
print(BSA_tukey)


#look at Chitin performance: 
#create a subset of data containing only Chitin treatment
Chitin_subset <- no_control_normalized_filtered %>%
  filter(Treatment == "Chitin")

#run a one-way ANOVA for species effect within Chitin
Chitin_model <- aov(mass_minus_plug_normalized ~ Species, data = Chitin_subset)
summary(Chitin_model)

#run post-hoc pairwise comparisons with Tukey adjustment
Chitin_tukey <- TukeyHSD(Chitin_model)
print(Chitin_tukey)
