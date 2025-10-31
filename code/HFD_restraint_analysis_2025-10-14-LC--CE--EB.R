#######################################
###                                 ###
###       Sciolino and Engburg      ###
### Statistical Consulting Services ###
###                                 ###
#######################################

# Load packages
library(dplyr)
library(tidyr)
library(glmmTMB)
library(emmeans)
library(readxl)
library(performance)
library(ggeffects)
library(sjPlot)
library(glmmTMB)
library(ggplot2)
library(bayesplot)
library(viridis)
library(kableExtra)
library(ggsignif)
library(ggpubr)
library(ggh4x)
library(extrafont)
library(tidyverse)
library(writexl)
library(cowplot)
library(patchwork)
library(remotes)

# remotes::install_version("Rttf2pt1", version = "1.3.8")
extrafont::font_import()

# Load fonts
loadfonts(device = "postscript")
fonts()

# Load and preprocess data
histology <- read.csv("SummarizedData_AllCellsAllRegions_8.21.25.csv", header = TRUE) %>%
  mutate(ID = factor(ID), # ID into factor
         Anatomy2 = Anatomy) %>% # Create duplicate for Anatomy (see later)
  mutate(
    # Extract Group as "HFD" or "SD"
    Group = case_when(
      str_detect(Condition, regex("HFD", ignore_case = TRUE)) ~ "HFD",
      str_detect(Condition, regex("SD", ignore_case = TRUE)) ~ "SD",
      TRUE ~ NA_character_
    ),
    # Extract Condition as "restraint" or "no restraint"
    Condition = case_when(
      str_detect(Condition, regex("no restraint", ignore_case = TRUE)) ~ "no restraint",
      str_detect(Condition, regex("restraint", ignore_case = TRUE)) ~ "restraint",
      TRUE ~ NA_character_
    ),
    # Group and Condition into factor
    Group = factor(Group, levels = c("SD", "HFD")), 
    Condition = factor(Condition, levels = c("no restraint", "restraint"))
  ) %>%
  separate(Anatomy2, c("Sample", "Side"), extra = "drop") %>% # Using the "Anatomy2" we created to split into Sample and Size
  filter(!(Sample %in% c("subC", "SubC")))  # Exclude any Sample containing "SubC", as there are not enough samples for stats

# Create proportions data set
proportions <- histology %>%
  group_by(ID, Sample, Group, Condition) %>%
  summarize(N = n(), 
            n = sum(ifelse(fos_num > 0, 1, 0)), 
            P = n/N) %>%
  ungroup() %>%
  arrange(Sample)

# Subset the ID and Sex
histology_sex <- histology %>%
  dplyr::select(ID, Sex) %>%
  distinct()

# Join the Sex data to the histology data set
proportions_sex <- proportions %>%
  left_join(histology_sex, by = "ID")

# Add colors for boxplots
proportions_sex_col <- proportions_sex %>%
  mutate(PointColors = ifelse(Group == "SD", "#414141", "#006EAF"),
         PointShapes = ifelse(Sex == "Female", 16, 1),
         PointFills = ifelse(Group == "SD", "#dbdbdb", "#79e9fb")) %>%
  mutate(Condition_Group = interaction(Condition, Group, sep = " + "))

# Condition labels, used for plotting
Condition_levels <- levels(proportions_sex_col$Condition)
Geno_labs <- c("no restraint + SD" = expression(italic(No)~italic(restraint)~italic("+")~italic(SD)), 
               "restraint + SD" = expression(italic(Restraint)~italic("+")~italic(SD)), 
               "no restraint + HFD" = expression(italic(No)~italic(restraint)~italic("+")~italic(HFD)), 
               "restraint + HFD" = expression(italic(Restraint)~italic("+")~italic(HFD)))


#############################################

########################
###                  ###
### With TH+ control ###
###                  ###
########################

# Merge proportions data to histology
histology_p <- histology %>%
  left_join(proportions, by = c("ID", "Sample", "Group", "Condition")) %>%
  mutate(n_no_dots = N - n, 
         Sample = as.factor(Sample), 
         Condition = as.factor(Condition), 
         Sex = as.factor(Sex)) %>%
  dplyr::select(n, n_no_dots, Condition, Group, Sample, Sex, ID) %>%
  distinct()

# Define contrasts
contrasts(histology_p$Condition) = contr.sum(levels(histology_p$Condition))
contrasts(histology_p$Sample) = contr.sum(levels(histology_p$Sample))
contrasts(histology_p$Group) = contr.sum(levels(histology_p$Group))

# Binomial mixed effects regression model with just Group as predictor -- DELTED

# Binomial mixed effects regression model with Condition as predictor

## Notice the mutate(factor) in glmmTMB has been dropped -EB

# C1/A1
histology_n <- histology %>%
  mutate(dots = ifelse(fos_num > 0, 1, 0), 
         mouse = ID) 

histology_n %>% dplyr::filter(Sample == "A1")

library(lme4)
glm.hist_n_A1 <- glmmTMB(dots ~ Condition * Group + (1|ID), 
                         family = binomial(link = "logit"), 
                         data = histology_n %>% dplyr::filter(Sample == "A1"))
summary(glm.hist_n_A1)

# Notice Condition + Group. Interaction removed. 
glm.hist_p_A1 <- glmmTMB(cbind(n, n_no_dots) ~ Condition + Group,# + (1|ID), 
                         family = binomial(link = "logit"), 
                         data = histology_p %>% dplyr::filter(Sample == "A1"))
summary(glm.hist_p_A1)
car::Anova(glm.hist_p_A1) # Chi-square values

# C2/A2
glm.hist_p_A2 <- glmmTMB(cbind(n, n_no_dots) ~ Condition * Group,# + (1|ID), 
                         family = binomial(link = "logit"), 
                         data = histology_p %>% dplyr::filter(Sample == "A2"))
car::Anova(glm.hist_p_A2) # Chi-square values


# LC
glm.hist_p_LC <- glmmTMB(cbind(n, n_no_dots) ~ Condition * Group,# + (1|ID), 
                         family = binomial(link = "logit"), 
                         data = histology_p %>% dplyr::filter(Sample == "LC"))
car::Anova(glm.hist_p_LC) # Chi-square values

### Summary of the model ###

# A1
summary(glm.hist_p_A1)
tab_model(glm.hist_p_A1, show.se = TRUE, show.stat = TRUE)
car::Anova(glm.hist_p_A1, type = "III")

# A2
summary(glm.hist_p_A2)
tab_model(glm.hist_p_A2, show.se = TRUE, show.stat = TRUE)
car::Anova(glm.hist_p_A2, type = "III")

# LC
summary(glm.hist_p_LC)
tab_model(glm.hist_p_LC, show.se = TRUE, show.stat = TRUE)
car::Anova(glm.hist_p_LC, type = "III")


### Contrasts ###

## No multiple comparisons adjustment

# # Group without sex
# 
# # C1/A1
# contrasts.hist_p_A1_grp <- emmeans(glm.hist_p_A1_grp, pairwise ~ Group, 
#                                    adjust = "none")
# contrasts.hist_p_A1_grp.df <- contrasts.hist_p_A1_grp$contrasts %>% 
#   as.data.frame()
# 
# # C2/A2
# contrasts.hist_p_A2_grp <- emmeans(glm.hist_p_A2_grp, pairwise ~ Group, 
#                                    adjust = "none")
# contrasts.hist_p_A2_grp.df <- contrasts.hist_p_A2_grp$contrasts %>% 
#   as.data.frame()
# 
# # LC
# contrasts.hist_p_LC_grp <- emmeans(glm.hist_p_LC_grp, pairwise ~ Group, 
#                                    adjust = "none")
# contrasts.hist_p_LC_grp.df <- contrasts.hist_p_LC_grp$contrasts %>% 
#   as.data.frame()


# Condition * Group

# C1/A1
contrasts.hist_p_A1 <- emmeans(glm.hist_p_A1, pairwise ~ Condition * Group, 
                               adjust = "none")
contrasts.hist_p_A1.df <- contrasts.hist_p_A1$contrasts %>% 
  as.data.frame()
ggemmeans(glm.hist_p_A1, terms = c("Condition", "Group")) %>% 
  plot() + 
  labs(y = "% TH+ Cells") +
  ggtitle("C1/A1") +
  scale_x_discrete(limits = Condition_levels, # Fixed this - EB
                   labels = Geno_labs) + # Fixed this - EB
  theme(text = element_text(size = 10, family = "sans"), 
        axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 0.8, 
                                   size = 10), 
        plot.margin = unit(c(5.5, 5.5, 5.5, 8.5), "points"))
ggsave(filename = "contrasts_hist_p_A1.png", 
       width = 4.5, height = 3, units = "in")
print(contrasts.hist_p_A1.df)  

# C2/A2
contrasts.hist_p_A2 <- emmeans(glm.hist_p_A2, pairwise ~ Condition * Group, 
                               adjust = "none")
contrasts.hist_p_A2.df <- contrasts.hist_p_A2$contrasts %>% 
  as.data.frame()
ggemmeans(glm.hist_p_A2, terms = c("Condition", "Group")) %>%
  plot() + 
  labs(y = "% TH+ Cells") +
  ggtitle("C2/A2") +
  scale_x_discrete(limits = Condition_levels, # Fixed this - EB
                   labels = Geno_labs) + # Fixed this - EB
  theme(text = element_text(size = 10, family = "sans"), 
        axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 0.8, 
                                   size = 10), 
        plot.margin = unit(c(5.5, 5.5, 5.5, 8.5), "points"))
ggsave(filename = "contrasts_hist_p_A2.png", 
       width = 4.5, height = 3, units = "in")
print(contrasts.hist_p_A2.df)  

# LC
contrasts.hist_p_LC <- emmeans(glm.hist_p_LC, pairwise ~ Condition * Group, 
                               adjust = "none")
contrasts.hist_p_LC.df <- contrasts.hist_p_LC$contrasts %>% 
  as.data.frame()
ggemmeans(glm.hist_p_LC, terms = c("Condition", "Group")) %>%
  plot() + 
  labs(y = "% TH+ Cells") +
  ggtitle("LC") +
  scale_x_discrete(limits = Condition_levels, # Fixed this - EB
                   labels = Geno_labs) + # Fixed this - EB
  theme(text = element_text(size = 10, family = "sans"), 
        axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 0.8, 
                                   size = 10), 
        plot.margin = unit(c(5.5, 5.5, 5.5, 8.5), "points"))
ggsave(filename = "contrasts_hist_p_LC.png", 
       width = 4.5, height = 3, units = "in")

print(contrasts.hist_p_LC.df)  


### Check Model Fit & Prediction ###

# C1/A1
# dev.new()
performance::check_model(glm.hist_p_A1)

sim_A1 <- DHARMa::simulateResiduals(glm.hist_p_A1)
plot(sim_A1)
pp_check(glm.hist_p_A1, type = "discrete_dots")

# C2/A2
# dev.new()
performance::check_model(glm.hist_p_A2)

sim_A2 <- DHARMa::simulateResiduals(glm.hist_p_A2)
plot(sim_A2)
pp_check(glm.hist_p_A2, type = "discrete_dots")

# LC
# dev.new()
performance::check_model(glm.hist_p_LC)

sim_LC <- DHARMa::simulateResiduals(glm.hist_p_LC)
plot(sim_LC)
pp_check(glm.hist_p_LC, type = "discrete_dots")


##########################################################

################
### Boxplots ###
################

# Create combined factor
proportions_sex_col <- proportions_sex_col %>%
  filter(Sample != "subC") %>% # REmoving Sub C
  group_by(Sample, Group, Condition) %>%
  mutate(NN = sum(N), 
         nn = sum(n), 
         PP = nn/NN) %>%
  ungroup()

histology %>%
  dplyr::select(ID, Sample, Group) %>%
  distinct() %>%
  group_by(Sample, Group) %>%
  summarize(N = n())

# Generate A1, A2 and LC plots by Group
SampleUnique <- unique(proportions_sex_col$Sample)
SampleTitle <- c("C1/A1", "C2/A2", "LC")

# Condition * Group

i <- 1
max_y <- 40 # Here we define maximum y. Useful in plotting
tick_int_small <- max_y/4
# Generate C1/A1, C2/A2 and LC plots by Condition * Group
for (samp in SampleUnique) {
  
  # Filter by sample
  proportions_sex_col_samp <- proportions_sex_col %>%
    filter(Sample == samp)
  
  # Create combined factor with desired order and more spacing
  proportions_sex_col_samp <- proportions_sex_col_samp %>%
    mutate(Condition_Group = factor(interaction(Condition, Group, sep = " + "),
                                    levels = c("no restraint + SD", "no restraint + HFD", 
                                               "restraint + SD", "restraint + HFD")),
           # Create numeric positions with gap between conditions
           x_position = case_when(
             Condition_Group == "no restraint + SD" ~ 1,
             Condition_Group == "no restraint + HFD" ~ 2,
             Condition_Group == "restraint + SD" ~ 3.5,
             Condition_Group == "restraint + HFD" ~ 4.5
           ))
  
  # Filter by proportion by Condition * Group, for the contingency plot
  proportions_sex_col_prop <- proportions_sex_col_samp %>%
    dplyr::select(Sample, Condition_Group, NN, nn, PP) %>%
    distinct() %>%
    arrange(Condition_Group)
  
  proportions_sex_col_prop <- 
    proportions_sex_col_prop %>%
    mutate(Condition_Group2 = paste(Condition_Group, "% Fos+"))
  
  # Also create a count table per condition-group combination
  proportions_sex_col_samp_n <- proportions_sex_col_samp %>%
    count(Condition_Group, name = "N")
  
  # Defining y-axis tick
  tick_y <- rep("", ceiling(max_y/tick_int_small * 1.5) + 1)
  tick_y[2*(seq(0, max_y/(2*tick_int_small))) + 1] <- as.character(seq(0, max_y, by = 2*tick_int_small))
  
  if (i == 1) {
    test_result1 <- contrasts.hist_p_A1.df[4,] # restraint SD - no restraint HFD
    test_result2 <- contrasts.hist_p_A1.df[3,] # no restraint SD - restraint HFD
    test_result3 <- contrasts.hist_p_A1.df[6,] # restraint HFD - no restraint HFD 
    test_result4 <- contrasts.hist_p_A1.df[2,] # no restraint SD - no restraint HFD
    test_result5 <- contrasts.hist_p_A1.df[5,] # restraint SD - restraint HFD
    test_result6 <- contrasts.hist_p_A1.df[1,] # restraint SD - no restraint SD
  }
  else if (i == 2) {
    test_result1 <- contrasts.hist_p_A2.df[4,]
    test_result2 <- contrasts.hist_p_A2.df[3,]
    test_result3 <- contrasts.hist_p_A2.df[6,]
    test_result4 <- contrasts.hist_p_A2.df[2,]
    test_result5 <- contrasts.hist_p_A2.df[5,]
    test_result6 <- contrasts.hist_p_A2.df[1,]
  }
  else {
    test_result1 <- contrasts.hist_p_LC.df[4,]
    test_result2 <- contrasts.hist_p_LC.df[3,]
    test_result3 <- contrasts.hist_p_LC.df[6,]
    test_result4 <- contrasts.hist_p_LC.df[2,]
    test_result5 <- contrasts.hist_p_LC.df[5,]
    test_result6 <- contrasts.hist_p_LC.df[1,]
  }
  
  # Boxplots (blocked off)
  pp <- proportions_sex_col_samp %>%
    ggplot(aes(x = Condition_Group, y = P * 100, color = Condition_Group, fill = Condition_Group)) +
    geom_boxplot(outliers = FALSE, size = 0.25, show.legend = FALSE, width = 0.6) +
    geom_point(aes(shape = Sex, stroke = ifelse(Sex == "Male", 0.5, 0),
                   size = ifelse(Sex == "Male", 1.0, 1.5)),
               position = position_jitter(width = 0.15, height = 0),
               show.legend = FALSE) +
    scale_shape_manual(values = c("Male" = 21, "Female" = 19)) +  # 21 = fillable circle, 19 = solid circle
    scale_size_identity() +  # Use the actual size values specified
    geom_col(data = data.frame(P = rep(0, 4),
                               Condition_Group = factor(c("no restraint + SD", "no restraint + HFD",
                                                          "restraint + SD", "restraint + HFD"),
                                                        levels = c("no restraint + SD", "no restraint + HFD",
                                                                   "restraint + SD", "restraint + HFD"))),
             aes(x = Condition_Group, y = P, color = Condition_Group, fill = Condition_Group),
             position = "dodge", linewidth = 0.25, show.legend = TRUE) +
    geom_text(data = proportions_sex_col_samp_n, color = "black",
              aes(x = Condition_Group, y = 5, label = N), size = 7/.pt) +
  geom_signif(comparisons = list(c("restraint + SD", "no restraint + HFD")),
              map_signif_level = TRUE,
              annotations = ifelse(test_result1$p.value < 0.001, "***",
                                   ifelse(test_result1$p.value < 0.01, "**",
                                          ifelse(test_result1$p.value < 0.05, "*",
                                                 ifelse(test_result1$p.value < 0.15,
                                                        paste("P =", round(test_result1$p.value, 2)),
                                                        "ns")))),
              size = 0.25,
              tip_length = 0,
              color = "black",
              textsize = ifelse(test_result1$p.value < 0.05, 16/.pt, 7/.pt),
              vjust = ifelse(test_result1$p.value < 0.05, 0.5, 0),
              y_position = max_y,
              fontface = ifelse(test_result1$p.value < 0.05, "bold", "plain")) +
    geom_signif(comparisons = list(c("no restraint + SD", "restraint + HFD")),
                map_signif_level = TRUE,
                annotations = ifelse(test_result2$p.value < 0.001, "***",
                                     ifelse(test_result2$p.value < 0.01, "**",
                                            ifelse(test_result2$p.value < 0.05, "*",
                                                   ifelse(test_result2$p.value < 0.15,
                                                          paste("P =", round(test_result2$p.value, 2)),
                                                          "ns")))),
                size = 0.25,
                tip_length = 0,
                color = "black",
                textsize = ifelse(test_result2$p.value < 0.05, 16/.pt, 7/.pt),
                vjust = ifelse(test_result2$p.value < 0.05, 0.5, 0),
                y_position = max_y*1.1,
                fontface = ifelse(test_result2$p.value < 0.05, "bold", "plain")) +
    geom_signif(comparisons = list(c("restraint + HFD", "no restraint + HFD")),
                map_signif_level = TRUE,
                annotations = ifelse(test_result3$p.value < 0.001, "***",
                                     ifelse(test_result3$p.value < 0.01, "**",
                                            ifelse(test_result3$p.value < 0.05, "*",
                                                   ifelse(test_result3$p.value < 0.15,
                                                          paste("P =", round(test_result3$p.value, 2)),
                                                          "ns")))),
                size = 0.25,
                tip_length = 0,
                color = "black",
                textsize = ifelse(test_result3$p.value < 0.05, 16/.pt, 7/.pt),
                vjust = ifelse(test_result3$p.value < 0.05, 0.5, 0),
                y_position = max_y*1.2,
                fontface = ifelse(test_result3$p.value < 0.05, "bold", "plain")) +
    geom_signif(comparisons = list(c("no restraint + SD", "no restraint + HFD")),
                map_signif_level = TRUE,
                annotations = ifelse(test_result4$p.value < 0.001, "***",
                                     ifelse(test_result4$p.value < 0.01, "**",
                                            ifelse(test_result4$p.value < 0.05, "*",
                                                   ifelse(test_result4$p.value < 0.15,
                                                          paste("P =", round(test_result4$p.value, 2)),
                                                          "ns")))),
                size = 0.25,
                tip_length = 0,
                color = "black",
                textsize = ifelse(test_result4$p.value < 0.05, 16/.pt, 7/.pt),
                vjust = ifelse(test_result4$p.value < 0.05, 0.5, 0),
                y_position = max_y*1.3,
                fontface = ifelse(test_result4$p.value < 0.05, "bold", "plain")) +
    geom_signif(comparisons = list(c("restraint + SD", "restraint + HFD")),
                map_signif_level = TRUE,
                annotations = ifelse(test_result5$p.value < 0.001, "***",
                                     ifelse(test_result5$p.value < 0.01, "**",
                                            ifelse(test_result5$p.value < 0.05, "*",
                                                   ifelse(test_result5$p.value < 0.15,
                                                          paste("P =", round(test_result5$p.value, 2)),
                                                          "ns")))),
                size = 0.25,
                tip_length = 0,
                color = "black",
                textsize = ifelse(test_result5$p.value < 0.05, 16/.pt, 7/.pt),
                vjust = ifelse(test_result5$p.value < 0.05, 0.5, 0),
                y_position = max_y*1.4,
                fontface = ifelse(test_result5$p.value < 0.05, "bold", "plain")) +
    geom_signif(comparisons = list(c("restraint + SD", "no restraint + SD")),
                map_signif_level = TRUE,
                annotations = ifelse(test_result6$p.value < 0.001, "***",
                                     ifelse(test_result6$p.value < 0.01, "**",
                                            ifelse(test_result6$p.value < 0.05, "*",
                                                   ifelse(test_result6$p.value < 0.15,
                                                          paste("P =", round(test_result6$p.value, 2)),
                                                          "ns")))),
                size = 0.25,
                tip_length = 0,
                color = "black",
                textsize = ifelse(test_result6$p.value < 0.05, 16/.pt, 7/.pt),
                vjust = ifelse(test_result6$p.value < 0.05, 0.5, 0),
                y_position = max_y*1.5,
                fontface = ifelse(test_result6$p.value < 0.05, "bold", "plain")) +
    # Custom axis lines
    geom_segment(x = 0.4, y = 0, xend = 0.4, yend = max_y,
                 inherit.aes = FALSE, linewidth = 0.25) + # y-axis
    geom_segment(x = 0.2, y = 0, xend = 5.3, yend = 0,
                 inherit.aes = FALSE, linewidth = 0.25) + # x-axis
    # Y-axis tick marks every 1/4 of the max_y
    geom_segment(x = 0.3, y = max_y/4, xend = 0.4, yend = max_y/4,
                 inherit.aes = FALSE, linewidth = 0.25) +
    geom_segment(x = 0.3, y = max_y/2, xend = 0.4, yend = max_y/2,
                 inherit.aes = FALSE, linewidth = 0.25) +
    geom_segment(x = 0.3, y = 3*max_y/4, xend = 0.4, yend = 3*max_y/4,
                 inherit.aes = FALSE, linewidth = 0.25) +
    geom_segment(x = 0.3, y = max_y, xend = 0.4, yend = max_y,
                 inherit.aes = FALSE, linewidth = 0.25) +
    # X-axis tick marks for conditions (with more spacing)
    geom_segment(x = 1.5, y = -2, xend = 1.5, yend = 0,
                 inherit.aes = FALSE, linewidth = 0.25) + # No restraint center
    geom_segment(x = 4, y = -2, xend = 4, yend = 0,
                 inherit.aes = FALSE, linewidth = 0.25) + # Restraint center
    theme_classic(base_size = 10) +
    scale_x_discrete("",
                     breaks = c("no restraint + SD", "no restraint + HFD",
                                "restraint + SD", "restraint + HFD"),
                     labels = c("", "", "", "")) +
    # Add condition labels manually
    annotate("text", x = 1.5, y = -5, label = "No Restraint",
             angle = 45, hjust = 1, vjust = 1, size = 10/.pt) +
    annotate("text", x = 4, y = -5, label = "Restraint",
             angle = 45, hjust = 1, vjust = 1, size = 10/.pt) +
    scale_y_continuous("Fos+ Cells \\
                       (% TH+ Cells)",
                       limits = c(0, tick_int_small*ceiling(max_y/tick_int_small*1.5)),
                       breaks = seq(0, tick_int_small*ceiling(max_y/tick_int_small*1.5), tick_int_small),
                       labels = tick_y) + # Where the y-tick comes in
    scale_color_manual(values = c("#6f6f6f", "#0698ff", "#6f6f6f", "#0698ff"),
                       labels = Geno_labs) +
    scale_fill_manual(values = c("#dbdbdb", "#79e9fb", "#dbdbdb", "#79e9fb"),
                      labels = Geno_labs) +
    theme(text = element_text(size = 10),
          axis.text = element_text(size = unit(10, "pt"), color = "black"),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 1.0, 
                                     size = 10), 
          # axis.text.x = element_blank(), # Remove default x-axis text since we're using custom annotations
          axis.title.y = ggtext::element_markdown(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, margin = margin(0, 0, 10, 0)),
          plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points"),
          legend.position = ifelse(i == 3, "right", "none"),
          legend.box.spacing = unit(0, "pt"),
          legend.key.width = unit(0.5, "lines"),
          legend.key.height = unit(0.5, "lines"),
          legend.text = element_text(margin = margin(r = 2)),
          axis.line = element_blank(),
          axis.ticks.length.y = unit(0, "pt"),
          axis.ticks.length.x = unit(0, "pt")) +    
    coord_cartesian(xlim = c(0.2, 5.3), 
                    ylim = c(0, tick_int_small*ceiling(max_y/tick_int_small*1.5)),
                    expand = FALSE, clip = "off") +
    force_panelsizes(rows = unit(1.5, "in"),
                     cols = unit(0.4, "in")) +
    ggtitle(SampleTitle[i])

  ggsave(paste0("Boxplot_Control_Conditions_", samp, ".png"), pp,
         width = ifelse(i == 3, 3, 1.5), height = 2.5, units = "in")
  ggsave(paste0("Boxplot_Control_Conditions_", samp, ".eps"), pp, device = cairo_ps,
         width = ifelse(i == 3, 3, 1.5), height = 2.5, units = "in")

  # Contingency plots
  pp2 <- proportions_sex_col_prop %>%
    ggplot(aes(x = Condition_Group, y = PP * 100, 
               color = Condition_Group, fill = Condition_Group)) +
    geom_col(show.legend = FALSE, width = 0.6) +
    geom_text(data = proportions_sex_col_samp_n, color = "black", inherit.aes = FALSE,
              aes(x = Condition_Group, y = max_y/12, label = N), size = 7/.pt) +
    geom_signif(comparisons = list(c("restraint + SD", "no restraint + HFD")),
                map_signif_level = TRUE,
                annotations = ifelse(test_result1$p.value < 0.001, "***",
                                     ifelse(test_result1$p.value < 0.01, "**",
                                            ifelse(test_result1$p.value < 0.05, "*",
                                                   ifelse(test_result1$p.value < 0.15,
                                                          paste("P =", round(test_result1$p.value, 2)),
                                                          "ns")))),
                size = 0.25,
                tip_length = 0,
                color = "black",
                textsize = ifelse(test_result1$p.value < 0.05, 16/.pt, 7/.pt),
                vjust = ifelse(test_result1$p.value < 0.05, 0.5, 0),
                y_position = max_y,
                fontface = ifelse(test_result1$p.value < 0.05, "bold", "plain")) +
    geom_signif(comparisons = list(c("no restraint + SD", "restraint + HFD")),
                map_signif_level = TRUE,
                annotations = ifelse(test_result2$p.value < 0.001, "***",
                                     ifelse(test_result2$p.value < 0.01, "**",
                                            ifelse(test_result2$p.value < 0.05, "*",
                                                   ifelse(test_result2$p.value < 0.15,
                                                          paste("P =", round(test_result2$p.value, 2)),
                                                          "ns")))),
                size = 0.25,
                tip_length = 0,
                color = "black",
                textsize = ifelse(test_result2$p.value < 0.05, 16/.pt, 7/.pt),
                vjust = ifelse(test_result2$p.value < 0.05, 0.5, 0),
                y_position = max_y*1.1,
                fontface = ifelse(test_result2$p.value < 0.05, "bold", "plain")) +
    geom_signif(comparisons = list(c("restraint + HFD", "no restraint + HFD")),
                map_signif_level = TRUE,
                annotations = ifelse(test_result3$p.value < 0.001, "***",
                                     ifelse(test_result3$p.value < 0.01, "**",
                                            ifelse(test_result3$p.value < 0.05, "*",
                                                   ifelse(test_result3$p.value < 0.15,
                                                          paste("P =", round(test_result3$p.value, 2)),
                                                          "ns")))),
                size = 0.25,
                tip_length = 0,
                color = "black",
                textsize = ifelse(test_result3$p.value < 0.05, 16/.pt, 7/.pt),
                vjust = ifelse(test_result3$p.value < 0.05, 0.5, 0),
                y_position = max_y*1.2,
                fontface = ifelse(test_result3$p.value < 0.05, "bold", "plain")) +
    geom_signif(comparisons = list(c("no restraint + SD", "no restraint + HFD")),
                map_signif_level = TRUE,
                annotations = ifelse(test_result4$p.value < 0.001, "***",
                                     ifelse(test_result4$p.value < 0.01, "**",
                                            ifelse(test_result4$p.value < 0.05, "*",
                                                   ifelse(test_result4$p.value < 0.15,
                                                          paste("P =", round(test_result4$p.value, 2)),
                                                          "ns")))),
                size = 0.25,
                tip_length = 0,
                color = "black",
                textsize = ifelse(test_result4$p.value < 0.05, 16/.pt, 7/.pt),
                vjust = ifelse(test_result4$p.value < 0.05, 0.5, 0),
                y_position = max_y*1.3,
                fontface = ifelse(test_result4$p.value < 0.05, "bold", "plain")) +
    geom_signif(comparisons = list(c("restraint + SD", "restraint + HFD")),
                map_signif_level = TRUE,
                annotations = ifelse(test_result5$p.value < 0.001, "***",
                                     ifelse(test_result5$p.value < 0.01, "**",
                                            ifelse(test_result5$p.value < 0.05, "*",
                                                   ifelse(test_result5$p.value < 0.15,
                                                          paste("P =", round(test_result5$p.value, 2)),
                                                          "ns")))),
                size = 0.25,
                tip_length = 0,
                color = "black",
                textsize = ifelse(test_result5$p.value < 0.05, 16/.pt, 7/.pt),
                vjust = ifelse(test_result5$p.value < 0.05, 0.5, 0),
                y_position = max_y*1.4,
                fontface = ifelse(test_result5$p.value < 0.05, "bold", "plain")) +
    geom_signif(comparisons = list(c("restraint + SD", "no restraint + SD")),
                map_signif_level = TRUE,
                annotations = ifelse(test_result6$p.value < 0.001, "***",
                                     ifelse(test_result6$p.value < 0.01, "**",
                                            ifelse(test_result6$p.value < 0.05, "*",
                                                   ifelse(test_result6$p.value < 0.15,
                                                          paste("P =", round(test_result6$p.value, 2)),
                                                          "ns")))),
                size = 0.25,
                tip_length = 0,
                color = "black",
                textsize = ifelse(test_result6$p.value < 0.05, 16/.pt, 7/.pt),
                vjust = ifelse(test_result6$p.value < 0.05, 0.5, 0),
                y_position = max_y*1.5,
                fontface = ifelse(test_result6$p.value < 0.05, "bold", "plain")) +
    # Custom axis lines
    geom_segment(x = 0.4, y = 0, xend = 0.4, yend = max_y,
                 inherit.aes = FALSE, linewidth = 0.25) + # y-axis
    geom_segment(x = 0.2, y = 0, xend = 5.3, yend = 0,
                 inherit.aes = FALSE, linewidth = 0.25) + # x-axis
    # Y-axis tick marks every 1/4 of the max_y
    geom_segment(x = 0.3, y = max_y/4, xend = 0.4, yend = max_y/4,
                 inherit.aes = FALSE, linewidth = 0.25) +
    geom_segment(x = 0.3, y = max_y/2, xend = 0.4, yend = max_y/2,
                 inherit.aes = FALSE, linewidth = 0.25) +
    geom_segment(x = 0.3, y = 3*max_y/4, xend = 0.4, yend = 3*max_y/4,
                 inherit.aes = FALSE, linewidth = 0.25) +
    geom_segment(x = 0.3, y = max_y, xend = 0.4, yend = max_y,
                 inherit.aes = FALSE, linewidth = 0.25) +
    # X-axis tick marks for conditions (with more spacing)
    geom_segment(x = 1.5, y = -2, xend = 1.5, yend = 0,
                 inherit.aes = FALSE, linewidth = 0.25) + # No restraint center
    geom_segment(x = 4, y = -2, xend = 4, yend = 0,
                 inherit.aes = FALSE, linewidth = 0.25) + # Restraint center
    theme_classic(base_size = 10) +
    scale_x_discrete("",
                     breaks = c("no restraint + SD", "no restraint + HFD",
                                "restraint + SD", "restraint + HFD"),
                     labels = c("No Stress", "", "", " Stress")) + 
    # Add condition labels manually
    scale_y_continuous("Fos+ Cells \\
                       (% TH+ Cells)",
                       limits = c(0, tick_int_small*ceiling(max_y/tick_int_small*1.5)),
                       breaks = seq(0, tick_int_small*ceiling(max_y/tick_int_small*1.5), tick_int_small),
                       labels = tick_y) + # Where the y-tick comes in
    scale_color_manual(values = c("#6f6f6f", "#0698ff", "#6f6f6f", "#0698ff"),
                       labels = Geno_labs) +
    scale_fill_manual(values = c("#dbdbdb", "#79e9fb", "#dbdbdb", "#79e9fb"),
                      labels = Geno_labs) +
    theme(text = element_text(size = 10),
          axis.text = element_text(size = unit(10, "pt"), color = "black"),
          axis.title = element_text(size = 10),
          axis.text.x = element_text(angle = 30, vjust = 0.9, hjust = 1.0, 
                                     size = 10), 
          # axis.text.x = element_blank(), # Remove default x-axis text since we're using custom annotations
          axis.title.y = ggtext::element_markdown(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, margin = margin(0, 0, 10, 0)),
          plot.margin = unit(c(5.5, 5.5, 5.5, 0), "points"),
          legend.position = "none",
          # legend.position = ifelse(i == 3, "right", "none"),
          legend.box.spacing = unit(0, "pt"),
          legend.key.width = unit(0.5, "lines"),
          legend.key.height = unit(0.5, "lines"),
          legend.text = element_text(margin = margin(r = 2)),
          axis.line = element_blank(),
          axis.ticks.length.y = unit(0, "pt"),
          axis.ticks.length.x = unit(0, "pt")) +
    # coord_cartesian(xlim = c(0.2, 5.3), ylim = c(0, 100), 
    #                 expand = FALSE, clip = "off") +
    coord_cartesian(xlim = c(0.2, 5.3), 
                    ylim = c(0, tick_int_small*ceiling(max_y/tick_int_small*1.5)),
                    expand = FALSE, clip = "off") +
    force_panelsizes(rows = unit(1.5, "in"),
                     cols = unit(0.4, "in")) +
    ggtitle(SampleTitle[i])
  
  if (i == 3) {
    pp_legend <- 
      get_legend(ggplot(data.frame(P = rep(0, 4), 
                                   Condition_Group = factor(proportions_sex_col_prop$Condition_Group,
                                                            levels = c("no restraint + SD", "no restraint + HFD", 
                                                                       "restraint + SD", "restraint + HFD"))),
                        aes(x = Condition_Group, y = P, 
                            color = Condition_Group, fill = Condition_Group)) +
                   geom_col(position = "dodge", linewidth = 0.25, show.legend = TRUE) +
                   theme_classic(base_size = 10) +
                   # Add condition labels manually
                   scale_color_manual(values = c("#6f6f6f", "#0698ff", "#6f6f6f", "#0698ff"),
                                      labels = Geno_labs) +
                   scale_fill_manual(values = c("#dbdbdb", "#79e9fb", "#dbdbdb", "#79e9fb"), 
                                     labels = Geno_labs) +    
                   theme(text = element_text(size = 10),
                         axis.title = element_blank(), 
                         plot.margin = unit(c(-5, -5, -5, -5), "points"),
                         legend.title = element_blank(),
                         legend.box.spacing = unit(0, "pt"),
                         legend.key.width = unit(0.5, "lines"),
                         legend.key.height = unit(0.5, "lines"),
                         legend.box.margin = margin(0, 20, 0, -20, unit = "pt"),
                         legend.text = element_text(margin = margin(r = 2))))
    ggsave("Contplot_Control_Conditions_legends.png", 
           pp_legend, 
           width = 2, height = 1, units = "in")
    ggsave("Contplot_Control_Conditions_legends.eps", 
           pp_legend, 
           device = cairo_ps, 
           width = 2, height = 1, units = "in")
  }
  
  ggsave(paste0("Contplot_Control_Conditions_", samp, ".png"), 
         pp2,
         width = 1.5, height = 2.5, units = "in")
  ggsave(paste0("Contplot_Control_Conditions_", samp, ".eps"), 
         pp2, 
         device = cairo_ps, 
         width = 1.5, height = 2.5, units = "in")
  i <- i + 1
}

dev.off()

# Export to Excel with 3 sheets
writexl::write_xlsx(
  list(
    A1 = contrasts.hist_p_A1.df,
    A2 = contrasts.hist_p_A2.df,
    LC = contrasts.hist_p_LC.df
  ),
  path = "contrasts_output.xlsx"
)

