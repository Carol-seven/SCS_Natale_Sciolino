library(ggplot2)
library(ggeffects)
library(glmmTMB)
library(gtsummary)
library(tidyr)
library(ggplot2)
library(sjPlot)
library(readxl)
library(janitor)
library(marginaleffects)
library(emmeans)
library(performance)
list.files()
df <- read_xlsx("Cricket Hunting Task Data 10-10-2025.xlsx")
df <- clean_names(df)
head(df)
str(df)
colnames(df)
table(df$stim_frequency)
df$stim_frequency <- as.character(df$stim_frequency)
table(df$gender)
table(df$ms_number)
df$grp <- df$group
table(df$grp)
df <- df %>% dplyr::select(-group)

table(df$grp, df$eaten_binary, df$stim_frequency)

df$killed_binary



# Full Model with  no interactions
m1_full <- glmmTMB(eaten_binary ~ gender + grp +
                            stim_frequency + implant_sites + (1|ms_number), 
              data = df, family = binomial())
# model summary
summary(m1_full)
tab_model(m1_full)
check_model(m1_full)
check_residuals(m1_full)

# marginal means and effects. NOTE: these are on the probability scale, so for this model they are the Pr(eaten)
# marginal means group
avg_predictions(m1_full,variables = "grp", re.form = NULL) # this accounts for RE variance

# marginal means group
avg_comparisons(m1_full,variables = "grp", re.form = NULL) # this accounts for RE variance

# marginal means sex
avg_predictions(m1_full,variables = "gender", re.form = NULL) # this accounts for RE variance

# marginal means sex
avg_comparisons(m1_full,variables = "gender", re.form = NULL) # this accounts for RE variance

# marginal means stim_frequency
avg_predictions(m1_full,variables = "stim_frequency", re.form = NULL) # this accounts for RE variance

# marginal means sex
avg_comparisons(m1_full,variables = list(stim_frequency = "pairwise"), re.form = NULL) # this accounts for RE variance


# Reduced Model with  interaction
m1_int <- glmmTMB(eaten_binary ~ grp*stim_frequency + implant_sites + (1|ms_number), 
                   data = df, family = binomial())
# model summary
summary(m1_int)
tab_model(m1_int)
check_model(m1_int)
check_residuals(m1_int)

newdata <- datagrid(model=m1_int, grp = c("EYFP", "ChR2"), stim_frequency = c("5", "10", "20"))

# marginal means and effects. NOTE: these are on the probability scale, so for this model they are the Pr(eaten)
# marginal means group by stim_frequency
avg_predictions(m1_int, newdata = newdata, by = c("grp", "stim_frequency"), re.form = NULL)

# marginal effects group by stim_frequency
avg_comparisons(m1_int, newdata = newdata, variables = list(grp = "pairwise"), by = c("stim_frequency"), re.form = NULL)

# plotting marginal effects
plot_comparisons(m1_int, newdata = newdata, variables = list(grp = "pairwise"), by = c("stim_frequency"), re.form = NULL, 
                 vcov = TRUE)+ scale_y_continuous(limits = c(-1,0.1))+ geom_hline(yintercept = 0, lty = 2)+
  theme_bw()+
  labs(x = "Stimulus Frequency",
       y = "Marginal Effect (EYFP - ChR2)")
ggsave("comparisons.png", height = 6, width = 6.5, units = "in", dpi = 400)


## additional models
table(df$implant_sites)


m2_BNST <- glmmTMB(eaten_binary ~ grp*stim_frequency + (1|ms_number), 
                  data = df %>% filter(implant_sites=="BNST"), family = binomial())
summary(m2_BNST)
tab_model(m2_BNST)
check_model(m2_BNST)
check_residuals(m2_BNST)


m2_DP <- glmmTMB(eaten_binary ~ grp*stim_frequency + (1|ms_number), 
                  data = df %>% filter(implant_sites=="DP"), family = binomial())
summary(m2_DP)
tab_model(m2_DP)
check_model(m2_DP)
check_residuals(m2_DP)

m2_LHA <- glmmTMB(eaten_binary ~ grp*stim_frequency + (1|ms_number), 
                  data = df %>% filter(implant_sites=="LHA"), family = binomial())
summary(m2_LHA)
tab_model(m2_LHA)
check_model(m2_LHA)
check_residuals(m2_LHA)
