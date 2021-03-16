# Analyzing intracellular pH dynamics in coral cells
# Luella Allen-Waller
# for Innis et al. "Marine heatwaves depress metabolic activity..."

library(dplyr)
library(ggpmisc)
library(ggplot2)
library(doBy)
library(csv)
library(tidyverse)
library(stats)
library(plotrix)
library(arsenal)
library(Rmisc)
library(ggpubr)
library(lme4)
library(lmerTest)

# Check out the timecourse data frame and format
O_frame <- read_csv("Oct19_pHi_timecourse.csv") %>% na.omit()
O_frame$timepoint <- as.numeric(O_frame$timepoint)
head(O_frame)

# read in acidification magnitude and pHi recovery rate data
pH_rates <- read.csv("Oct19_pHi_summary.csv")
pH_rates$celltype <- as.factor(pH_rates$celltype)
pH_rates$hist <- as.factor(pH_rates$hist)
pH_rates$celltype = factor(pH_rates$celltype,levels(pH_rates$celltype)[c(2,1)]) %>%
  revalue(c("nonsymb"="Nonsymbiocytes", "symbiocyte" = "Symbiocytes"))
pH_rates$hist = factor(pH_rates$hist,levels(pH_rates$hist)[c(2,1)]) %>%
  revalue(c("B" = "Bleached", "NB" = "Nonbleached"))

symb_frame <- subset(O_frame, celltype == "symbiocytes")
nonsymb_frame <- subset(O_frame, celltype == "nonsymbiocytes")

####################################################################################

#### Visualize
# Timecourses grouped by history within each cell typ
pHsymbiocytes <- ggplot(symb_frame, aes(x = timepoint - 2, y = pH, color = history, 
                                        group = history)) +
  stat_summary(fun = mean, geom = "path", size = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = "triangle", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", size = 0.8, width=3) +
  ylab(expression(paste("Intracellular pH"))) + xlab("Time post-acidification (min)") +
  scale_color_manual(values = c("Bleached" = "gray", "Nonbleached" = "black")) +
  theme_bw() + ggtitle("Symbiocytes") +
  theme(panel.grid = element_blank(), legend.position = "none", text = element_text(size=16)) +
  coord_cartesian(ylim = c(7.0,8.1))
pHsymbiocytes

pHnonsymbiocytes <- ggplot(nonsymb_frame, aes(x = timepoint - 2, y = pH, color = history,
                                              group = history)) +
  stat_summary(fun = mean, geom = "path", size = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = "triangle", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", size = 0.8, width=3) +
  ylab(expression(paste("Intracellular pH"))) + xlab("Time post-acidification (min)") +
  scale_color_manual(values = c("Bleached" = "gray", "Nonbleached" = "black")) +
  theme_bw() + ggtitle("Nonsymbiocytes") + labs(color = "Bleaching history") +
  theme(panel.grid = element_blank(), legend.position = c(0.75,0.75), legend.text.align = 1, 
        legend.title.align = 1, text = element_text(size=16)) +
  coord_cartesian(ylim = c(7.0,8.1))
pHnonsymbiocytes

# setpoint pH (basal) boxplot
basal_frame <- filter(O_frame, timepoint < 5)
basal_frame$history <- as.factor(basal_frame$history)
basal_frame$celltype <- factor(basal_frame$celltype, levels = basal_frame$celltype[c(2,1)]) %>%
  revalue(c("nonsymbiocytes"="Nonsymbiocytes", "symbiocytes" = "Symbiocytes"))
basal_frame$history <- factor(basal_frame$history, levels = c("Nonbleached", "Bleached"))
pHBasal <- ggplot(basal_frame, aes(x = celltype, y = pH, color = history,
                                   group = interaction(celltype, history))) +
  geom_boxplot(size = 0.7) +
  ylab(expression(paste("Basal pH"[i]))) + 
  xlab("") + 
  scale_color_manual(values = c("Nonbleached" = "black", "Bleached" = "gray")) +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "none", text = element_text(size=14))
pHBasal

# acidification boxplot
pHAcidification <- ggplot(pH_rates, aes(x = celltype, y = acid_mag, color = hist,
                                        group = interaction(celltype, hist))) +
  geom_boxplot(size = 0.7) +
  ylab(expression(paste("Initial acidification ("* Delta*"pH"[i]*")"))) + 
  xlab("") +
  scale_color_manual(values = c("Nonbleached" = "black", "Bleached" = "gray")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), legend.position = "none", text = element_text(size = 14)) +
  geom_hline(yintercept = 0, linetype = "dashed")
pHAcidification

# recovery rate boxplot
pHRecovery <- ggplot(pH_rates, aes(x = celltype, y = recovery_rate, color = hist,
                                   group = interaction(celltype, hist))) +
  geom_boxplot(size = 0.7) +
  ylab(expression(paste("pH"[i]*" recovery rate (pH min"^-1*")"))) + xlab("") +
  scale_color_manual(values = c("Nonbleached" = "black", "Bleached" = "gray")) +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "none", text = element_text(size = 14)) +
  geom_hline(yintercept = 0, linetype = "dashed")
pHRecovery

### July controls
control_rates <- read.csv("Jul19control_pHi_summary.csv")

# setpoint pH (basal) boxplot
c_pHBasal <- ggplot(control_rates, aes(x = celltype, y = basal, color = hist,
                                   group = interaction(celltype, hist))) +
  geom_boxplot(size = 0.7) +
  ylab(expression(paste("Basal pH"[i]))) + 
  xlab("") + 
  scale_color_manual(values = c("Resistant" = "black", "Susceptible" = "gray")) +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "none", text = element_text(size=14))
c_pHBasal

# acidification boxplot
c_pHAcidification <- ggplot(control_rates, aes(x = celltype, y = acid_mag, color = hist,
                                        group = interaction(celltype, hist))) +
  geom_boxplot(size = 0.7) +
  ylab(expression(paste("Initial acidification ("* Delta*"pH"[i]*")"))) + 
  xlab("") +
  scale_color_manual(values = c("Resistant" = "black", "Susceptible" = "gray")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), legend.position = "none", text = element_text(size = 14)) +
  geom_hline(yintercept = 0, linetype = "dashed")
c_pHAcidification

# recovery rate boxplot
c_pHRecovery <- ggplot(control_rates, aes(x = celltype, y = recovery_rate, color = hist,
                                   group = interaction(celltype, hist))) +
  geom_boxplot(size = 0.7) +
  ylab(expression(paste("pH"[i]*" recovery rate (pH min"^-1*")"))) + xlab("") +
  scale_color_manual(values = c("Resistant" = "black", "Susceptible" = "gray")) +
  theme_bw() + 
  theme(panel.grid = element_blank(), legend.position = "none", text = element_text(size = 14)) +
  geom_hline(yintercept = 0, linetype = "dashed")
c_pHRecovery

####################################
### Stats 

# Setpoint/basal pH
# Peak bleaching (October) experimentals:
# Check normality: 
basal_B_symb <- subset(basal_frame, history == "B" & celltype == "symbiocyte")
shapiro.test(basal_B_symb$pH)
basal_NB_symb <- subset(basal_frame, history == "NB" & celltype == "symbiocyte")
shapiro.test(basal_NB_symb$pH)
basal_B_nonsymb <- subset(basal_frame, history == "B" & celltype == "nonsymbiocyte")
shapiro.test(basal_B_nonsymb$pH)
basal_NB_nonsymb <- subset(basal_frame, history == "NB" & celltype == "nonsymbiocyte")
shapiro.test(basal_NB_nonsymb$pH)
# anova:
basal.aov <- aov(pH ~ history * celltype, data = basal_frame)
summary(basal.aov)
# post-hoc testing:
TukeyHSD(basal.aov, ordered = T, conf.level = 0.99)

# July controls:
# Check normality: 
c_basal_B_symb <- subset(control_rates, hist == "Susceptible" & celltype == "symbiocytes")
shapiro.test(c_basal_B_symb$basal)
# Not normal
basal_NB_symb <- subset(control_rates, hist == "Resistant" & celltype == "symbiocytes")
shapiro.test(basal_NB_symb$basal)
basal_B_nonsymb <- subset(control_rates, hist == "Susceptible" & celltype == "nonsymbiocytes")
shapiro.test(basal_B_nonsymb$basal)
basal_NB_nonsymb <- subset(control_rates, hist == "Resistant" & celltype == "nonsymbiocytes")
shapiro.test(basal_NB_nonsymb$basal)
# Kruskal-Wallis
kruskal.test(basal ~ hist, data = control_rates)
# No significant effect of bleaching history (X2 = 0.45833, df = 1, p-value = 0.4984)
kruskal.test(basal ~ celltype, data = control_rates)
# No significant effect of celltype (X2 = 3.6402, df = 1, p-value = 0.0564)

#### Linear mixed effects models on whole timecourses
# H/T Craig Nelson, Kristin Brown, and Teegan Innis for advice

### Peak bleaching (October 2019)
qqnorm(O_frame$pH)

# LMM on all cells together, with time (numeric), cell type and history as fixed effects
# and pair as a random effect
pair.lmm.timecourse <- lmer(pH ~ timepoint * celltype * history + (1|pair), data = O_frame)
summary(pair.lmm.timecourse)
anova(pair.lmm.timecourse, type="III")

# LMM on acidification magnitude, with fixed effects of bleaching history and cell type
pair.lmm.acidmag <- lmer(acid_mag ~  celltype * hist + (1|pair), data = pH_rates)
anova(pair.lmm.acidmag, type="III")

# for basal pHi: take only 0 timepoint, see whether cell type etc. affects
basal_frame <- filter(O_frame, timepoint < 3)
pair.lmm.basal <- lmer(pH ~ celltype * history + (1|pair), data = basal_frame)
anova(pair.lmm.basal, type="III")

# for recovery rate
pair.lmm.recovrate <- lmer(recovery_rate ~  celltype * hist + (1|pair), data = pH_rates)
anova(pair.lmm.recovrate, type="III")


### Controls (July 2019)
shapiro.test(control_rates$basal)
qqnorm(control_rates$basal)
# p = 0.013 - data are non-normal...slightly
# for basal pHi:
c.lmm.basal <- lmer(basal ~ celltype * hist + (1|pair), data = control_rates)
anova(c.lmm.basal, type="III")

# LMM on acidification magnitude, with fixed effects of bleaching history and cell type
shapiro.test(control_rates$acid_mag) # normal
c.lmm.acidmag <- lmer(acid_mag ~  celltype * hist + (1|pair), data = control_rates)
anova(c.lmm.acidmag, type="III")

# for recovery rate
shapiro.test(control_rates$recovery_rate) # normal
c.lmm.recovrate <- lmer(recovery_rate ~  celltype * hist + (1|pair), data = control_rates)
anova(c.lmm.recovrate, type="III")
