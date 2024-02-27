# Alpha diversity analysis

- [A) Shannon Diversity](#a--shannon-diversity)
  * [1. Load libraries](#1-load-libraries)
  * [2. Load the data](#2-load-the-data)
  * [3. Plot shannon distribution](#3-plot-shannon-distribution)
  * [4. Model  sex + location](#4-model--sex---location)
  * [5. Model SMI + location](#5-model-smi---location)
  * [6. Model shannon + infection status](#6-model-shannon---infection-status)
  * [7. Model infection + location](#7-model-infection---location)
  * [8. Adjust p values BH correction](#8-adjust-p-values-bh-correction)
  * [9. Model interaction SMI*location](#9-model-interaction-smi-location)
  * [10. Model interaction location*infection](#10-model-interaction-location-infection)
- [B) Faith PD](#b--faith-pd)
  * [1. Plot faith distribution](#1-plot-faith-distribution)
  * [2. Model  sex + location](#2-model--sex---location)
  * [3. Model SMI + location](#3-model-smi---location)
  * [4. Model faith + infection status](#4-model-faith---infection-status)
  * [5. Model infection + location](#5-model-infection---location)
  * [6. Adjust p values BH correction](#6-adjust-p-values-bh-correction)
  * [7. Model interaction SMI*location](#7-model-interaction-smi-location)
  * [8. Model interaction location*infection](#8-model-interaction-location-infection)

# A) Shannon Diversity



## 1. Load libraries

```R
library(qiime2R)
library(lubridate)
library(phyloseq)
library(tidyverse)
library(microbiome)
library(lme4)
library(MuMIn)
library (performance)
library(datawizard)
library(car)
library(effects)
library(openxlsx) 
```



## 2. Load the data

```R
metadata <- read.xlsx("sparrows_meta.xlsx", detectDates = T)

str(metadata)
metadata <- janitor::clean_names(metadata)
metadata <- clean_names(metadata)
metadata$location <- gsub("IBISS, yard", "IBISS", metadata$location)
metadata$location <- gsub("Kalemegdan, ZOO", "ZOO", metadata$location)
metadata$identifier <- as.factor(metadata$identifier)
metadata$ring_number <- as.factor(metadata$ring_number)
metadata$location <- as.factor(metadata$location)
metadata$sex <- as.factor(metadata$sex)
metadata$m_infection <- as.factor(metadata$m_infection)
metadata$SMI <- as.numeric(metadata$SMI)
metadata$std_SMI <- scale(metadata$SMI) # scale body mass idex
metadata$std_SMI <- as.numeric(metadata$std_SMI)
metadata$parasetimia <- as.numeric(metadata$parasetimia)

saveRDS(metadata, "metadata_sparrows.rds")
metadata <- readRDS("metadata_sparrows.rds")

```



## 3. Plot shannon distribution

```R
hist(metadata$shannon_entropy)

# test for normal distribution
shapiro.test(metadata$shannon_entropy)

#Transform shannon (after running all models with non transformed shannon, model diagnostics showed that residuals were not normal distributed)
metadata$log_shannon <- log(metadata$shannon_entropy)

#Check for outliers shannon
check_outliers(metadata$log_shannon)

plot(check_distribution(metadata$log_shannon)
```



## 4. Model  sex + location

```R
model_shannon_sex <- lm(log_shannon ~ location + sex, data = metadata)

# model diagnostics
plot(check_distribution(model_shannon_m_infection))#distribution of response seems to be gamma but normal is also a good fit
check_normality(model_shannon_sex)
check_model(model_shannon_sex)
summary(model_shannon_sex)

#get p values and R2
p_sex_shannon <- Anova(model_shannon_sex) 
r.squaredGLMM(model_shannon_sex)
plot(allEffects(model_shannon_sex)) 

#plot shannon vs sex
shannon_sex_plot <- ggplot(metadata, aes(x = sex, y = log_shannon)) +
  geom_boxplot(fill = "white", width=0.5) +  # No fill color for the boxplot
  geom_jitter(position = position_jitter(0.2), aes(color = location), alpha = 1, size=2.5) +
  #stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.03) +
  labs(x = "Sex", y = "Log(Shannon diversity)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  # Adjust the size as needed
        axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18),  # Adjust the size as needed
        axis.title.y = element_text(size = 18))+
  theme( panel.grid.major = element_blank(),  # Remove major grid lines
         panel.grid.minor = element_blank())+ # Remove minor grid lines
  theme(text = element_text(family = "Arial"))+
  scale_color_manual(values = c("#4CAF50", "#D32F2F"))+
  theme(legend.position = "none")
shannon_sex_plot


```



## 5. Model SMI + location

```R
model_shannon_SMI <- lm(log_shannon ~ location + std_SMI, data = metadata)

#Model diagnostics
plot(check_distribution(model_shannon_SMI))
check_normality(model_shannon_SMI)
check_model(model_shannon_SMI)
summary(model_shannon_SMI)

#Get p values and R2
p_SMI_shannon <- Anova(model_shannon_SMI) #
r.squaredGLMM(model_shannon_SMI)

# plot shannon vs location
shannon_location_plot <- ggplot(metadata, aes(x = location, y = log_shannon)) +
  geom_boxplot(fill = "white", width=0.5) +  # No fill color for the boxplot
  geom_jitter(position = position_jitter(0.2), aes(color = location), alpha = 1, size=2.5) +
  #stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.03) +
  labs(x = "Sampling location", y = "Log(Shannon index)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  # Adjust the size as needed
        axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18),  # Adjust the size as needed
        axis.title.y = element_text(size = 18))+
  theme( panel.grid.major = element_blank(),  # Remove major grid lines
         panel.grid.minor = element_blank())+ # Remove minor grid lines
  theme(text = element_text(family = "Arial"))+
  scale_color_manual(values = c("#4CAF50", "#D32F2F"))+
  theme(legend.position = "none")
shannon_location_plot

#plot shannon vs SMI
shannon_SMI_plot <- ggplot(metadata, aes(x = std_SMI, y = log_shannon)) +
  geom_point(position = position_jitter(0.2), aes(color = location), alpha = 1, size=2.5) +
  geom_smooth(se = TRUE, method = 'lm') +
  labs(x = " Std. Scale mass index", y = "Log(Shannon diversity)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  # Adjust the size as needed
        axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18),  # Adjust the size as needed
        axis.title.y = element_text(size = 18))+
  theme( panel.grid.major = element_blank(),  # Remove major grid lines
         panel.grid.minor = element_blank())+ # Remove minor grid lines
  theme(text = element_text(family = "Arial"))+
  scale_color_manual(values = c("#4CAF50", "#D32F2F"))+
  theme(legend.position = "none")
shannon_SMI_plot

plot(predictorEffect("std_SMI", model_shannon_SMI, residuals=TRUE),partial.residuals=list(smooth=FALSE))
```



## 6. Model shannon + infection status

```R
model_shannon_infection <- lm(log_shannon ~ location + m_infection ,data = metadata)

#Model diagnostics
plot(check_distribution(model_shannon_infection))
check_normality(model_shannon_infection)
check_model(model_shannon_infection)
summary(model_shannon_infection)

#Get p values and R2
p_infection_shannon <- Anova(model_shannon_infection) # p> 0.05 no direct effect of body condition
r.squaredGLMM(model_shannon_infection)

plot(allEffects(model_shannon_infection))

#plot shannon vs infection status

shannon_infection_plot <- ggplot(metadata, aes(x = m_infection, y = log_shannon)) +
  geom_boxplot(fill = "white", width=0.5) +  # No fill color for the boxplot
  geom_jitter(position = position_jitter(0.2), aes(color = location), alpha = 1, size=2.5) +
  #stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.03) +
  labs(x = "Infection status", y = "Log(Shannon diversity)") +
  scale_x_discrete(labels = c("0" = "Non-infected", "1" = "Infected")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  # Adjust the size as needed
        axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18),  # Adjust the size as needed
        axis.title.y = element_text(size = 18))+
  theme( panel.grid.major = element_blank(),  # Remove major grid lines
         panel.grid.minor = element_blank())+ # Remove minor grid lines
  theme(text = element_text(family = "Arial"))+
  scale_color_manual(values = c("#4CAF50", "#D32F2F")) + 
  guides(fill = guide_legend(title = "location"))
shannon_infection_plot
```



## 7. Model infection + location

```R
model_shannon_infection <- lm(log_shannon ~ location + m_infection ,data = metadata)

#Model diagnostics
plot(check_distribution(model_shannon_infection))
check_normality(model_shannon_infection)
check_model(model_shannon_infection)
summary(model_shannon_infection)

#Get p values and R2
p_infection_shannon <- Anova(model_shannon_infection) # p> 0.05 no direct effect of body condition
r.squaredGLMM(model_shannon_infection)

#plot shannon vs infection status
shannon_infection_plot <- ggplot(metadata, aes(x = m_infection, y = log_shannon)) +
  geom_boxplot(fill = "white", width=0.5) +  # No fill color for the boxplot
  geom_jitter(position = position_jitter(0.2), aes(color = location), alpha = 1, size=2.5) +
  #stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.03) +
  labs(x = "Infection status", y = "Log(Shannon diversity)") +
  scale_x_discrete(labels = c("0" = "Non-infected", "1" = "Infected")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  # Adjust the size as needed
        axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18),  # Adjust the size as needed
        axis.title.y = element_text(size = 18))+
  theme( panel.grid.major = element_blank(),  # Remove major grid lines
         panel.grid.minor = element_blank())+ # Remove minor grid lines
  theme(text = element_text(family = "Arial"))+
  scale_color_manual(values = c("#4CAF50", "#D32F2F")) + 
  guides(fill = guide_legend(title = "location"))
shannon_infection_plot
```



## 8. Adjust p values BH correction

```R
p_values_shannon <- c(p_sex_shannon$`Pr(>F)`[1], p_SMI_shannon$`Pr(>F)`[1], p_infection_shannon$`Pr(>F)`[1])

adjusted_p_values_shannon <- p.adjust(p_values_shannon, method = "BH")
adjusted_p_values_shannon <- format(round(adjusted_p_values_shannon, digits = 3), scientific = FALSE)
adjusted_p_values_shannon

p_values_shannon2 <- c(p_sex_shannon$`Pr(>F)`[2], p_SMI_shannon$`Pr(>F)`[2], p_infection_shannon$`Pr(>F)`[2])

adjusted_p_values_shannon2 <- p.adjust(p_values_shannon2, method = "BH")
adjusted_p_values_shannon2 <- format(round(adjusted_p_values_shannon2, digits = 3), scientific = FALSE)
adjusted_p_values_shannon2
```



## 9. Model interaction SMI*location

```R
model_interaction <- lm(log_shannon ~ location*std_SMI, data = metadata)

#Model diagnostics
plot(check_distribution(model_interaction))
check_normality(model_interaction)
check_model(model_interaction)
summary(model_interaction)

#Get p values and R2
p_interaction <- Anova(model_interaction) # p> 0.05 no direct effect of body condition only location
r.squaredGLMM(model_interaction)

#plot interaction
plot(allEffects(model_interaction, residuals=TRUE),partial.residuals=list(smooth=FALSE))

effect_interaction<- plot (predictorEffect("std_SMI", model_interaction, residuals=TRUE),partial.residuals=list(smooth=FALSE,col="black"),
                           xlab = "Std SMI", ylab = "Log (Shannon)", main = "")
par(family = "Arial", cex.lab = 1.4, cex.axis = 1.4, cex.main = 1.4, cex.sub = 1.4)

effect_interaction 
```



## 10. Model interaction location*infection

```R
model_interaction2 <- lm(log_shannon ~ location*m_infection, data = metadata)

#Model daignostics
plot(check_distribution(model_interaction2))
check_normality(model_interaction2)
check_model(model_interaction2)
summary(model_interaction2)

#Get p values and R2
p_interaction2 <- Anova(model_interaction2) # p> 0.05 no direct effect of body condition only location
r.squaredGLMM(model_interaction2)

#plot interaction 
plot(allEffects(model_interaction, residuals=TRUE),partial.residuals=list(smooth=FALSE))

effect_interaction2<- plot (predictorEffect("std_SMI", model_interaction, residuals=TRUE),partial.residuals=list(smooth=FALSE,col="black"),
                            xlab = "Std SMI", ylab = "Log (Shannon)", main = "")
par(family = "Arial", cex.lab = 1.4, cex.axis = 1.4, cex.main = 1.4, cex.sub = 1.4)
effect_interaction2 
```





# B) Faith PD



## 1. Plot faith distribution

```R
hist(metadata$faith_pd)

# test for normal distribution
shapiro.test(metadata$faith_pd)
plot(check_distribution(metadata$faith_pd))

##Check for outliers
check_outliers(metadata$faith_pd)
metadata_out <- metadata[-c(14),] #remove outliers

#Transform faith (after running all models with non transformed shannon, model diagnostics showed that residuals were not normal distributed)
metadata$log_faith <- log(metadata$faith_pd)
plot(check_distribution(metadata$log_faith))
shapiro.test(metadata$log_faith)
```



## 2. Model  sex + location

```R
model_faith_sex <- lm(log_faith ~ location + sex, data = metadata)

# model diagnostics
plot(check_distribution(model_faith_m_infection))#distribution of response seems to be gamma but normal is also a good fit
check_normality(model_faith_sex)
check_model(model_faith_sex)
summary(model_faith_sex)

#get p values and R2
p_sex_faith <- Anova(model_faith_sex) 
r.squaredGLMM(model_faith_sex)
plot(allEffects(model_faith_sex)) 

#plot faith vs sex
faith_sex_plot <- ggplot(metadata, aes(x = sex, y = log_faith)) +
  geom_boxplot(fill = "white", width=0.5) +  # No fill color for the boxplot
  geom_jitter(position = position_jitter(0.2), aes(color = location), alpha = 1, size=2.5) +
  #stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.03) +
  labs(x = "Sex", y = "Log(faith diversity)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  # Adjust the size as needed
        axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18),  # Adjust the size as needed
        axis.title.y = element_text(size = 18))+
  theme( panel.grid.major = element_blank(),  # Remove major grid lines
         panel.grid.minor = element_blank())+ # Remove minor grid lines
  theme(text = element_text(family = "Arial"))+
  scale_color_manual(values = c("#4CAF50", "#D32F2F"))+
  theme(legend.position = "none")
faith_sex_plot


```



## 3. Model SMI + location

```R
model_shannon_SMI <- lm(log_shannon ~ location + std_SMI, data = metadata)

#Model diagnostics
plot(check_distribution(model_shannon_SMI))
check_normality(model_shannon_SMI)
check_model(model_shannon_SMI)
summary(model_shannon_SMI)

#Get p values and R2
p_SMI_shannon <- Anova(model_shannon_SMI) #
r.squaredGLMM(model_shannon_SMI)

# plot shannon vs location
shannon_location_plot <- ggplot(metadata, aes(x = location, y = log_shannon)) +
  geom_boxplot(fill = "white", width=0.5) +  # No fill color for the boxplot
  geom_jitter(position = position_jitter(0.2), aes(color = location), alpha = 1, size=2.5) +
  #stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.03) +
  labs(x = "Sampling location", y = "Log(Shannon index)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  # Adjust the size as needed
        axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18),  # Adjust the size as needed
        axis.title.y = element_text(size = 18))+
  theme( panel.grid.major = element_blank(),  # Remove major grid lines
         panel.grid.minor = element_blank())+ # Remove minor grid lines
  theme(text = element_text(family = "Arial"))+
  scale_color_manual(values = c("#4CAF50", "#D32F2F"))+
  theme(legend.position = "none")
shannon_location_plot

#plot shannon vs SMI
shannon_SMI_plot <- ggplot(metadata, aes(x = std_SMI, y = log_shannon)) +
  geom_point(position = position_jitter(0.2), aes(color = location), alpha = 1, size=2.5) +
  geom_smooth(se = TRUE, method = 'lm') +
  labs(x = " Std. Scale mass index", y = "Log(Shannon diversity)") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  # Adjust the size as needed
        axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18),  # Adjust the size as needed
        axis.title.y = element_text(size = 18))+
  theme( panel.grid.major = element_blank(),  # Remove major grid lines
         panel.grid.minor = element_blank())+ # Remove minor grid lines
  theme(text = element_text(family = "Arial"))+
  scale_color_manual(values = c("#4CAF50", "#D32F2F"))+
  theme(legend.position = "none")
shannon_SMI_plot

plot(predictorEffect("std_SMI", model_shannon_SMI, residuals=TRUE),partial.residuals=list(smooth=FALSE))
```



## 4. Model faith + infection status

```R
model_faith_infection <- lm(log_faith ~ location + m_infection ,data = metadata)

#Model diagnostics
plot(check_distribution(model_faith_infection))
check_normality(model_faith_infection)
check_model(model_faith_infection)
summary(model_faith_infection)

#Get p values and R2
p_infection_faith <- Anova(model_faith_infection) # p> 0.05 no direct effect of body condition
r.squaredGLMM(model_faith_infection)

plot(allEffects(model_faith_infection))

#plot faith vs infection status

faith_infection_plot <- ggplot(metadata, aes(x = m_infection, y = log_faith)) +
  geom_boxplot(fill = "white", width=0.5) +  # No fill color for the boxplot
  geom_jitter(position = position_jitter(0.2), aes(color = location), alpha = 1, size=2.5) +
  #stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.03) +
  labs(x = "Infection status", y = "Log(faith diversity)") +
  scale_x_discrete(labels = c("0" = "Non-infected", "1" = "Infected")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  # Adjust the size as needed
        axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18),  # Adjust the size as needed
        axis.title.y = element_text(size = 18))+
  theme( panel.grid.major = element_blank(),  # Remove major grid lines
         panel.grid.minor = element_blank())+ # Remove minor grid lines
  theme(text = element_text(family = "Arial"))+
  scale_color_manual(values = c("#4CAF50", "#D32F2F")) + 
  guides(fill = guide_legend(title = "location"))
faith_infection_plot
```



## 5. Model infection + location

```R
model_faith_infection <- lm(log_faith ~ location + m_infection ,data = metadata)

#Model diagnostics
plot(check_distribution(model_faith_infection))
check_normality(model_faith_infection)
check_model(model_faith_infection)
summary(model_faith_infection)

#Get p values and R2
p_infection_faith <- Anova(model_faith_infection) # p> 0.05 no direct effect of body condition
r.squaredGLMM(model_faith_infection)

#plot faith vs infection status
faith_infection_plot <- ggplot(metadata, aes(x = m_infection, y = log_faith)) +
  geom_boxplot(fill = "white", width=0.5) +  # No fill color for the boxplot
  geom_jitter(position = position_jitter(0.2), aes(color = location), alpha = 1, size=2.5) +
  #stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.03) +
  labs(x = "Infection status", y = "Log(faith diversity)") +
  scale_x_discrete(labels = c("0" = "Non-infected", "1" = "Infected")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14),  # Adjust the size as needed
        axis.text.y = element_text(size = 14))+
  theme(axis.title.x = element_text(size = 18),  # Adjust the size as needed
        axis.title.y = element_text(size = 18))+
  theme( panel.grid.major = element_blank(),  # Remove major grid lines
         panel.grid.minor = element_blank())+ # Remove minor grid lines
  theme(text = element_text(family = "Arial"))+
  scale_color_manual(values = c("#4CAF50", "#D32F2F")) + 
  guides(fill = guide_legend(title = "location"))
faith_infection_plot
```



## 6. Adjust p values BH correction

```
p_values_faith <- c(p_sex_faith$`Pr(>F)`[1], p_SMI_faith$`Pr(>F)`[1], p_infection_faith$`Pr(>F)`[1])

adjusted_p_values_faith <- p.adjust(p_values_faith, method = "BH")
adjusted_p_values_faith <- format(round(adjusted_p_values_faith, digits = 3), scientific = FALSE)
adjusted_p_values_faith

p_values_faith2 <- c(p_sex_faith$`Pr(>F)`[2], p_SMI_faith$`Pr(>F)`[2], p_infection_faith$`Pr(>F)`[2])

adjusted_p_values_faith2 <- p.adjust(p_values_faith2, method = "BH")
adjusted_p_values_faith2 <- format(round(adjusted_p_values_faith2, digits = 3), scientific = FALSE)
adjusted_p_values_faith2
```



## 7. Model interaction SMI*location

```R
model_interaction <- lm(log_faith ~ location*std_SMI, data = metadata)

#Model diagnostics
plot(check_distribution(model_interaction))
check_normality(model_interaction)
check_model(model_interaction)
summary(model_interaction)

#Get p values and R2
p_interaction <- Anova(model_interaction) # p> 0.05 no direct effect of body condition only location
r.squaredGLMM(model_interaction)

#plot interaction
plot(allEffects(model_interaction, residuals=TRUE),partial.residuals=list(smooth=FALSE))

effect_interaction<- plot (predictorEffect("std_SMI", model_interaction, residuals=TRUE),partial.residuals=list(smooth=FALSE,col="black"),
                           xlab = "Std SMI", ylab = "Log (faith)", main = "")
par(family = "Arial", cex.lab = 1.4, cex.axis = 1.4, cex.main = 1.4, cex.sub = 1.4)

effect_interaction 
```



## 8. Model interaction location*infection

```R
model_interaction2 <- lm(log_faith ~ location*m_infection, data = metadata)

#Model daignostics
plot(check_distribution(model_interaction2))
check_normality(model_interaction2)
check_model(model_interaction2)
summary(model_interaction2)

#Get p values and R2
p_interaction2 <- Anova(model_interaction2) # p> 0.05 no direct effect of body condition only location
r.squaredGLMM(model_interaction2)

#plot interaction 
plot(allEffects(model_interaction, residuals=TRUE),partial.residuals=list(smooth=FALSE))

effect_interaction2<- plot (predictorEffect("std_SMI", model_interaction, residuals=TRUE),partial.residuals=list(smooth=FALSE,col="black"),
                            xlab = "Std SMI", ylab = "Log (faith)", main = "")
par(family = "Arial", cex.lab = 1.4, cex.axis = 1.4, cex.main = 1.4, cex.sub = 1.4)
effect_interaction2
```

