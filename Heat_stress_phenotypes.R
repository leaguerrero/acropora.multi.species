# Testing the correlation between heat stress data
# - Color analysis
# - Symbiont counts (normalized by protein content)
# - ED50

library(ggplot2)
library(tidyverse)
library(plotly)
library(ggbeeswarm)
library(beeswarm)
library(cowplot)
library(plyr)
library(performance)
library(qqplotr)
library(ggpubr)

# Official species colors
official.species.color.palette<- c("#F7C11E", # Acropora tutuilensis
                             "#FF9A17", # Acropora lutkeni
                             "#fb5607", # Acropora retusa
                             "#ff006e", # Acropora pulchra
                             "#c11cad", # Acrpora spp
                             "#8338ec", # Acropora robusta
                             "#3E6FCB", # Acropora hyacinthus
                             "#3a86ff") # Acropora abrotanoides 


# Load the symbiont count data
cyt_data <- read.csv("~/Lab Notebook/Chapter3/Thermal_stress_analysis/Symbiont_counts/Data/moorea 2021 flow cytometry data.xlsx - data.csv")
names(cyt_data)
dim(cyt_data) # 228  25
cyt_data_trim <- cyt_data[,-25]

# Update species ID based on Tag-Seq Phylogeny
my.metadata <- read.csv("/Users/tillandsia/Lab Notebook/Chapter3/Tag_seq_analysis/tagseq_meta_data.csv", header = TRUE)
cyt_data$species<-my.metadata[match(cyt_data$colony, my.metadata$colony.id),2]
## Remove data rows where species ID is NA
sum(is.na(cyt_data$species)) # 12
cyt_data_clean_species <- cyt_data[complete.cases(cyt_data[ , 2]),]
sum(is.na(cyt_data_clean_species$species)) 

# Add proportion of symbiont:coral cells (P1:P3)
cyt_data_clean_species$proportion.p1 <- cyt_data_clean_species$events.p1/cyt_data_clean_species$events.p3

## Summary statistics of proportion.p1
# Explain the central tendency and dispersion of this column
mean(cyt_data_clean_species$proportion.p1) # 0.4739381
median(cyt_data_clean_species$proportion.p1) # 0.3863381
sd(cyt_data_clean_species$proportion.p1) #0.5245752
min(cyt_data_clean_species$proportion.p1) # 0.02968878
max(cyt_data_clean_species$proportion.p1) # 4.815854

## Visually explore the proportion.p1 variable:
ggplot(cyt_data_clean_species, aes(x = proportion.p1)) + 
  theme_classic() + 
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(aes(xintercept = mean(proportion.p1)), color = "red", linetype = "dashed") +
  geom_vline(aes(xintercept = mean(proportion.p1) - sd(proportion.p1)), color = "blue", linetype = "dotted") +
  geom_vline(aes(xintercept = mean(proportion.p1) + sd(proportion.p1)), color = "blue", linetype = "dotted") +
  labs(title = "Density Plot with Standard Deviation Bands", x = "Proportion P1", y = "Density")
# Standard deviation seems high.

species.prop1.dist <- ggplot(cyt_data_clean_species, aes(x = treatment, y = proportion.p1, fill = treatment)) +
  geom_quasirandom(fill = "lightseagreen", 
                   color = "lightseagreen") +
  theme_classic()+
  ylim(0,6.5) +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("C", "H")), label = "p.signif")+
  facet_wrap(~species)
species.prop1.dist

species.p1 <- ggplot(cyt_data_clean_species, aes(x = treatment, y = events.p1, fill = treatment)) +
  geom_jitter(fill = "#506e40", 
              colour = "#506e40",
              na.rm = TRUE,
              position = position_jitter(height = 0, width = .1),
              alpha = .5)  +
  ylab("Symbiont cell count") +
  ylim(0,80000)+
  theme_classic()+
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("C", "H")), label = "p.signif")+
  facet_wrap(~species)
species.p1


species.p3 <- ggplot(cyt_data_clean_species, aes(x = treatment, y = events.p3, fill = treatment)) +
  geom_jitter(fill = "#f08080", 
              colour = "#f08080",
              na.rm = TRUE,
              position = position_jitter(height = 0, width = .1),
              alpha = .5)  +
  ylab("Coral cell count") +
  ylim(0,150000) +
  theme_classic()+
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("C", "H")), label = "p.signif") +
  facet_wrap(~species) 
species.p3

combined_plot <- plot_grid(species.p1, species.p3, ncol = 1)
combined_plot


# All outliers are from A. abrotanoides colony 76
# Any notes from colony 76 flow cytometry run?

# Add column grouping by species & treatment to make stats easier
cyt_data_clean_species$spec_trt <- paste(as.character(cyt_data_clean_species$species), as.character(cyt_data_clean_species$treatment))
head(cyt_data_clean_species)

# Summary of NAs in the dataset
sapply(cyt_data_clean_species, function(x) sum(is.na(x))) 
cyt_data_na_omit <- na.omit(cyt_data_clean_species)

# Save this dataset to work from:
# write.csv(cyt_data_na_omit, 
          # file = "~/Lab Notebook/Chapter3/Thermal_stress_analysis/cleaned.data.symbiont.counts.csv", 
          # quote = FALSE,
          # row.names = FALSE)

# Clean-up
rm(cyt_data, cyt_data_clean_species)
cyt_data <- cyt_data_na_omit
rm(cyt_data_na_omit)

#####
# Summary Statistics
#####
# now look at summary statistics
# by treatment
# all events
tapply(cyt_data$events.all, cyt_data$treatment, summary)

# P1 events by treatment
tapply(cyt_data$events.p1, cyt_data$treatment, summary)

# P3 events by treatment
tapply(cyt_data$events.p3, cyt_data$treatment, summary)

# proportion by treatment
tapply(cyt_data$proportion.p1, cyt_data$treatment, summary)

# by individual & treatment
# all events
tapply(cyt_data$events.all, cyt_data$colony.id, summary)

# P1 events
tapply(cyt_data$events.p1, cyt_data$colony.id, summary)

# P3 events
tapply(cyt_data$events.p3, cyt_data$colony.id, summary)

# P1 proportion
tapply(cyt_data$proportion.p1, cyt_data$colony.id, summary)

# by species & treatment
# all events
tapply(cyt_data$events.all, cyt_data$spec_trt, summary)

# P1 events
tapply(cyt_data$events.p1, cyt_data$spec_trt, summary)

# P3 events
tapply(cyt_data$events.p3, cyt_data$spec_trt, summary)

# P1 proportion
tapply(cyt_data$proportion.p1, cyt_data$spec_trt, summary)

# Variance in the data
# by treatment
var_trt <- group_by(cyt_data, treatment) %>% summarize(var(events.all),var(events.p1),var(events.p3))
std_trt <- sqrt(var_trt[,1:3])

# by individual
var_ind <- group_by(cyt_data, colony.id, treatment) %>% summarize(var(events.all),var(events.p1),var(events.p3))
std_ind <- sqrt(var_ind[,1:3])
ind_stat <- cbind(var_ind, std_ind) 

# by species
var_spec <- group_by(cyt_data, species, treatment) %>% summarize(var(events.all),var(events.p1),var(events.p3))
std_spec <- sqrt(var_spec[,1:3]) ## Double check this if we ever use these results

#####
# Plotting P1, P3, & P1:P3 data
#####
# looking by treatment
# find overall average P1 events by treatment
mean.p1.events.trt <- tapply(cyt_data$events.p1, cyt_data$treatment, mean)
mean.p1.events.trt <- as.data.frame(mean.p1.events.trt)
mean.p1.events.trt$trtmt <- rownames(mean.p1.events.trt)
colnames(mean.p1.events.trt) <- c('mean_P1_events', 'treatment')

# and plot
trt.plot.p1events <- ggplot(mean.p1.events.trt, aes(x=treatment, y=mean_P1_events, fill=treatment))+
  theme_bw()+
  geom_col()+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+
  ggtitle('P1 Events per Treatment')
ggplotly(trt.plot.p1events)

# find overall average P3 events by treatment
mean.p3.events.trt <- tapply(cyt_data$events.p3, cyt_data$treatment, mean)
mean.p3.events.trt <- as.data.frame(mean.p3.events.trt)
mean.p3.events.trt$trtmt <- rownames(mean.p3.events.trt)
colnames(mean.p3.events.trt) <- c('mean_P3_events', 'treatment')

# and plot
trt.plot.p3events <- ggplot(mean.p3.events.trt, aes(x=treatment, y=mean_P3_events, fill=treatment))+
  theme_bw()+
  geom_col()+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+
  ggtitle('P3 Events by Treatment')
ggplotly(trt.plot.p3events)

# looking by individual
# find overall average P1 events by individual
mean.p1.events.ind <- tapply(cyt_data$events.p1, cyt_data$colony.id, mean)
mean.p1.events.ind <- as.data.frame(mean.p1.events.ind)
mean.p1.events.ind$colony.id <- rownames(mean.p1.events.ind)
mean.p1.events.ind$colony.id <- gsub("[CH]$", "", mean.p1.events.ind$colony.id)

mean.p1.events.ind <-mean.p1.events.ind %>%
  tibble::rownames_to_column(var = "Sample") %>%
  mutate(treatment = ifelse(grepl("C$", Sample), "C", "H")) %>%
  dplyr::select(-Sample)

p1.limits <- mean.p1.events.ind$colony.id
ind.plot.p1events <- ggplot(mean.p1.events.ind, aes(x=colony.id, y=mean.p1.events.ind, fill=treatment))+
  theme_bw()+
  geom_col()+
  scale_x_discrete(limits = p1.limits)+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+
  ggtitle('P1 Events per Colony')
ggplotly(ind.plot.p1events)

# let's look just at paired comparisons between individuals we have both heated and control for
ind.plot.p1events.paired <- ggplot(mean.p1.events.ind, aes(x=colony.id, y=mean.p1.events.ind, fill=treatment))+
  theme_bw()+
  geom_col()+
  scale_x_discrete(limits = p1.limits)+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+
  ggtitle('P1 Events per Colony')
ggplotly(ind.plot.p1events.paired)

# find overall average P3 events by individual
mean.p3.events.ind <- tapply(cyt_data$events.p3, cyt_data$colony.id, mean)
mean.p3.events.ind <- as.data.frame(mean.p3.events.ind)
mean.p3.events.ind$colony.id <- rownames(mean.p3.events.ind)
mean.p3.events.ind$colony.id <- gsub("[CH]$", "", mean.p3.events.ind$colony.id)

mean.p3.events.ind <-mean.p3.events.ind %>%
  tibble::rownames_to_column(var = "Sample") %>%
  mutate(treatment = ifelse(grepl("C$", Sample), "C", "H")) %>%
  dplyr::select(-Sample)

mean.p3.events.ind

# and plot 
p3.limits<- mean.p3.events.ind$colony.id
ind.plot.p3events <- ggplot(mean.p3.events.ind, aes(x=colony.id, y=mean.p3.events.ind, fill=treatment))+
  theme_bw()+
  geom_col()+
  scale_x_discrete(limits = p3.limits)+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+
  ggtitle('P3 Events per Colony')
ggplotly(ind.plot.p3events)



# look at just paired individuals
ind.plot.p3events.paired <- ggplot(mean.p3.events.ind, aes(x=colony.id, y=mean.p3.events.ind, fill=treatment))+
  theme_bw()+
  geom_col()+
  scale_x_discrete(limits = p3.limits)+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+
  ggtitle('P3 Events per Colony')
ggplotly(ind.plot.p3events.paired)

# looking by species (& treatment)
# find overall average P1 events by species & treatment
mean.p1.events.spec <- tapply(cyt_data$events.p1, cyt_data$spec_trt, mean)
mean.p1.events.spec <- as.data.frame(mean.p1.events.spec)
mean.p1.events.spec$species <- rownames(mean.p1.events.spec)
mean.p1.events.spec$species <- gsub("[CH]$", "", mean.p1.events.spec$species)
mean.p1.events.spec <-mean.p1.events.spec %>%
  tibble::rownames_to_column(var = "Sample") %>%
  mutate(treatment = ifelse(grepl("C$", Sample), "C", "H")) %>%
  dplyr::select(-Sample)

# and plot
p1.species <-mean.p1.events.spec$species
p1.species.colors <- c("C" = "#028090", "H" = "#e29578")
spec.plot.p1events <- ggplot(mean.p1.events.spec, aes( x = species, y=mean.p1.events.spec, fill=treatment))+
  theme_bw()+
  geom_col()+
  scale_x_discrete(limits = p1.species)+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+
  scale_fill_manual(values = p1.species.colors)+
  ggtitle('P1 Events by Species')
ggplotly(spec.plot.p1events)

# find overall average P3 percentage by species and treatment
mean.p3.events.spec <- tapply(cyt_data$events.p3, cyt_data$spec_trt, mean)
mean.p3.events.spec <- as.data.frame(mean.p3.events.spec)
mean.p3.events.spec$species <- rownames(mean.p3.events.spec)
mean.p3.events.spec$species <- gsub("[CH]$", "", mean.p3.events.spec$species)
mean.p3.events.spec <-mean.p3.events.spec %>%
  tibble::rownames_to_column(var = "Sample") %>%
  mutate(treatment = ifelse(grepl("C$", Sample), "C", "H")) %>%
  dplyr::select(-Sample)

# and plot
p3.species <-mean.p3.events.spec$species 
p3.species.colors <- c("C" = "#1a5e63", "H" = "#f05a29") 
spec.plot.p3events <- ggplot(mean.p3.events.spec, aes(x=species, y=mean.p3.events.spec, fill=treatment))+
  theme_bw()+
  geom_col()+
  scale_x_discrete(limits = p3.species)+
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))+
  scale_fill_manual(values = p3.species.colors)+
  ggtitle('P3 Events per Species')
ggplotly(spec.plot.p3events)

# Proportion of P1:P3 (symbiont:coral)
p1_proportion_plot<- ggplot(cyt_data, aes(x = sample, y = proportion.p1, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE)) +
  facet_wrap(~species, scales = 'free_x', ncol = 3) +
  ylim(0, 5)

p1_proportion_plot
ggplotly(p1_proportion_plot)

#####
# Normalizing P1 with protein data
#####




#####
# Correlation between species level ED50 and other heat stress metrics
#####
# Read in ED50 data
ED50.data <- read.csv("~/Lab Notebook/Chapter3/Thermal_stress_analysis/Phylogeny resolved ED50 curves/individual.specied.ED50.csv", 
                       row.names = NULL)

# Merge symbiont count and ED50 dataframes
cyt_data$mean.ed50 <- ED50.data[match(cyt_data$colony, ED50.data$colony.id),2]
cyt_data$ind.ed50 <- ED50.data[match(cyt_data$colony, ED50.data$colony.id),3]

# Average replicates across individuals
cyt_data <- cyt_data %>%
  group_by(colony.id) %>%
  mutate(avg.proportion.p1 = mean(proportion.p1))

# Condensed cyt_data
cyt_data_propP1_C_H <- cyt_data[,c(1,2,4,5,6,8,9,27,28,29,30)]
collapsed.cyt.data <- cyt_data_propP1_C_H %>%
  group_by(species, colony.id, colony, treatment, mean.ed50, ind.ed50) %>%
  dplyr::summarize(
    mean.events.p1 = mean(events.p1),
    mean.events.p3 = mean(events.p3)
  ) %>%
  ungroup()
collapsed.cyt.data$proportion.p1 <- collapsed.cyt.data$mean.events.p1/collapsed.cyt.data$mean.events.p3

collapsed.cyt.data <- collapsed.cyt.data %>%
  group_by(colony) %>%
  dplyr::mutate(diff_proportion = proportion.p1[treatment == "C"] - proportion.p1[treatment == "H"])


# Plot the relationship between P1:P3 and ind rel therm tol
ggplot(collapsed.cyt.data, aes(x = ind.ed50, y = proportion.p1)) +
  geom_point() +  # Add points
  geom_smooth(method = "lm", se = FALSE, formula = y~x) +
  labs(x = "Individual Relative Thermal Tolerance",
       y = "Proportion P1",
       title = "Proportion P1 vs Individual Relative Thermal Tolerance") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 36.5, label.y = 3, size = 3)


filtered.data.propP1 <- collapsed.cyt.data %>%
  filter(proportion.p1 < 2)


ggplot(filtered.data.propP1, aes(x = ind.ed50, y = proportion.p1, color = treatment)) +
  geom_point() +  # Add points
  geom_smooth(method = "lm", se = FALSE, formula = y~x) +
  labs(x = "Individual Relative Thermal Tolerance",
       y = "Proportion P1",
       title = "Proportion P1 vs Individual Relative Thermal Tolerance") +
 # facet_wrap(~species) #+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 36.5, label.y = 3, size = 3)
mod.prop.p1<-lm(proportion.p1 ~ ind.ed50 + species + treatment, data = filtered.data.propP1)
summary(mod.prop.p1) # significant relationship between the percent methylation mean and the log10 base means

correlation.ind.ED50 <- cor.test(filtered.data.propP1$proportion.p1, filtered.data.propP1$ind.ed50)
correlation.mean.ED50 <- cor.test(filtered.data.propP1$proportion.p1, filtered.data.propP1$mean.ed50)


# Difference
ggplot(collapsed.cyt.data, aes(x = ind.ed50, y = proportion.p1)) +
  geom_point() +  # Add points
  geom_smooth(method = "lm", se = FALSE, formula = y~x) +
  labs(x = "Individual Relative Thermal Tolerance",
       y = "Proportion P1",
       title = "Proportion P1 vs Individual Relative Thermal Tolerance") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 36.5, label.y = 3, size = 3)
# Remove the outlier
filtered.data <- collapsed.cyt.data %>%
  dplyr::filter(diff_proportion > -2)

library(ggpubr)
library(ggpmisc)
ggplot(filtered.data, aes(x = ind.ed50, y = diff_proportion)) +
  geom_point() +  # Add points
  geom_smooth(method = "lm", formula= y~x) +
  labs(x = "Individual Relative Thermal Tolerance",
       y = "Proportion P1",
       title = "Proportion P1 vs Individual Relative Thermal Tolerance")  +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 36.5, label.y = 0.5, size = 3)


#####
# Testing Cor between ED50 and R channel intensity
#####
# load in channel intensity data
color.data <- read.csv("~/Lab Notebook/Chapter3/Thermal_stress_analysis/Symbiont_counts/Data/Color Data.csv")

setdiff(unique(color.data$Colony_no),
        unique(filtered.data.propP1$colony))

# Filter by colonies included in all analyses (including gene expression)
filtered.color.data <-color.data %>%
                      filter(Colony_no %in% filtered.data.propP1$colony) 

setdiff(unique(filtered.data.propP1$colony),
        unique(filtered.color.data$Colony_no))
# Colonies 90 and 94 are missing from the color data

# Find the average R intensity for each fragment
R.avg.df <- filtered.color.data  %>% 
  group_by(Colony_no, Tank_no) %>% 
  summarize_at(vars(R), list(R.avg = mean))

difference.R.avg <- R.avg.df %>%
     filter(Tank_no %in% c(1, 3)) %>%
     group_by(Colony_no) %>%
     reframe(diff_R_avg = diff(R.avg))
colnames(difference.R.avg) <- c("colony.id", "diff_R_avg")




# Merge the colony (individual) ED50 an R channel Dataframes
ED50.R.channel.diff <- inner_join(ED50.data, difference.R.avg, by = "colony.id")

# write.csv(ED50.R.channel.diff,
#           file = "~/Lab Notebook/Chapter3/Thermal_stress_analysis/ED50.R.channel.diff.csv",
#           quote = FALSE,
#           row.names = FALSE)

ggplot(ED50.R.channel.diff, aes(x = ind.rel.therm.tol, y = diff_R_avg)) +
  geom_point() +  
  geom_smooth(method = "lm", se = FALSE, formula = y~x) +
  labs(x = "Individual Relative Thermal Tolerance",
       y = "Difference in R.avg") +
  ggtitle("Relationship between Individual Relative Thermal Tolerance and Difference in R.avg") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 36.5, label.y = 0.5, size = 3)

# I just want to see the Tank 3 R channel intensity plotted against the individual ED50
tank3.R.intensity <- R.avg.df %>%
  filter(Tank_no == 3)
colnames(tank3.R.intensity) <- c("colony.id","tank.no", "R.avg")

ED50.R.channel.tank3 <- inner_join(ED50.data, tank3.R.intensity, by = "colony.id")
ggplot(ED50.R.channel.tank3, aes(x = ind.rel.therm.tol, y = R.avg)) +
  geom_point() +  
  geom_smooth(method = "lm", se = FALSE, formula = y~x) +
  labs(x = "Individual Relative Thermal Tolerance",
       y = " R.avg") +
  ggtitle("Relationship between Individual Relative Thermal Tolerance and R.avg")

# Lets just look at Red channel intensity:
head(filtered.color.data)
names(filtered.color.data)[names(filtered.color.data) == "Colony_no"] <- "colony.id"
color.data.ed50 <- inner_join(ED50.data, filtered.color.data, by = "colony.id")
color.data.ed50$treatment <- ifelse(color.data.ed50$Tank_no == 1, 29,
                       ifelse(color.data.ed50$Tank_no == 2, 34,
                              ifelse(color.data.ed50$Tank_no == 3, 36,
                                     ifelse(color.data.ed50$Tank_no == 4, 38, NA))))
color.data.ed50$treatment <- as.factor(color.data.ed50$treatment)

# write.csv(color.data.ed50,
# file = "~/Lab Notebook/Chapter3/Thermal_stress_analysis/color.data.ed50.csv",
# quote = FALSE,
# row.names = FALSE)

# library(lme4)
# library(lmerTest)
# # Fit linear mixed-effects model with lmer
# lmer_model <- lmer(R ~ as.numeric(treatment) + (1|species), data = color.data.ed50)
# 
# # Predict values using the lmer model
# color.data.ed50$predicted_R <- predict(lmer_model, newdata = color.data.ed50)
# 
# species.red.channel <- ggplot(color.data.ed50, aes(x = treatment, y = R, color = species)) +
#   geom_jitter(na.rm = TRUE,
#               position = position_jitter(height = 0, width = .1),
#               alpha = .5)  +
#   geom_smooth(method = "lm", se = FALSE, size = 1) +  # Add line representing fitted values
#   stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = "right", label.y = "top", size = 4 ) +
#   ylab("R channel intensity") +
#   xlab("Temperature") +
#   theme_classic() +
#   scale_color_manual(values = rev(official.species.color.palette))
#   facet_wrap(~species)
# 
# # Print summary of the lmer model
# summary(lmer_model)
# 
# species.red.channel

official.species.colors <- c("Acropora tutuilensis" = "#F7C11E",
                    "Acropora lutkeni" = "#FF9A17",
                    "Acropora pulchra" = "#fb5607",
                    "Acropora retusa" = "#ff006e",
                    "Acropora nasuta" = "#c11cad",
                    "Acropora robusta" = "#8338ec",
                    "Acropora hyacinthus" = "#3a86ff",
                    "Acropora abrotanoides" = "#3E6FCB")

species.order <- c("Acropora tutuilensis", "Acropora lutkeni",
                   "Acropora pulchra", "Acropora retusa","Acropora nasuta","Acropora robusta","Acropora hyacinthus","Acropora abrotanoides")
color.data.ed50$species <- as.character(color.data.ed50$species)
color.data.ed50$species[color.data.ed50$species == "Acropora sp"] <- "Acropora nasuta"
color.data.ed50$species <- factor(color.data.ed50$species, levels = species.order)

species.red.channel <- ggplot(color.data.ed50, aes(x = treatment, y = R, color = species)) +
  geom_jitter(na.rm = TRUE,
              position = position_jitter(height = 0, width = .1),
              alpha = .5)  +
  geom_smooth(method = lm, aes(x = as.numeric(treatment)), se = FALSE) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = "right", label.y = "top", size = 4 ) +
  ylab("R channel intensity") +
  xlab("Temperature") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values =official.species.colors ) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 35, label.y = 250, size = 6 ) +
  facet_wrap(~species)

#   
# 
species.red.channel

# Slopes of the lines:

# Acropora abrotanoides
a.abr.model <- lm(R ~ as.numeric(treatment), data=subset(color.data.ed50, species=="Acropora abrotanoides"))
a.abro.slope <- coef(a.abr.model)[2]


# Acropora hycinthus
a.hya.model <- lm(R ~ as.numeric(treatment), data=subset(color.data.ed50, species=="Acropora hyacinthus"))
a.hya.slope <- coef(a.hya.model)[2]

# Acropora lutkeni
a.lut.model <- lm(R ~ as.numeric(treatment), data=subset(color.data.ed50, species=="Acropora lutkeni"))
a.lut.slope <- coef(a.lut.model)[2]

# Acropora pulchra
a.pul.model <- lm(R ~ as.numeric(treatment), data=subset(color.data.ed50, species=="Acropora pulchra"))
a.pul.slope <- coef(a.pul.model)[2]

# Acropora retusa
a.ret.model <- lm(R ~ as.numeric(treatment), data=subset(color.data.ed50, species=="Acropora retusa"))
a.ret.slope <- coef(a.ret.model)[2]

# Acropora robusta
a.rob.model<- lm(R ~ as.numeric(treatment), data=subset(color.data.ed50, species=="Acropora robusta"))
a.rob.slope <- coef(a.rob.model)[2]

# Acropora nasuta
a.nas.model<- lm(R ~ as.numeric(treatment), data=subset(color.data.ed50, species=="Acropora nasuta"))
a.nas.slope <- coef(a.nas.model)[2]

# Acropora tutuilensis
a.tut.model<- lm(R ~ as.numeric(treatment), data=subset(color.data.ed50, species=="Acropora tutuilensis"))
a.tut.slope <- coef(a.tut.model)[2]

# Make a data frame for all the slopes
slopes.df <- data.frame(
  species = c("Acropora abrotanoides", "Acropora hyacinthus", "Acropora lutkeni",
              "Acropora pulchra", "Acropora retusa", "Acropora robusta",
              "Acropora nasuta", "Acropora tutuilensis"),
  slope = c(a.abro.slope, a.hya.slope, a.lut.slope,
            a.pul.slope, a.ret.slope, a.rob.slope,
            a.nas.slope, a.tut.slope)
)


# Merge slopes with ED.50 data
color.data.ed50 <- merge(color.data.ed50, slopes.df , by = "species")
condensed.color.data.ed50 <- unique(color.data.ed50[,c(1,2,16)])
condensed.color.data.ed50$species <- factor(condensed.color.data.ed50$species, levels = species.order)
# # Plot
# spearman_cor <- cor.test(condensed.color.data.ed50$mean.therm.tol, condensed.color.data.ed50$slope, method = "spearman")
# library(ggpubr)
# r.chan.slope.ed50.plot <- ggplot(data = condensed.color.data.ed50, aes(x = mean.therm.tol, y = slope)) +
#   geom_point(aes(color = species), size = 5) +
#   theme_classic() +
#   ylim(0, 35) +
#   labs(x = "Mean Thermal Tolerance", y = "Slope") +
#   ggtitle("Relationship between Slope and Mean Thermal Tolerance") +
#   stat_cor(
#     aes(label = paste("Spearman: r = ", round(spearman_cor$estimate, 2), ", p = ", round(spearman_cor$p.value, 3), sep = "")),
#     label.x = 35, label.y = 32, size = 5
#   ) +
#   scale_color_manual(values = official.species.colors)
# 
# # Display the plot
# r.chan.slope.ed50.plot

spearman_result <- cor.test(condensed.color.data.ed50$mean.therm.tol, condensed.color.data.ed50$slope, method = "spearman")

# Extract the correlation coefficient and p-value
spearman_cor <- round(spearman_result$estimate, 2)
spearman_p <- format.pval(spearman_result$p.value, digits = 2)

r.chan.slope.ed50.plot <- ggplot(condensed.color.data.ed50, aes(x = mean.therm.tol, y = slope)) +
  theme_classic(base_size = 15) +
  ylim(10,29) +
  labs(title = "Correlation between\nRed Channel Slope and ED50",
       x = "Estimated species thermal tolerance (ED50)",
       y = "Red Channel Slope") +
  geom_point(aes(color = species), size = 5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = official.species.colors) +  
  labs(color = "Species") +
  annotate("text", x = Inf, y = Inf, label = paste("Spearman's ρ =", spearman_cor, "\np-value =", spearman_p), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black", fontface = "italic")

# Print the plot
print(r.chan.slope.ed50.plot)

library(patchwork)
species.red.channel + r.chan.slope.ed50.plot

## Remove robusta as an outlier
spearman_result <- cor.test(condensed.color.data.ed50$mean.therm.tol[-6], condensed.color.data.ed50$slope[-6], method = "spearman")

# Extract the correlation coefficient and p-value
spearman_cor <- round(spearman_result$estimate, 2)
spearman_p <- format.pval(spearman_result$p.value, digits = 2)

r.chan.slope.ed50.plot <- ggplot(condensed.color.data.ed50[-6,], aes(x = mean.therm.tol, y = slope)) +
  theme_classic(base_size = 15) +
   labs(title = "Correlation between\nRed Channel Slope and ED50",
       x = "Estimated species thermal tolerance (ED50)",
       y = "Red Channel Slope") +
  geom_point(aes(color = species), size = 5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = official.species.colors[-6]) +  
  labs(color = "Species") +
  annotate("text", x = Inf, y = Inf, label = paste("Spearman's ρ =", spearman_cor, "\np-value =", spearman_p), 
           hjust = 1.1, vjust = 1.1, size = 5, color = "black", fontface = "italic")

# Print the plot
print(r.chan.slope.ed50.plot)


# Test the correlation between slope and ed50:
cor.test(color.data.ed50$slope, color.data.ed50$mean.therm.tol, method = "pearson")
# 
# Let's collate slopes of the individuals
# Define a vector of colonies
colonies <- unique(color.data.ed50$colony.id)

# Create an empty list to store the slopes
slopes.list <- vector("numeric", length(colonies))

# Loop through each individual
for (i in seq_along(colonies)) {
  # Extract the current colony ID
  current.colony <- colonies[i]
  
  # Fit linear model and extract slope
  model <- lm(R ~ as.numeric(treatment), data = subset(color.data.ed50, colony.id == current.colony))
  slopes.list[i] <- coef(model)[2]
}

# Create a data frame for all the slopes
colony.slopes.df <- data.frame(
  colony.id = colonies,
  slope = slopes.list
)

lookup_table <- unique(ED50.data[, c("colony.id", "species","ind.rel.therm.tol")])

colony.slopes.df <- merge(colony.slopes.df, lookup_table, by = "colony.id")


colony.scatter.plot <- ggplot(data = colony.slopes.df, aes(x = ind.rel.therm.tol, y = slope)) +
  geom_point(aes(color = species)) +
  geom_smooth(method = "lm") +  # Add linear regression line
  ylim(-5,35) +
  labs(x = "Individual ED50", y = "Slope") +  # Labels
  ggtitle("Relationship between Individual Slope and ED50") + # Title
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 35, label.y = 32, size = 5 )


# modeling the red channel intensity for each species with individual as a random effect
# Red~Temp + (1|individual)

library(lme4)

# Acropora abrotanoides
a.abr.ind.re.model <- lmer(R ~ as.numeric(treatment) + (1|colony.id), data=subset(color.data.ed50, species=="Acropora abrotanoides"))
coef(a.abr.ind.re.model)
check_model(a.abr.ind.re.model)
a.abro.diagnostic.plots <- plot(check_model(a.abr.ind.re.model))
# Looks like data are left skewed, potentially suggesting there 
# could be a better model for this data. Small amount of data could be
# influencing this though. Take a look at other species. 5 colonies here

# Acropora hyacinthus
a.hya.ind.re.model <- lmer(R ~ as.numeric(treatment) + (1|colony.id), data=subset(color.data.ed50, species=="Acropora hyacinthus"))
coef(a.hya.ind.re.model)
check_model(a.hya.ind.re.model)
a.hya.diagnostic.plots <- plot(check_model(a.hya.ind.re.model)) 
# Model works better for A. hyacinthus. 5 colonies

# Acropora lutkeni
a.lut.ind.re.model <- lmer(R ~ as.numeric(treatment) + (1|colony.id), data=subset(color.data.ed50, species=="Acropora lutkeni"))
coef(a.lut.ind.re.model)
check_model(a.lut.ind.re.model)
a.lut.diagnostic.plots <- plot(check_model(a.lut.ind.re.model)) 
# Model doesn't work as well for A. lutkeni. There seems to be some data points that 
# high higher residuals. Only 3 colonies here

# Acropora pulchra
a.pul.ind.re.model <- lmer(R ~ as.numeric(treatment) + (1|colony.id), data=subset(color.data.ed50, species=="Acropora pulchra"))
coef(a.pul.ind.re.model)
check_model(a.pul.ind.re.model)
a.pul.diagnostic.plots <- plot(check_model(a.pul.ind.re.model)) 
# Model works well for pulchra. 5 pulchra colonies

# Acropora retusa
a.ret.ind.re.model <- lmer(R ~ as.numeric(treatment) + (1|colony.id), data=subset(color.data.ed50, species=="Acropora retusa"))
coef(a.ret.ind.re.model)
check_model(a.ret.ind.re.model)
a.ret.diagnostic.plots <- plot(check_model(a.ret.ind.re.model)) 
# model doesn't work as well for retusa. There is less homogeneity of variance 
# across the data points. 3 colonies

# Acropora robusta
a.rob.ind.re.model <- lmer(R ~ as.numeric(treatment) + (1|colony.id), data=subset(color.data.ed50, species=="Acropora robusta"))
coef(a.rob.ind.re.model)
check_model(a.rob.ind.re.model)
a.rob.diagnostic.plots <- plot(check_model(a.rob.ind.re.model))
# model doesn't work as well here. Residuals don't look randomly distributed.
# 3 colonies

# Acropora sp 
a.sp.ind.re.model <- lmer(R ~ as.numeric(treatment) + (1|colony.id), data=subset(color.data.ed50, species=="Acropora sp"))
coef(a.sp.ind.re.model)
check_model(a.sp.ind.re.model)
a.sp.diagnostic.plots <- plot(check_model(a.sp.ind.re.model))
# 3 colonies. Similar non-randomness in the residuals

# Acropora tutuilensis
a.tut.ind.re.model <- lmer(R ~ as.numeric(treatment) + (1|colony.id), data=subset(color.data.ed50, species=="Acropora tutuilensis"))
coef(a.tut.ind.re.model)
tut.slope <- summary(a.tut.ind.re.model)$coefficients[2]
check_model(a.tut.ind.re.model)
a.tut.diagnostic.plots <- plot(check_model(a.tut.ind.re.model))
# Model seems ok. Bimodality in the data? What might drive this?
# 7 colonies

# Assumption for normality of random effects seem to check out.

# Plot the species slope against species ED50:

# Define a vector of species
species <- unique(color.data.ed50$species)
# Create an empty list to store the slopes
slopes.list <- vector("numeric", length(species))

# Loop through each species
for (i in seq_along(species)) {
  # Extract the current species
  current.species <- species[i]
  # Fit linear model and extract slope
  model <- lmer(R ~ as.numeric(treatment) + (1|colony.id) , data = subset(color.data.ed50, species == current.species))
  slopes.list[i] <- summary(model)$coefficients[2]
}

# Create a data frame for all the slopes
species.slopes.df <- data.frame(
  species = species,
  slope = slopes.list
)

lookup_table <- unique(ED50.data[, c("species","mean.therm.tol")])

species.slopes.df <- merge(species.slopes.df, lookup_table, by = "species")

species.slopes.df$species <- as.character(species.slopes.df$species)

species.slopes.df$species[species.slopes.df$species == "Acropora sp"] <- "Acropora nasuta"

scatter.plot.RE.ind  <- ggplot(data = species.slopes.df  , aes(x = mean.therm.tol, y = slope)) +
  geom_point(aes(color = species), size = 5) +
  theme_classic() + 
  theme(
    plot.title = element_text(size = 16),  # Increase title size
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14)  # Increase legend text size
  ) +
  geom_smooth(method = "lm") +  # Add linear regression line
  ylim(0,35) +
  labs(x = "Mean Thermal Tolerance", y = "Red Channel Slope") +  # Labels
  ggtitle("Relationship between Slope and Mean Thermal Tolerance") + # Title
  scale_color_manual(values = rev(official.species.color.palette)) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 35, label.y = 32, size = 5 )
scatter.plot.RE.ind 

# Remove robusta
species.slopes.no.robusta.df <- species.slopes.df[-which(species.slopes.df$species == "Acropora robusta"), ]
species.slopes.no.robusta.df$species <- factor(species.slopes.no.robusta.df$species, levels = species.order[-6]) 
scatter.plot.RE.ind.no.outlier  <- ggplot(data = species.slopes.no.robusta.df, aes(x = mean.therm.tol, y = slope)) +
  geom_point(aes(color = species), size = 5) +
  geom_smooth(method = "lm") +  # Add linear regression line
  theme_classic() +
  theme(
    plot.title = element_text(size = 16),  # Increase title size
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 16),  # Increase legend title size
    legend.text = element_text(size = 14)  # Increase legend text size
  ) +
  ylim(0,35) +
  labs(x = "Mean Thermal Tolerance", y = " Red Channel Slope") +  # Labels
  ggtitle("Relationship between Slope and Mean Thermal Tolerance") + # Title
  scale_color_manual(values = official.species.colors) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 35, label.y = 32, size = 5 )
scatter.plot.RE.ind.no.outlier 

# Trend is moving in the general direction that we would expect. 

# Adding protein data to normalize control and heated cell counts
prot.data <- read.csv("~/Lab Notebook/Chapter3/Thermal_stress_analysis/Symbiont_counts/Coral Protein Concentration - Sheet1.csv", header = TRUE)
prot.data.avg <- prot.data[,c(1,2,6)]
prot.data.avg <- prot.data.avg[!duplicated(prot.data.avg), ]

# Merge protein data with cell count data.

heat.prot.concentrations <- subset(prot.data.avg, treatment == "heat")
control.prot.concentrations <- subset(prot.data.avg, treatment == "control")

heat.cell.data <- merge(subset(cyt_data_clean_species, treatment == 'H'), heat.prot.concentrations, by = 'colony')
heat.cell.data$norm.p1 <- heat.cell.data$events.p1/heat.cell.data$avg.rep.prot.conc
heat.cell.data$norm.p3 <- heat.cell.data$events.p3/heat.cell.data$avg.rep.prot.conc

control.cell.data <- merge(subset(cyt_data_clean_species, treatment == 'C'), control.prot.concentrations, by = 'colony')
control.cell.data$norm.p1 <- control.cell.data$events.p1/control.cell.data$avg.rep.prot.conc
control.cell.data$norm.p3 <- control.cell.data$events.p3/control.cell.data$avg.rep.prot.conc


norm.cell.data <- rbind(heat.cell.data, control.cell.data)

# Add proportion of symbiont:coral cells (P1:P3)
norm.cell.data$norm.proportion.p1 <- norm.cell.data$norm.p1/norm.cell.data$norm.p3
norm_p1_proportion_plot<- ggplot(norm.cell.data, aes(x = sample, y = norm.proportion.p1, fill = treatment.y)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE)) +
  facet_wrap(~species, scales = 'free_x', ncol = 3) +
  ylim(0, 5)

norm_p1_proportion_plot
ggplotly(norm_p1_proportion_plot)

species.prop1.dist <- ggplot(norm.cell.data, aes(x = treatment.y, y = norm.proportion.p1, fill = treatment.y)) +
  geom_quasirandom(fill = "lightseagreen", 
                   color = "lightseagreen") +
  theme_classic()+
  ylim(0,6.5) +
  theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("control", "heat")), label = "p.signif")+
  facet_wrap(~species)
species.prop1.dist

# Something probs didn't go well with the cell counting. May not include this data 
# in the manuscript.
