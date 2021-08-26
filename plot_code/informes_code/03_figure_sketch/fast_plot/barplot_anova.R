# bar plots and anova from the brady + 

setwd("/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_plot")
#rawdata <-read.csv("long_format_forR.csv") #210521 peak minus baseline
#rawdata <-read.csv("long_format_forR_slopemean210521_400.csv") #210521 peak minus baseline
# rawdata <-read.csv("long_format_forR_wmean210521_400datacirc_filt309.csv") 
#rawdata <-read.csv("long_format_forR_wmean210521_400_datafilt512_allR.csv")
rawdata <-read.csv("plotE_long_format_forR_wmean_between_F11c_p015ms_and_F3c_t10ms_datafilt512.csv") 
rawbradi  <-read.csv("plotE_long_format_forR_bradi_between_F11c_p015ms_and_F3c_t10ms_datafilt512.csv") 

tit= "210521: window 0-250ms mean and brady plots"
file= 'RESULTS 210521 wmean 11p015 vs 3t10ms.txt' # to sink the data in

# Three-way mixed ANOVA: 2 between- and 1 within-subjects factors
# https://www.datanovia.com/en/lessons/mixed-anova-in-r/#three-way-bww-b

# LOAD PACKAGES
library(tidyverse) #for data manipulation and visualization
library(ggpubr) #for creating easily publication ready plots. contains ggqqplot
library(rstatix) #provides pipe-friendly R functions for easy statistical analyses. contains: shapiro_test and levene_test
library("scales")
library (emmeans)
library(car)

# Change 'tibble' options so all columns and rows are printed:
options(tibble.width = Inf)
options(tibble.print_min = Inf)

# LOAD AND RENAME DATA (wmean)
rawdata$cond  <-as.factor(rawdata$Rmean_long1)
rawdata$roi  <- as.factor(rawdata$Rmean_long2)
rawdata$mean  <-rawdata$Rmean_long3

data  <- rawdata[c(4,5,6)]

# LOAD AND RENAME DATA (bradi)

rawbradi$cond  <-as.factor(rawbradi$Rbradi_long1)
rawbradi$bradi  <- rawbradi$Rbradi_long2

dataBradi  <- rawbradi[c(3,4)]

# BARPLOT
library(ggplot2)
library(dplyr)
library(FSA)


# BARPLOT FOR THE MEASURE WMEAN

# First summarize the data 
sum = Summarize(mean ~ roi + cond,
                data=data)
sum$se = sum$sd / sqrt(sum$n)

# Then plot

pd = position_dodge(0.9)    ### How much to jitter the bars on the plot,
###   0.5 for overlapping bars

barplot_mean <-  ggplot(sum,                ### The data frame to use.
       aes(x= cond,
           y= mean,
           fill= roi)) +
  
  geom_bar(stat= "identity",
           color= "black",
           position= pd) +
  
  geom_errorbar(aes(ymin  = mean - se,
                    ymax  = mean + se),
                width = 0.1,
                size  = 0.5,
                position = pd,
                color = "black"
  ) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold")) +
  theme(legend.position = "bottom")+
  ylab("wmean")

barplot_mean

# BARPLOT BRADI
# First summarize the data 

sumBradi = Summarize(bradi ~ cond ,
                     data=dataBradi)
sumBradi$se = sumBradi$sd / sqrt(sumBradi$n)

# Then plot
barplot_bradi <-  ggplot(sumBradi,                ### The data frame to use.
                        aes(x= cond,
                            y= mean)) +
  
  geom_bar(stat= "identity",
           color= "black",
           position= pd) +
  
  geom_errorbar(aes(ymin  = mean - se,
                    ymax  = mean + se),
                width = 0.1,
                size  = 0.5,
                position = pd,
                color = "black"
  ) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold")) +
  ylab("bradi (ibi0)")

barplot_bradi

# BOXPLOT
bxp <- ggboxplot(
  data, x = "cond", y = "mean",
  color = "roi", 
  ggtheme = ggthemes::theme_fivethirtyeight(), # choosing a different theme
) +
  ggtitle('last acq and first reacq')

bxp

# PANEL
library(gridExtra)
grid.arrange(barplot_bradi, barplot_mean, ncol = 1)
barplot_bradi
barplot_mean
#########################################################################
## 2-WAY ANOVA: FORMATED AND SINKING INTO OUTPUT.txt #### 
#########################################################################
library(ggthemes) #https://rafalab.github.io/dsbook/ggplot2.html#add-on-packages
library(wesanderson)

# 1. SUMMARY STATISTICS ####

sink(file)
print(tit)
print("SUMMARY STATISTICS")
print("")
print("Summary grouped by factor ' roi'")
data %>%
  group_by(roi) %>%
  get_summary_stats(value, type = "mean_sd") #value cannot be a factor

print("")
print("Summary grouped by factor 'mA'")
data %>%
  group_by(mA) %>%
  get_summary_stats(value, type = "mean_sd") #value cannot be a factor

print("")
print("Summary grouped by both factors")
data%>%
  group_by(roi,mA) %>%
  get_summary_stats(value, show = c("mean", "sd", "se", "median"))


sink()
# 2. BOXPLOT
# change and save manually:
colormap = wes_palette(name= "Darjeeling1")

bxp <- ggboxplot(
  data, x = "roi", y = "value",
  color = "mA", 
  ggtheme = ggthemes::theme_fivethirtyeight(), # choosing a different theme
  palette = colormap,
) +
  ggtitle(tit)

bxp



#graphics.off()

# 3 . ANOVA: 2-way anova, not repeated measures
sink(file)

print("ANOVA: ANOVA SUMMARY")
anova_summary(aov(value ~roi*mA, data=data))
res.aov <- anova_test(data=data, formula = value ~roi*mA, # + Error(model)
                      #error function: we want to control for that between-participant variation over all of our within-subjects variables.
                      effect.size='pes', type=3)

print("ANOVA: ANOVA TABLE")
#get_anova_table(res.aov)
knitr::kable(get_anova_table(res.aov), caption = 'ANOVA Table (type III tests) - lo mismo, mas bonito')

print("Complete ANOVA object to check Mauchly and Sphericity corrections") --???
  res.aov  

# 4. POST-HOC ANALYSIS

print("POST-HOC with Bonferroni adjustment")

# 4.1: PAIRWISE COMPARISONS each level of the within-s factors

# 4.1.1. Each level of block
# Paired t-test is used because we have repeated measures by time


print("POST-HOC: PAIRWISE COMPARISONS BETWEEN ROI LEVELS  (non paired)")
pwc1 <- data %>%
  group_by(roi) %>%
  pairwise_t_test(
    formula = value ~ mA, paired = FALSE, 
    p.adjust.method = "bonferroni"
  )
#  %>% select(-df, -statistic, -p) # Remove details
#pwc1

knitr::kable(pwc1, caption = 'POST-HOC: PAIRWISE COMPARISONS BETWEEN ROI LEVELS')


print("POST-HOC: PAIRWISE COMPARISONS BETWEEN mA LEVELS (non paired)")
# non-pair to compare EC? because is not a rm?
pwc2 <- data %>%
  group_by(mA) %>%
  pairwise_t_test(
    formula = value ~ roi, paired = FALSE, 
    p.adjust.method = "bonferroni"
  )
#%>%   select(-df, -statistic, -p) # Remove details
#pwc2

knitr::kable(pwc2, caption = 'POST-HOC: PAIRWISE COMPARISONS BETWEEN mA LEVELS')

sink ()



# 5. TEST ASSUMPTIONS  ####

# 5.1: PRESENCE OF OUTLIERS 

print("ASSUMPTIONS: PRESENCE OF OUTLIERS")
data %>%
  group_by(roi, mA) %>%
  identify_outliers(value)

data %>%
  group_by(mA) %>%
  identify_outliers(value)


# 5.2: NORMALITY ASSUMPTION
# The normality assumption can be checked by computing Shapiro-Wilk test for each time point.
# p > 0.05: normal distribution (H0)

print("ASSUMPTIONS: NORMALITY - SHAPIRO-WILK TEST")
data %>%
  group_by(roi, mA) %>%
  shapiro_test(value)

# CREATE QQ plot for each cell of design: (for bigger sample sizes)
#ggqqplot(dm2, "brady", ggtheme = theme_bw()) +
#facet_grid(block ~ EC, labeller = "label_both")
#If all the points fall approx. along the ref.line, for each cell, we can assume normality of the data.


# 5.3. HOMOCEDASTICITY: not needed (not between subjects comparisons)

# 5.4: SPERICITY ASSUMPTION: done and corrected internally by the functions anova_test() and get_anova_table()
# Sphericity is an important assumption of a repeated-measures ANOVA 

print("ASSUMPTIONS: SPERICITY ASSUMPTION is performed and applied internally")

# 5.5 : HOMOGENEITY OF  COVARIANCES  !!! sale diferente
print("ASSUMPTIONS: EQUALITY OF COVARIANCEs - Box's Text")
box_m(data[, "value", drop = FALSE], data$roi)
# The assumption of homogeneity of variance is an assumption of the independent samples t-test and ANOVA stating
#that all comparison groups have the same variance.  The independent samples t-test and ANOVA utilize the t and F statistics respectively,
# which are generally robust to violations of the assumption as long as group sizes are equal
sink()

#######################################3

# BOXPLOT
bxp <- ggboxplot(
  data, x = "mA", y = "value",
  color = "roi", palette = c("jco"),
  #facet.by =  "mA"
  title =tit 
)
bxp

bxp <- ggboxplot(
  data, x = "roi", y = "value",
  color = "mA", palette = c("jco"),
  title =tit 
  
)
bxp

# CHECK ASSUMPTIONS

# A1: OUTLIERS
data %>%
  group_by(roi, mA) %>%
  identify_outliers(value)


# A2: NORMALITY ASSUMPTION
data %>%
  group_by(roi, mA) %>%
  shapiro_test(value)
# p < 0.05: non-normal distribution

# CREATE QQ plot for each cell of design:
ggqqplot(data, "value", ggtheme = theme_bw()) +
  facet_grid(mA ~ roi, labeller = "label_both")

# Create QQ plot for each cell of design:
ggqqplot(data, "value", ggtheme = theme_bw()) +
  facet_grid(roi~mA , labeller = "label_both")

#facet_grid(EC+group ~ sesion, labeller = "label_both")

#If all the points fall approx. along the ref.line, for each cell, we can assume normality of the data.

# A3: HOMOGENEITY OF VARIANCE ASSUMPTION (Levene's test at each level of the within-subjects factor)
data %>%
  group_by(roi) %>%
  levene_test(value ~ mA)

data %>%
  group_by(mA) %>%
  levene_test(value ~ roi)



# Note that, if you do not have homogeneity of variances, you can try to transform the outcome (dependent) variable to correct for the unequal variances.
# If group sample sizes are (approximately) equal, run the three-way mixed ANOVA anyway because it is somewhat robust to heterogeneity of variance in these circumstances.
# It's also possible to perform robust ANOVA test using the WRS2 R package.'''


# A4: SPERICITY ASSUMPTION: done and corrected internally by the functions anova_test() and get_anova_table()

# COMPUTATION OF ANOVA
#anova_summary(aov(value ~ roi*mA + Error(idx/(mA)), data=data))
anova_summary(aov(value ~ roi*mA, data=data))

# TWO-WAY ANOVA
res.aov <- anova_test(data = data, dv = value,
                      between = roi, within = mA)

get_anova_table(res.aov)

# extracted from inside: wid = idx,

# COMPUTATION ONE WAY FOR mA
model <- lm(value ~ mA*roi, data = data)
# effect on each level of mA
data %>%
  group_by(roi) %>%
  anova_test(value ~ mA, error = model)
# effect on each level of roi
data %>%
  group_by(mA) %>%
  anova_test(value ~ roi, error = model)

# TWO-WAY ANOVA
anova_2way <- aov(value~roi*mA, data = data)
summary(anova_2way)

get_anova_table(anova_2way)
## Contrasts set to contr.sum for the following variables: Btw_Cond
knitr::kable(nice(anova_2way)) ##call for formatted ANOVA table using knitr

# summary
require("dplyr")
group_by(data, roi, mA) %>%
  summarise(
    count = n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE)
  )
# POST-HOC
# Using Turkey
TukeyHSD(anova_2way, which = c("roi" , "mA"))


##Main effects
mixed_fitted_roi<-emmeans(anova_2way, ~roi)
mixed_fitted_roi

mixed_fitted_mA<-emmeans(anova_2way, ~mA)
mixed_fitted_mA


mixed_fitted_interactions1 <-emmeans(anova_2way, ~roi|mA)
mixed_fitted_interactions1
pairs(mixed_fitted_interactions1)


##pairwise comparison with no correction.
pairs(mixed_fitted_interactions1)

# LEVENE TEST
test_levene(anova_2way) #lawstat package


# 6. PLOT

# OLD VISUALIZATION : VISUALIZATION: BOXPLOT WITH p-values (not working - because it's two way???)

# visualizacion grouped by mA
pwc1 <- pwc1%>% add_xy_position(x = "mA")

bxp <- ggboxplot(
  data, x = "mA", y = "value",
  color = "roi", 
  ggtheme = ggthemes::theme_fivethirtyeight(), # choosing a different theme
  package = "wesanderson", # package from which color palette is to be taken
  palette = "Darjeeling1" # choosing a different color palette
  # palette = "jco",
) +
  ggtitle('anova result')+
  
  stat_pvalue_manual(pwc1, tip.length = 0, hide.ns = TRUE)+
  labs(
    subtitle = get_test_label(pwc1, detailed = TRUE),
    caption = get_pwc_label(pwc1)
  )
bxp




# link: https://rpkgs.datanovia.com/rstatix/reference/anova_test.html#details

#########################################################################
## ROBUST 2-WAY ANOVA: WRS2
#########################################################################
library(WRS2)
# Median-based (more robust to outliers)
set.seed(123)
#check unused factors levels
b1= nlevels(data$mA)
b2= length(unique(data$mA)) # if b1 and b2 are different, drop unused levels:

# we drop unused levels (otherwise med2way will raise an error)
data$mA <- droplevels(data$mA)
data$roi <- droplevels(data$roi)

med2way(value~roi*mA, data = data)


postHoc1 <- mcp2a(value~roi*mA, data = data,)
postHoc1

# see levels to be able to understand the posthoc
levels(data$roi)
levels(data$roi)

# DUNN's TEST (PAIWISE COMPARISONS) for each level of mA
data  %>%
  group_by(mA) %>%
  dunn_test(value ~ roi)

data  %>%
  group_by(roi) %>%
  dunn_test(value ~ mA)
