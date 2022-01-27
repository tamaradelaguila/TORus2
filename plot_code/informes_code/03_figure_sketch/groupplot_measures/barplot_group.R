# bar plots and anova from the brady + 

setwd("/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures")
rawdata <-read.xls("long_format_forR_group__17filt5.xls") 

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