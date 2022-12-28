# FROM 'regression_R2.R
# REGRESSION LINES (LINEAR MODEL): BRADYCARDIA vs Dm2 ACTIVITY
rm(list = ls())
graphics.off()
# Restart to clear also loaded packages: Ctrl + Shift + (Fn) +  F10 >>> RESTARTS

# SETTINGS
sourcefile ="GROUPbrady_intensity.xls"

sinkfile =  'GROUPbrady_intensity'  #
infotitle= sinkfile

sinkfile = paste(sinkfile)
sinktests = paste('TESTS_' , sinkfile, '.txt')

# GROUP BARPLOT - MEASURES 
library("gdata")

setwd("/home/tamara/Documents/MATLAB/VSDI/TORus/plot/brady2")

rawdata <-read.xls(sourcefile)

# LOAD PACKAGES
library(tidyverse) #for data manipulation and visualization
library(ggpubr) #for creating easily publication ready plots. contains ggqqplot
library(rstatix) #provides pipe-friendly R functions for easy statistical analyses. contains: shapiro_test and levene_test
library("scales")
library (emmeans)
library(car)
library(ggplot2)

# Change 'tibble' options so all columns and rows are printed:
options(tibble.width = Inf)
options(tibble.print_min = Inf)

# LOAD AND RENAME DATA (wmean) - AND TURN INTO FACTORS
data<-rawdata 

names(data)[1]="fish" ;
names(data)[2]="ibi0";
names(data)[3]="mA";

data$fish <- as.factor(data$fish)


library(tidyr)

### A: TODOS LOS PUNTOS (continuo, non-binned)

#% A1.Cálculo del modelo de regresión lineal simple
# 
# modelo_lineal1 <- lm(data$ibi00 ~ data$roiact, data)
# # lm() devuelve el valor de la variable y para x=0 (intersección) junto 
# # con la pendiente de la recta.
# # Para ver la información del modelo se requiere summary().
# summary(modelo_lineal1)
# shapiro.test(modelo_lineal1$residuals) # if p>0 homocedasticity
# # confint(modelo_lineal1) # get intercept
# 
# modelo_lineal2 <- lm(data$roiact~ data$ibi00, data)
# # lm() devuelve el valor de la variable y para x=0 (intersección) junto 
# # con la pendiente de la recta.
# # Para ver la información del modelo se requiere summary().
# summary(modelo_lineal2)
# shapiro.test(modelo_lineal2$residuals) # if p>0 homocedasticity

# A2. GRADO CORRELACIÓN ENTRE VARIABLES

testP<- cor.test(x = data$ibi0, y = data$mA, method = "pearson")
testS<- cor.test(x = data$ibi0, y = data$mA, method = "spearman", exdm2= FALSE)


sink(sinktests)
print('TESTS FOR CORRELATION...........................')
print(testP)
print(testS)

sink()

#% A.3.1 REPRESENTACIÓN GRÁFICA DEL MODELO

tit = paste( 'Diagrama de dispersión + regresión:',  infotitle)
rho <-round(testS$estimate, 2)
p <- round(testS$p.value ,5)
labelcor<-  paste('rhoS=', rho, '(', p, ')')



A2<-ggplot(data = data, aes(x = mA, y = ibi0, colour= fish)) +
  geom_point( size = 2) + #"firebrick"     show.legend = TRUE
  labs(title  =  tit, x  =  "mA") + #'se = TRUE' for plotting intervals # y~model
  geom_smooth(formula = y ~ x, method = "lm", se = FALSE, color = "black") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate(geom= "text", x= 0, y=100,hjust=0,vjust=0 ,label =  labelcor , color='red')
#
svg(paste('REGRESSION_mA_vs_brady_', sinkfile,'.svg'), width = 5, height = 5)
A2
dev.off()
A2


# ### C: 2 REGIONS (non-binned)
# 
# #% C1.Cálculo del modelo de regresión lineal simple
# 
# modelo_lineal1 <- lm(data$ibi0 ~ data$act, data_long)
# # lm() devuelve el valor de la variable y para x=0 (intersección) junto 
# # con la pendiente de la recta.
# # Para ver la información del modelo se requiere summary().
# summary(modelo_lineal1)
# shapiro.test(modelo_lineal1$residuals) # if p>0 homocedasticity
# # confint(modelo_lineal1) # get intercept
# 
# modelo_lineal2 <- lm(data_long$act ~ data_long$ibi0, data_long)
# # lm() devuelve el valor de la variable y para x=0 (intersección) junto 
# # con la pendiente de la recta.
# # Para ver la información del modelo se requiere summary().
# summary(modelo_lineal1)
# shapiro.test(modelo_lineal1$residuals) # if p>0 homocedasticity
# # confint(modelo_lineal1) # get intercept
# 
# 
# # C2. GRADO CORRELACIÓN ENTRE VARIABLES (GET VECTORS FROM data', NOT FROM 'data_long'!)
# 
# testP2_dm2 <- cor.test(x = data$ibi0, y = data$dm2, method = "pearson")
# testS2_dm2 <- cor.test(x = data$ibi0, y = data$dm2, method = "spearman", exdm2= FALSE)
# testS2_dm4 <- cor.test(x = data$ibi0, y = data$dm4, method = "spearman", exdm2= FALSE)
# 
# #% C3.1 REPRESENTACIÓN GRÁFICA DEL MODELO: x=dm2y = BINNED ibi0
# 
# tit = paste( 'ibi0:',  infotitle, '(' ,region ,')')
# rho2 <-round(testS2_dm2$estimate, 2)
# p2 <- round(testS2_dm2$p.value ,3)
# 
# rho4 <-round(testS2_dm4$estimate, 2)
# p4 <- round(testS2_dm4$p.value ,3)
# 
# labelcor3<-  paste('rho 2:', rho2, '(', p2, ')- rho4:', rho4, '(', p4 , ')' )
# 
# C1<- ggplot(data = data_long, mapping = aes(x = ibi0, y = act, color = region, shape= region)) +
#   geom_point(size = 2) + #color = "firebrick", 
#   labs(title  = tit, x  =  'ibi0') + #'se = TRUE' for plotting intervals
#   geom_smooth(method = "lm", se = TRUE, color = "black", fullrange = TRUE, aes(fill=region)) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   annotate(geom= "text", x= -Inf, y=Inf,hjust=-0.1,vjust=2 ,label =  labelcor3 , color='red')
# 
# 
# # SAVE PLOTS 
# 
# jpeg(paste('REGRESSION_2lines_', sinkfile,'.jpg'), width = 500, height = 400)
# C1
# dev.off()
# 
# C1 
# 

# BARPLOTS OF VALUES 
# BARPLOT
library(ggplot2)
library(dplyr)
library(FSA)

# First summarize the data
sum = Summarize(ibi0 ~ mA,
                data=data) 
sum$se = sum$sd / sqrt(sum$n)

# Then plot

pd = position_dodge(0.9)    ### How much to jitter the bars on the plot,
###   0.5 for overlapping bars

printplot <-  barplot(height= sum$mean, main= paste(infotitle),
                      ylab="ibi0", xlab="mA", names.arg= sum$mA)


sum_dm2 = Summarize(dm2 ~ mA,
                    data=data) 
sum_dm2$se = sum_dm2$sd / sqrt(sum_dm2$n)

barplot(height= sum_dm2$mean, main= paste(infotitle),
        ylab="dm2", xlab="mA", names.arg= sum$mA, horiz = TRUE)

sum_dm4 = Summarize(dm4 ~ mA,
                    data=data) 
sum_dm4$se = sum_dm4$sd / sqrt(sum_dm4$n)

barplot(height= sum_dm4$mean, main= paste(infotitle),
        ylab="dm4", xlab="mA", names.arg= sum_dm4$mA,horiz = TRUE)

# barplot_mean <-  ggplot(sum ,                ### The data frame to use.
#                         aes(x= mA,
#                             y= mean)) +
# 
#   geom_bar(stat= "identity",
#            color= "black",
#            position= pd) +
#   
#   geom_errorbar(aes(ymin  = mean - sd,
#                     ymax  = mean + sd),
#                 width = 0.1,
#                 size  = 0.5,
#                 position = pd,
#                 color = "black"
#   ) +
#   theme_grey() +
#   theme(axis.title = element_text(face = "bold")) +
#   theme(legend.position = "bottom")+
#   ylab(paste(measure , "+- sd")) +
#   
#   ggtitle(tit)
# barplot_mean
# 
# barplot_mean

# STATISTICAL COMPARISSON OF LINES: 

lmplot <- ggplot(data,aes(x=ibi0,y=act,col=roi)) + geom_point() + geom_smooth(method="lm")

sink(paste('LMresult_' , sinkfile, '.txt'))
summary(lm(ibi0 ~ act + roi + act:roi, data=data))
sink()

jpeg(paste('LM_', sinkfile,'.jpg'), width = 500, height = 400)
lmplot
dev.off()
lmplot
