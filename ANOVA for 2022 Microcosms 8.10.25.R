# Title: ANOVA for 2022 Microcosms 8.10.25.R
# Author: Nettie Calvin
# Date: August 30, 2025
# Purpose: The script provides code to conduct an analysis of variance (ANOVA) on potential mercury methylation rate constants (kmeth) and potential methylmercury demethylation rate constants (kdemeth) from microcosm experiments conducted in October 2022.  
# The script contains the code to calculate the kdemeth values while the kmeth values are already calculated  prior to the data being added to kmeth Microcosms 2022 8.10.25_for_code.csv.  
#  The script and statistical analysis were prepared to accompany the publication 
# "Nitrate Treatment Suppresses Mercury Demethylation in Coastal Estuarine Sediment" (Calvin et al. in prep)
# R version 4.5.0 (2025-04-11)


setwd("")
#make specific to your working directory

#make sure to save copies of 
#kmeth Microcosms 2022 8.10.25_for_code.csv
#kdemeth Microcosms 2022 8.9.25_for_code.csv
#in your working directory

install.packages(c("ggplot2", "ggpubr", "car"))
#load libraries 
library(ggplot2)
library(ggpubr)
library(car)

#theme for formatting plots
Nettie_Theme = theme( legend.position = "none",
                      title = NULL,
                     axis.title.x = element_text(size = 14, face = "bold"),
                     axis.text.x = element_text(size = 14, face = "bold"),
                     axis.text.y = element_text(size = 14,  face = "bold"),
                     axis.title.y = element_text(size = 14, face = "bold"), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_line(colour = "black"))



#read in kmeth and kdemeth data frames
redo.rates.kmeth <- read.csv("kmeth Microcosms 2022 8.10.25_for_code.csv")
redo.rates.kdemeth = read.csv("kdemeth Microcosms 2022 8.9.25_for_code.csv")

#make the treatment a factor variable for kmeth
redo.rates.kmeth$Treatment = as.factor(redo.rates.kmeth$Treatment)

#then specify the order of factor levels 
redo.rates.kmeth$Treatment <- factor(redo.rates.kmeth$Treatment, levels = c("Control","Low Nitrate", "High Nitrate"))    
redo.rates.kmeth$Treatment

#make the treatment a factor variable for kdemeth
redo.rates.kdemeth$Treatment = as.factor(redo.rates.kdemeth$Treatment)

#then specify the order of factor levels
redo.rates.kdemeth$Treatment <- factor(redo.rates.kdemeth$Treatment, levels = c("Control", "Time0", "Low Nitrate", "High Nitrate"))    
redo.rates.kdemeth$Treatment

#inspect the kmeth df
head(redo.rates.kmeth)

#create a variable for 48 hr control group kmeth 
#(already in the form of potential methylation rate constants with units of percent per day)
control.group = redo.rates.kmeth$percent.per.day[redo.rates.kmeth$Treatment=="Control"]
control.group

#calculate summary stats for 48 hr control group kmeth
median(control.group) #4.015
mean(control.group) #4.079
sd(control.group) #1.671491

#create variables for Low and High Nitrate Treatment groups kmeth
micromolar.group = redo.rates.kmeth$percent.per.day[redo.rates.kmeth$Treatment=="Low Nitrate"]
millimolar.group = redo.rates.kmeth$percent.per.day[redo.rates.kmeth$Treatment=="High Nitrate"]

#calculate summary stats for for Low Nitrate Treatment group kmeth
median(micromolar.group) #3.667
mean(micromolar.group) #3.674667
sd(micromolar.group) #0.5855376

#calculate summary stats for High Nitrate Treatment group kmeth
median(millimolar.group) #3.129
mean(millimolar.group) #2.9708
sd(millimolar.group) # 0.3453559


#potential MeHg demethylation rate constants need to be calculated based on the 
#ng/g concentrations of 198MMHg measured in the t=0 (insta) and t=48 (control) or t=48 High or Low Nitrate (treatment) microcosms
# now inspect the demethylation data frame
head(redo.rates.kdemeth)

#calculate the potential MeHg demethylation rate constants by steps
instas = redo.rates.kdemeth$X198MeHg.ng.g[redo.rates.kdemeth$Treatment =="Time0"]

instas

#calculate an average time 0 to use to divide all other values by
mean.t0=mean(instas)
mean.t0
#1.06016 (ng/g) average t=0 concentration of 198MMHg from spike, not yet a rate constant!!!!
sd(instas) #0.1714231; std. deviation of the instantaneous control samples

#divide the column of ng.g conc. by the avg. t0 and assign a new vector to the df
redo.rates.kdemeth$ratios = redo.rates.kdemeth$X198MeHg.ng.g/mean.t0
redo.rates.kdemeth$ratios

#now take the negative natural log of the ratio
redo.rates.kdemeth$nat.log = -log(redo.rates.kdemeth$ratios)
redo.rates.kdemeth$nat.log 

#divide by 2 days (incubation time) and multiply by 100 to get rate constants in units of percent per day
redo.rates.kdemeth$percent.per.day =  (redo.rates.kdemeth$nat.log/2)*100

#modify the dataframe so it is just the part that is the control and treatment, no time0's
redo.rates.kdemeth.smaller =redo.rates.kdemeth[6:18,]

# drop unused levels from dataframe (Time = 0 condition)
redo.rates.kdemeth.smaller = droplevels(redo.rates.kdemeth.smaller)

#get just the rate constants (kdemeth) from the (48 hour) control group 
control.group.kdemeth = redo.rates.kdemeth.smaller$percent.per.day[redo.rates.kdemeth.smaller$Treatment=="Control"]
control.group.kdemeth

control.median.kdemeth = median(control.group.kdemeth)
control.median.kdemeth #35.68948
control.mean.kdemeth =mean(control.group.kdemeth) 
control.mean.kdemeth #31.70064
sd(control.group.kdemeth) #11.28917

#get the rate constants (kdemeth) from millimolar (High Nitrate) treatment 
millimolar.group.kdemeth = redo.rates.kdemeth.smaller$percent.per.day[redo.rates.kdemeth.smaller$Treatment=="High Nitrate"]
millimolar.group.kdemeth
median(millimolar.group.kdemeth) #12.76616
mean(millimolar.group.kdemeth) #13.16634
sd(millimolar.group.kdemeth) #5.14298

#get the rate constants (kdemeth) from micromolar (Low Nitrate) treatment 
micromolar.group.kdemeth = redo.rates.kdemeth.smaller$percent.per.day[redo.rates.kdemeth.smaller$Treatment=="Low Nitrate"]
micromolar.group.kdemeth
median(micromolar.group.kdemeth) #15.2384
mean(micromolar.group.kdemeth) #14.66646
sd(micromolar.group.kdemeth) #2.054594

#produce the boxplot for kmeth

bxplot_NO3_kmeth =
ggplot(redo.rates.kmeth, aes(x =Treatment, y=percent.per.day, fill = Treatment))+
  geom_hline(yintercept=4.015,linetype=3, col = "black", lwd = 2)+
  stat_boxplot(geom ='errorbar', width = 0.2)+ #this supposedly plots error bars first and then we will plot the boxplot atop them??
  geom_boxplot()+   
  scale_fill_manual(values=c("#FF0000", "#00A08A", "#5BBCD6"))+
  
  #add color here and also overplot see-> # shorthand for  stat_boxplot(geom='boxplot')
  #the width specification inside stat_boxpot is to control the errorbar endcap width
  # from : https://stackoverflow.com/questions/12993545/put-whisker-ends-on-boxplot
  labs( x ="Treatment Group", y="kmeth (percent per day)")+
  Nettie_Theme

#export the plot to the window
bxplot_NO3_kmeth

#create the kdemeth boxplot
#
bxplot_NO3_kdemeth =
ggplot(redo.rates.kdemeth.smaller, aes(x =Treatment, y=percent.per.day, fill = Treatment))+
  geom_hline(yintercept=35.69,linetype=3, col = "black", lwd = 2)+
  stat_boxplot(geom ='errorbar', width = 0.2)+ #this supposedly plots error bars first and then we will plot the boxplot atop them??
  geom_boxplot()+   
  scale_fill_manual(values=c("#FF0000", "#00A08A", "#5BBCD6"))+
  #add color here and also overplot see-> # shorthand for  stat_boxplot(geom='boxplot')
  #the width specification inside stat_boxpot is to control the errorbar endcap width
  # from : https://stackoverflow.com/questions/12993545/put-whisker-ends-on-boxplot
  labs( x ="Treatment Group", y="kdemeth (percent per day)")+
  Nettie_Theme

#export to window
bxplot_NO3_kdemeth

#group both plots together for producing a figure
figure.2022 = ggarrange(bxplot_NO3_kmeth, bxplot_NO3_kdemeth, 
          labels = c("A", "B"), 
          common.legend = TRUE, 
          legend = "none",
          ncol = 2, nrow = 1)

#produce and export a 2 panel figure
png("8.28.25 Boxplot kmeth and kdemeth BZ 2022.png", units='in',width=14,height=6,res=300)
figure.2022
dev.off()


#do a Shapiro-Wilk test to check for normal distribution, kmeth
#
#do a Shapiro-Wilk test to check for normal distribution, kdemeth
#

#run an unbalanced design, one way ANOVA on the kdemeth data
#table of data, is it balanced?  (we know it is not)
table(redo.rates.kdemeth.smaller$Treatment)

#Control  Low Nitrate High Nitrate 
#5            3            5 

#run an unbalanced design, one way ANOVA on the kdemeth data
anova_kdemeth <- aov(percent.per.day ~ Treatment , data = redo.rates.kdemeth.smaller)
Anova(anova_kdemeth, type = "III")
#                Sum Sq Df F value   Pr(>F)    
#(Intercept)    5024.7  1 80.5200 4.25e-06 ***
#  Treatment    998.0  2  7.9966 0.008428 ** 
 # Residuals    624.0 10   
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#produce the Q-Q plots for the kdemeth values
par(mfrow=c(2,2))
plot(anova_kdemeth)
par(mfrow=c(1,1))

#run post hoc tests for the kdemeth ANOVA, which had a significant p value
#the post hoc test will tell us which group(s) are different from the control
tukey.two.way<-TukeyHSD(anova_kdemeth)
tukey.two.way

#diff       lwr       upr     p adj
#Low Nitrate-Control      -17.034185 -32.84872 -1.219654 0.0353486
#High Nitrate-Control     -18.534300 -32.23009 -4.838515 0.0102741
#High Nitrate-Low Nitrate  -1.500115 -17.31465 14.314415 0.9635470


#what do we know?
#treatment is significant p = 0.008428 ** from ANOVA
#0.0353486 Low Nitrate-Control 
#0.0102741 High Nitrate-Control 
#High Nitrate-Low Nitrate  0.9635470
#both NO3 treatments are different than control, 
#nitrate treatments are not diff. than each other



#let's do a one way ANOVA on potential methylation rate constant (kmeth) data 
##table of data, is it balanced?  (we know it is not)
table(redo.rates.kmeth$Treatment)
#Control  Low Nitrate High Nitrate 
#5            3            5 

#run an unbalanced design, one way ANOVA on the kmeth data
anova_kmeth <- aov(percent.per.day ~ Treatment , data = redo.rates.kmeth)
Anova(anova_kmeth, type = "III")

#               Sum Sq Df F value    Pr(>F)    
#   (Intercept) 83.191  1 67.4251 9.365e-06 ***
 # Treatment    3.122  2  1.2652    0.3237  #not a significant difference for Treatment
#  Residuals   12.338 10  

#no sig. differences from ANOVA, no need for post hoc test

#do testing following the ANOVAs to check the residuals' 
#normality of distribution and homogeneity of variance
## let's grab the residuals from our model for kmeth
kmeth.anova.residuals <- residuals( object = anova_kmeth ) # extract the residuals
# A simple histogram
hist( x = kmeth.anova.residuals) 
png("8.30.25.edit kmeth residuals histogram 2022.png", units='in',width=37,height=24,res=300)
hist( x = kmeth.anova.residuals) # another way of seeing them
dev.off()

#Let’s run a test of the normality assumption for kmeth ANOVA model
#the Kolmogorov-Smirnov test (against a normal distribution).
ks.test(kmeth.anova.residuals, "pnorm", mean(kmeth.anova.residuals), sd(kmeth.anova.residuals) )
#	Exact one-sample Kolmogorov-Smirnov test
#data:  kmeth.anova.residuals
#D = 0.22787, p-value = 0.4442
#alternative hypothesis: two-sided
##ok normal enough
##
#The final test is homogeneity of variance also known as (homoscedasticity).
leveneTest(anova_kmeth, center = median) 
# technically a Brown Forsyth test with median as center (a version of Levene's test)
#Levene's Test for Homogeneity of Variance (center = median)
#      Df F value Pr(>F)
#group  2  1.4363 0.2829
#      10  
#interpretation:
#p > 0.05 → Fail to reject the null hypothesis → Variances are equal (homogeneity holds).
# can assume equal variances, ok that ANOVA was used

############
## let's grab the residuals from our model for kdemeth
kdemeth.anova.residuals <- residuals( object = anova_kdemeth ) # extract the residuals
# A simple histogram for kdemeth residuals
hist( x = kdemeth.anova.residuals) 
png("8.30.25.edit kdemeth residuals histogram 2022.png", units='in',width=37,height=24,res=300)
hist( x = kdemeth.anova.residuals) # another way of seeing them
dev.off()
#Let’s run a test of the normality assumption for kdemeth ANOVA model
#the Kolmogorov-Smirnov test (against a normal distribution).
ks.test(kdemeth.anova.residuals, "pnorm", mean(kdemeth.anova.residuals), sd(kdemeth.anova.residuals) )
#	Exact one-sample Kolmogorov-Smirnov test
#data:  kdemeth.anova.residuals
#D = 0.17018, p-value = 0.7878
#alternative hypothesis: two-sided
##ok normal enough
##
#The final is homogeneity of variance also known as (homoscedasticity).
leveneTest(anova_kdemeth, center = median) 
# technically a Brown Forsyth test with median as center (a version of Levene's test)
#Levene's Test for Homogeneity of Variance (center = median)
#      Df F value Pr(>F)
#group  2  0.6798 0.5287
#      10     
#interpretation:
#p > 0.05 → Fail to reject the null hypothesis → Variances are equal (homogeneity holds).
# can assume equal variances, ok that ANOVA was used

