# Title: ANOVA for 2020 Microcosms 8.26.25.R
# Author: Nettie Calvin
# Date: August 29, 2025
# Purpose: The script conducts an analysis of variance (ANOVA) on potential mercury methylation rate constants from microcosm experiments conducted in October 2020.  
# The script also contains code to produce a boxplot figure to visualize the kmeth values.  The script and statistical analysis are prepared to accompany the publication 
# "Nitrate Treatment Suppresses Mercury Demethylation in Coastal Estuarine Sediment" (Calvin et al. in prep)
# R version 4.5.0 (2025-04-11)

setwd("")  #make specific to your working directory

#make sure to save copies of Beach Zone kmeth Data 2020 8.25.25.csv 
#and East Fork kmeth Data 2020 8.25.25.csv
#in your working directory

#install required packages
install.packages(c( "wesanderson", "rstatix", "car", "ggplot2", "ggpubr", "DescTools"))

#load libraries
library(DescTools)
library(ggpubr)
library(ggplot2)
library(car)
library(rstatix)
library(wesanderson)

#read in data frames of kmeth values
bz.kmeth.df=read.csv("Beach Zone kmeth Data 2020 8.25.25.csv")
#Beach Zone kmeth Data 2020 8.25.25.csv
ef.kmeth.df=read.csv("East Fork kmeth Data 2020 8.25.25.csv")
#East Fork kmeth Data 2020 8.25.25.csv


#c("48 hour Control", "Nitrate", "Ammonium", "Molybdate", "Molybdate + Nitrate")
#quickly reassign treatment group labels to the group column
bz.kmeth.df$group=c(rep("48 hour Control", 3), rep("Nitrate", 3), rep("Ammonium", 3), rep("Molybdate", 2), rep("Molybdate + Nitrate", 3))

bz.kmeth.df$group
#woohoo it worked
#
#do the same for EF
#quickly reassign treatment group labels to the group column
ef.kmeth.df$group=c(rep("48 hour Control", 3), rep("Nitrate", 3), rep("Ammonium", 3), rep("Molybdate", 3), rep("Molybdate + Nitrate", 3))

ef.kmeth.df$group

##test for normal distribution
###use wilk-shapiro test for normality
###http://www.sthda.com/english/wiki/normality-test-in-r


#Density plot: the density plot provides a visual judgment about whether the distribution is bell shaped.
ggdensity(bz.kmeth.df$kmeth.percent.per.day, 
          main = "Density plot of BZ potential methylation rate constants",
          xlab = "Percent of Spike Methylated per Day")

#distribution of samples at beach zone is not bell shaped

#Density plot: the density plot provides a visual judgment about whether the distribution is bell shaped.


ggdensity(ef.kmeth.df$kmeth.percent.per.day, 
          main = "Density plot of EF potential methylation rate constants",
          xlab = "Percent of Spike Methylated per Day")

#dev.off()

# looks lognormal at BZ 

##make the q-q and 1:1 line plot
ggqqplot(bz.kmeth.df$kmeth.percent.per.day)
ggqqplot(ef.kmeth.df$kmeth.percent.per.day)

#look pretty good, better at EF

#Note that, normality test is sensitive to sample size. Small samples most often pass normality tests.
##do Wilk-Shapiro
shapiro.test(bz.kmeth.df$kmeth.percent.per.day)
#W = 0.83208, p-value = 0.0128 not normally distributed

shapiro.test(ef.kmeth.df$kmeth.percent.per.day)
#W = 0.93848, p-value = 0.3636
#ef values are normally distibuted enough to pass the shapiro-wilk test

## try a natural log transformation on the BZ data
ln.bz.kmeth =log(bz.kmeth.df$kmeth.percent.per.day) #calculate natural log
ln.bz.kmeth

#assign it back into a column of the bz dataframe
bz.kmeth.df$nat.log.bz.kmeth = ln.bz.kmeth

#now retry the density plots and norm tests

#Density plot: the density plot provides a visual judgment about whether the distribution is bell shaped.
ggdensity(ln.bz.kmeth, 
          main = "Density plot of ln-transformed BZ potential methylation rate constants",
          xlab = "Percent of Spike Methylated per Day")

#now looks bimodal at bz

ggqqplot(ln.bz.kmeth)

#looks pretty good do shapiro test

#Note that, normality test is sensitive to sample size. Small samples most often pass normality tests.
##do Shapiro-Wilk
shapiro.test(ln.bz.kmeth)
#W = 0.90053, p-value = 0.1148

#OK SO NOW PASSES Shapiro-Wilk TEST 
#
#SO ATTEMPT AN ANOVA? ON LN TRANS DATA for BZ

#make the group names factor variables and 
#specify their order as listed in original dataframe (so they won't be alphebetized)
bz.group.factor = factor(bz.kmeth.df$group,
                    levels=unique(bz.kmeth.df$group))
ef.group.factor = factor(ef.kmeth.df$group,
                    levels=unique(ef.kmeth.df$group))


#https://www.scribbr.com/statistics/anova-in-r/
#A ONE-WAY ANOVA (means you have 1 IV)
#run an unbalanced design, one way ANOVA on the bz kmeth data 
#(as molybdate has only 2 groups after removing outlier with dixon's q test)
anova.kmeth.bz <- aov(nat.log.bz.kmeth ~ bz.group.factor , data = bz.kmeth.df)
Anova(anova.kmeth.bz, type = "III")


#nat.log.bz.kmeth
                  #Sum Sq Df F value  Pr(>F)  
#(Intercept)      0.2700  1  0.3637 0.56136  
#bz.group.factor 11.2793  4  3.7987 0.04468 *
#  Residuals     6.6808  9                  
 
#barely significant do post hoc test
#
DunnettTest(x=bz.kmeth.df$nat.log.bz.kmeth, g=bz.group.factor)
#                                     pval    
#Nitrate-48 hour Control             0.3845    
#Ammonium-48 hour Control            0.9304      
#Molybdate-48 hour Control           0.3194    
#Molybdate + Nitrate-48 hour Control 0.4067   

#not sig. try tukeys HSD test

bz.tukey.two.way<-TukeyHSD(anova.kmeth.bz)
bz.tukey.two.way

                                    #p adj
#Nitrate-48 hour Control             0.5429841
#Ammonium-48 hour Control            0.9706583
#Molybdate-48 hour Control           0.4687297
#Molybdate + Nitrate-48 hour Control 0.5672492
#Ammonium-Nitrate                    0.8602707
#Molybdate-Nitrate                   0.0700157
#Molybdate + Nitrate-Nitrate         0.0718531
#Molybdate-Ammonium                  0.2410764
#Molybdate + Nitrate-Ammonium        0.2825215
#Molybdate + Nitrate-Molybdate       0.9963125

#no significant contrasts at BZ

#no log transform at EF run anova on regular data
anova.ef <- aov(ef.kmeth.df$kmeth.percent.per.day ~ ef.group.factor)
summary(anova.ef)
#Based on the ANOVA model the p-value is statistically significant (p<0.05), indicate that each group does not have the same average values.

#                   Df  Sum Sq  Mean Sq F value Pr(>F)  
#ef.group.factor    4 0.07390 0.018475    3.68 0.0431 *
#  Residuals       10 0.05021 0.005021                 
    
#the ANOVA result suggests significance (p<0.05)  but the post hoc test does not indicate any 
#groups are different from control
#dunnett's post hoc test at EF to reveal which groups have a significant difference
DunnettTest(x=ef.kmeth.df$kmeth.percent.per.day, g=ef.group.factor)

#The Pr(>F) column is the p-value of the F-statistic. 
##This shows how likely it is that the F-value calculated 
##from the test would have occurred if 
###the null hypothesis of no difference among group means were true.

                                    ##pval    
#Nitrate-48 hour Control             0.9689    
#Ammonium-48 hour Control             0.9995    
#Molybdate-48 hour Control           0.0782 .  
#Molybdate + Nitrate-48 hour Control 0.0986 .

#try tukey's hsd for comparison
ef.tukey.two.way<-TukeyHSD(anova.ef)
ef.tukey.two.way

#                                    p adj
#Nitrate-48 hour Control             0.9885188
#Ammonium-48 hour Control            0.9998456
#Molybdate-48 hour Control           0.1403638
#Molybdate + Nitrate-48 hour Control 0.1734266
#Ammonium-Nitrate                    0.9676199
#Molybdate-Nitrate                   0.2737281
#Molybdate + Nitrate-Nitrate         0.3307302
#Molybdate-Ammonium                  0.1112750
#Molybdate + Nitrate-Ammonium        0.1381000
#Molybdate + Nitrate-Molybdate       0.9998846


#also not significantly different at EF

#do some diagnostic plots to look at the model fit
par(mfrow=c(2,2))
plot(anova.kmeth.bz)
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(anova.ef)
par(mfrow=c(1,1))


##report anova stuff with a plot and :
##A brief description of the variables you tested
#The f-value, degrees of freedom, and p-values for each independent variable
#What the results mean.


#plot residuals
## let's grab the residuals from our models, first BZ
bz.anova.residuals <- residuals( object = anova.kmeth.bz ) # extract the residuals
# A simple histogram
hist( x = bz.anova.residuals) # another way of seeing them
png("8.26.25.edit. BZ residuals histogram.png", units='in',width=37,height=24,res=300)
hist( x = bz.anova.residuals) 
dev.off()

#@EF
## let's grab the residuals from our model at EF
ef.anova.residuals <- residuals( object = anova.ef ) # extract the residuals
# A simple histogram
hist( x = ef.anova.residuals) 
png("8.26.25.edit. EF residuals histogram.png", units='in',width=37,height=24,res=300)
hist( x = ef.anova.residuals) # another way of seeing them
dev.off()

#Let’s also run the two most common tests of the normality assumption (as usual there are others). 
#First the Shapiro-Wilk test and then the Kolmogorov-Smirnov test (against a normal distribution).
shapiro.test( x =ef.anova.residuals ) # run Shapiro-Wilk test

#data:  ef.anova.residuals
#W = 0.9136, p-value = 0.1537
##ok cool residuals have a normal distribution
###One-sample Kolmogorov-Smirnov test

ks.test(ef.anova.residuals, "pnorm", mean(ef.anova.residuals), sd(ef.anova.residuals) )
#data:  ef.anova.residuals
#DD = 0.19188, p-value = 0.5738
#alternative hypothesis: two-sided
##ok also normal enough
##
#The final is homogeneity of variance also known as (homoscedasticity).

leveneTest(anova.ef, center = median) 
# technically a Brown Forsyth test with center as median ( a version of Levene's test)

#        Df F value Pr(>F)
#group   4  0.6955  0.612
#       10    

#interpretation:
   #p > 0.05 → Fail to reject the null hypothesis → Variances are equal (homogeneity holds).
# can assume equal variances, ok that standard ANOVA was used

#shapiro test for BZ
shapiro.test( x =bz.anova.residuals ) # run Shapiro-Wilk test
#data:  bz.anova.residuals
#W = 0.94033, p-value = 0.4226
#normally distributed 

#Kolgomorov-Smirnov test for BZ
ks.test(bz.anova.residuals, "pnorm", mean(bz.anova.residuals), sd(bz.anova.residuals) )
#data:  bz.anova.residuals
#D = 0.16045, p-value = 0.8102
#alternative hypothesis: two-sided
#normally distributed 

#The final test is homogeneity of variance also known as (homoscedasticity).
leveneTest(anova.kmeth.bz, center = median) 
# # technically a Brown Forsyth test with center as median ( a version of Levene's test)
#   Df F value  Pr(>F)  
#group    4  0.4864  0.746
#         9     
# p>0.05, can assume equal variances, ok that standard ANOVA was used

#create a theme and then produce boxplots
Nettie_Theme = theme( legend.position = "none",
                      title = NULL,
                      axis.title.x = element_text(size = 14, face = "bold"),
                      axis.text.x = element_text(size = 10, face = "bold"),
                      axis.text.y = element_text(size = 14,  face = "bold"),
                      axis.title.y = element_text(size = 14, face = "bold"), 
                      panel.grid.major = element_blank(), 
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"))



png("8.27.25 Boxplot Potential Methylation Rate Constants BZ.png", units='in',width=6,height=8,res=300)

bz.boxplot = ggplot(bz.kmeth.df, aes(x =bz.group.factor, y=kmeth.percent.per.day))+
  geom_hline(yintercept = 0.9966 ,linetype=3, col = "black", lwd = 2)+ #dotted line for median of control
  stat_boxplot(geom ='errorbar', width = 0.2) + #this supposedly plots error bars first and then we will plot the boxplot atop them??
  geom_boxplot(color = "black", fill =wes_palette("Darjeeling1"))+ #add color here and also overplot see-> # shorthand for  stat_boxplot(geom='boxplot')
  #the width specification inside stat_boxpot is to control the errorbar endcap width
  # from : https://stackoverflow.com/questions/12993545/put-whisker-ends-on-boxplot
  labs( x ="Treatment Group", y="kmeth (percent per day)") +
  Nettie_Theme
  
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        #panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title=element_text( face="bold", size=12),   axis.text=element_text(size=8, face = "bold") ,  axis.title=element_text(size=12, face = "bold")) 

dev.off()



#produce the boxplot for EF
#
png("8.27.25 Boxplot Potential Methylation Rate Constants EF.png", units='in',width=6,height=8,res=300)
ef.boxplot =ggplot(ef.kmeth.df, aes(x =ef.group.factor, y=kmeth.percent.per.day))+
  geom_hline(yintercept = 0.2758, linetype=3, col = "black", lwd = 2)+ #dotted line for median of control
  stat_boxplot(geom ='errorbar', width = 0.2) + #this supposedly plots error bars first and then we will plot the boxplot atop them??
  geom_boxplot(color = "black", fill =wes_palette("Darjeeling1"))+ #add color here and also overplot see-> # shorthand for  stat_boxplot(geom='boxplot')
  #the width specification inside stat_boxpot is to control the errorbar endcap width
  # from : https://stackoverflow.com/questions/12993545/put-whisker-ends-on-boxplot
  labs( x ="Treatment Group", y="kmeth (percent per day)") +
  Nettie_Theme

dev.off()

#group both plots together for producing a figure
figure.2020 = ggarrange(bz.boxplot, ef.boxplot, 
          labels = c("A", "B"), 
          common.legend = TRUE, 
          legend = "none",
          ncol = 2, nrow = 1)

png("8.27.25 Boxplot kmeth BZ and EF.png", units='in',width=14,height=6,res=300)
figure.2020
dev.off()

#summary stats for bz
#
bz.control.group.kmeth=bz.kmeth.df$kmeth.percent.per.day[bz.kmeth.df$group=="48 hour Control"]
bz.control.group.kmeth
median(bz.control.group.kmeth) #0.9966
mean(bz.control.group.kmeth) #1.014467
sd(bz.control.group.kmeth) #0.7997497

bz.nitrate.group.kmeth=bz.kmeth.df$kmeth.percent.per.day[bz.kmeth.df$group=="Nitrate"]
bz.nitrate.group.kmeth
median(bz.nitrate.group.kmeth) #1.9918
mean(bz.nitrate.group.kmeth) # 2.435967
sd(bz.nitrate.group.kmeth) #1.220643

bz.ammonium.group.kmeth=bz.kmeth.df$kmeth.percent.per.day[bz.kmeth.df$group=="Ammonium"]
bz.ammonium.group.kmeth
median(bz.ammonium.group.kmeth) #1.6694
mean(bz.ammonium.group.kmeth) #1.7054
sd(bz.ammonium.group.kmeth) #1.448136

bz.molybdate.group.kmeth=bz.kmeth.df$kmeth.percent.per.day[bz.kmeth.df$group=="Molybdate"]
bz.molybdate.group.kmeth
median(bz.molybdate.group.kmeth) #0.20655
mean(bz.molybdate.group.kmeth) #0.20655
sd(bz.molybdate.group.kmeth) #0.109248

bz.mandn.group.kmeth=bz.kmeth.df$kmeth.percent.per.day[bz.kmeth.df$group=="Molybdate + Nitrate"]
bz.mandn.group.kmeth
median(bz.mandn.group.kmeth) #0.2889
mean(bz.mandn.group.kmeth) #0.2665
sd(bz.mandn.group.kmeth) #0.1021587

#summary stats at East Fork
ef.control.group.kmeth=ef.kmeth.df$kmeth.percent.per.day[ef.kmeth.df$group=="48 hour Control"]
ef.control.group.kmeth
median(ef.control.group.kmeth) #0.2758
mean(ef.control.group.kmeth) #0.2333
sd(ef.control.group.kmeth) #0.1059486

ef.nitrate.group.kmeth=ef.kmeth.df$kmeth.percent.per.day[ef.kmeth.df$group=="Nitrate"]
ef.nitrate.group.kmeth
median(ef.nitrate.group.kmeth) #0.2006
mean(ef.nitrate.group.kmeth) #0.2061333
sd(ef.nitrate.group.kmeth) #0.01019624

ef.ammonium.group.kmeth=ef.kmeth.df$kmeth.percent.per.day[ef.kmeth.df$group=="Ammonium"]
ef.ammonium.group.kmeth
median(ef.ammonium.group.kmeth) #0.2785
mean(ef.ammonium.group.kmeth) #0.2422667
sd(ef.ammonium.group.kmeth) # 0.09819849

ef.molybdate.group.kmeth=ef.kmeth.df$kmeth.percent.per.day[ef.kmeth.df$group=="Molybdate"]
ef.molybdate.group.kmeth
median(ef.molybdate.group.kmeth) #0.1011
mean(ef.molybdate.group.kmeth) #0.08196667
sd(ef.molybdate.group.kmeth) #0.05926359

ef.mandn.group.kmeth=ef.kmeth.df$kmeth.percent.per.day[ef.kmeth.df$group=="Molybdate + Nitrate"]
ef.mandn.group.kmeth
median(ef.mandn.group.kmeth) #0.0907
mean(ef.mandn.group.kmeth) #0.0903
sd(ef.mandn.group.kmeth) #0.02490241

