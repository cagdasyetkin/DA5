################################################################################
# HEALTH EXPENDITURES AND LIFE EXPECTANCY
################################################################################

library(ggplot2)
library(gridExtra)  # to combine graphs
library(plm)  # for panels
library(data.table)
library(lmtest)
library(sandwich)

#setwd("/home/zsuzsa/Documents/DA5")


# Load data --------------------------------------------------------------------

health <- fread('22_WB_data_on_health.csv')

original_names <- names(health)
names(health) <- c(
    'year', 'year_code', 'country_name', 'country_code', 
    'hexppc', 'gdppc', 'gdp_growth', 
    'lifeexp', 'infmortality', 'u5mortality', 'pop'
)

#gdppc "GDP per capita (PPP constant 2005 $)"
#lifeexp "Life expectancy at birth"
#hexppc "Health expenditure per capita (PPP constant 2005 $)"
#pop "Population"

health <- health[ (year < 2011) & (!is.na(hexppc)) & (!is.na(lifeexp)) 
        & (!is.na(pop)) & (!is.na(gdppc)) ]
health[,pop := pop/10^6]
health[,country_code := factor(country_code)]

# Know your data ---------------------------------------------------------------

# nonmissing years for each country
#160 countries with complete data
health[, COUNT := .N , by = country_code]
unique(health[COUNT!=16,list(country_code,COUNT)])

# statistics from 2010
#all highly skewed except life expectancy - take logs
summary(health[year == 2010, list(hexppc, lifeexp, gdppc, pop)])

#histograms
ggplot(aes(x=lifeexp), data=health) +
  geom_histogram(bins = 35)+
  xlab('Life expectancy at birth')

ggplot(aes(x=log(hexppc)), data=health) +
  geom_histogram(bins = 40)+
  xlab('ln health expenditure per capita (PPP constant 2005$)')

# Cross country OLS between life expectancy and log health expenditures per capita
ggplot(aes(log(hexppc),lifeexp ), data=health[year == 2003] ) + 
  geom_point(size = 3, color = 'orange', shape=4 ) +
  geom_smooth(method = "lm",se=FALSE  ) +
  ylab('Life expectancy at birth') +
  xlab('ln health expenditure per capita (PPP constant 2005$)')

ccols1 <- lm(formula= lifeexp~ log(hexppc), data=health[year==2003])
coeftest(ccols1, vcov=sandwich)

ccols2 <- lm(formula= lifeexp~ log(hexppc) + log(gdppc) + log(pop), data=health[year==2003])
coeftest(ccols2, vcov=sandwich)

#World trends of life expectancy, health expenditures per capita GDP per capita 
#Avg of all countries weighted by population, normalized by subtracting 1995 values

world_trend <- health[,lapply(mget(c("hexppc","lifeexp","gdppc")) ,weighted.mean,w=pop), 
         by = year]

world_trend[order(year),`:=`(rellnhexppc = log(hexppc) - first(log(hexppc)),
                  rellngdp = log(gdppc) - first(log(gdppc))) ]


p1 <- ggplot(aes(year, lifeexp), data = world_trend) + 
    geom_line(size = 1, color = 'darkgreen') + 
    ylab('average life expectancy')

p2 <- ggplot(aes(x = year), data = world_trend) +
    geom_line(aes(y = rellnhexppc), size = 1, linetype = 'dotted', color = 'blue') +
    geom_line(aes(y = rellngdp), size = 1, linetype = 'dashed', color = 'firebrick') +
    geom_text(x = 2009, y = 0.7, label = 'health expenditure', color = 'blue') +
    geom_text(x = 2009, y = 0.25, label = 'GDP', color = 'firebrick') +
    ylab('log change from 1995')

grid.arrange(p1, p2)

# Analysis ---------------------------------------------------------------------

health[health$hexppc==0, hexppc := NA]
panel_health <- pdata.frame(health, index = c('country_name', 'year'))
head(panel_health)

# First differences
diff1 <- plm(
    diff(lifeexp) ~ diff(log(hexppc)) + year, 
    data = panel_health, model = 'pooling'
)
#obtain clustered standard error
coeftest(diff1, vcov=vcovHC(diff1,type="HC0",cluster="group"))
#logs- if we have 10% increase in health exp in a given country the increase should be 10% expen above overall change accross countries
#0.01 year increase in life exp in that country.
#control for year. log diff on right hand side. divide by 100.
#additional to this average yearly change, how much it will increase. annual trend accross countries are the dummies.
#accross countries and all years.

#this doesnt take the lag into account. expend will effect 5 years later. not now. we will do it in the following part:

diff2 <- plm(
    diff(lifeexp) ~ diff(log(hexppc)) + lag(diff(log(hexppc)), 1:2) + year,
    data = panel_health, model = 'pooling'
)
coeftest(diff2, vcov=vcovHC(diff2,type="HC0",cluster="group"))
#2 lags, original veriable not significant any more. but the second one is significant.

diff3 <- plm(
    diff(lifeexp) ~ diff(log(hexppc)) + lag(diff(log(hexppc)), 1:4) + year, 
    data = panel_health, model = 'pooling'
)
coeftest(diff3, vcov=vcovHC(diff3,type="HC0",cluster="group"))
diff4 <- plm(
    diff(lifeexp) ~ diff(log(hexppc)) +  lag(diff(log(hexppc)), 1:6) + year,
    data = panel_health, model = 'pooling'
)
coeftest(diff4, vcov=vcovHC(diff4,type="HC0",cluster="group"))
#second adn third lags in all the cases. we have a lagged effect. for 2-3 years. which lag to include?
#we will try an approach to decide that.

#How many lags? Focus on cumulative associations
#Sum the lagged coeffs and stop when that remains same

diff5 <- plm(diff(lifeexp) ~ lag(diff(log(hexppc)), 2) + lag(diff(diff(log(hexppc))), 0:1) + year,
    data = panel_health, model = 'pooling')
coeftest(diff5, vcov=vcovHC(diff5,type="HC0",cluster="group"))
#obtain a significance level for all the lagged effect. cummulatie as one. instead of haingf.
#second lag of the first diff. then we have second diff of original veriable.
#first variable will have the same coeff 0.355 sum of coefs in previous model with first lags.
#here we wanna have a cummulative effect. 
summary(diff2) #check it out. see 0.35 - 0.18 gives coef of the previous one.
#result of this model is 0.36 and it is significant at 5% level


diff6 <- plm(diff(lifeexp) ~ lag(diff(log(hexppc)), 4) + lag(diff(diff(log(hexppc))), 0:3) + year,
    data = panel_health, model = 'pooling')
coeftest(diff6, vcov=vcovHC(diff6,type="HC0",cluster="group"))
diff7 <- plm(diff(lifeexp) ~ lag(diff(log(hexppc)), 6) + lag(diff(diff(log(hexppc))), 0:5) + year,
    data = panel_health, model = 'pooling')
coeftest(diff7, vcov=vcovHC(diff7,type="HC0",cluster="group"))
diff8 <- plm(diff(lifeexp) ~ lag(diff(log(hexppc)), 8) + lag(diff(diff(log(hexppc))), 0:7) + year,
    data = panel_health, model = 'pooling')
coeftest(diff8, vcov=vcovHC(diff8,type="HC0",cluster="group"))

#4 lag: cumm eff. 0.509 is signif at 5%level
#6 lag still significant
#8 lags lost significance. this additional lafgs re not needed. conclude that upto 6 we have
#but after 6 it doest, so stop at 6 and use it.



# Cumulative associations with controls
diff4_1 <- plm(
    diff(lifeexp) ~ lag(diff(log(hexppc)), 6) + 
        lag(diff(log(gdppc)), 6) + 
        lag(diff(log(pop)), 6) +
        lag(diff(diff(log(hexppc))), 0:5) + 
        lag(diff(diff(log(gdppc))), 0:5) + 
        lag(diff(diff(log(pop))), 0:5) + year, 
    data = panel_health, model = 'pooling'
)
coeftest(diff4_1, vcov=vcovHC(diff4_1,type="HC0",cluster="group"))
#control for gdp. if the significance is decreasing or not. include population for technical reasons.
#because it effects both of your heath exp and RH. include pop. avoid effects which might both life exp and health exp.
#because pop effects. added with the same number of lags. 
#health exp still significant in the result. pop has high significant coeff here. we wanna inclide that.
#but gdp doesnt seem to have an effect on life exp. including gdp doesnt really change the reuslt.
#so we can exclude that.

diff4_2 <- plm(
    diff(lifeexp) ~ lag(diff(log(hexppc)), 6) + 
        lag(diff(log(pop)), 6) +
        lag(diff(diff(log(hexppc))), 0:5) + 
        lag(diff(diff(log(pop))), 0:5) + year, 
    data = panel_health, model = 'pooling'
)
coeftest(diff4_2, vcov=vcovHC(diff4_2,type="HC0",cluster="group"))
#precision increased.

# FE > use within, country fixed effect and year fixed effect > twoways
fe1 <- plm( 
    lifeexp ~ log(hexppc), data = panel_health, 
    model = 'within', effect = 'twoways'
)
coeftest(fe1, vcov=vcovHC(fe1,type="HC0",cluster="group"))
#cant see the country and year dummies here in the result. but similar interpretation.
#diff in fe and fd estimation its a bit different. here we have,
#if you have between countries or between years 10% higher h exp in a given country
#will increase the life ext by 0.07 years compared to the average in life exp for that country.
#its always compared to the average for that given coutnry. basically, estimating a fe model
#is the same like OLS. effect is always compared to the mean of that country. demeaned OLS.

fe2 <- plm( 
    lifeexp ~ log(hexppc) + log(pop), data = panel_health, 
    model = 'within', effect = 'twoways'
)
coeftest(fe2, vcov=vcovHC(fe2,type="HC0",cluster="group"))
#only we need to specifiy, having the same population, you compare the expenditure.
#no laggs anymore. we dont include lags. why. because, intuitively, using fixed effect already takes it into account
#think of OLS demeaned regression. demeaned, overall change in health expenditure and life exp. manual control not needed.
#never include laggs in fe models. 


# FE with explicit time dummies -same results
fe1b <- plm( 
    lifeexp ~ log(hexppc) + year, data = panel_health, 
    model = 'within'
)
fe2b <- plm( 
    lifeexp ~ log(hexppc) + log(pop) + year, data = panel_health, 
    model = 'within'
)

#Long run . long run change in life exp
panel_health2 <- pdata.frame(health[(year==2010)|(year==1995)], #first and the last year
                             index = c('country_name', 'year'))

lr1 <- plm(
  diff(lifeexp) ~ diff(log(hexppc)) + diff(log(pop)), 
  data = panel_health2, model = 'pooling'
)
coeftest(lr1, vcov=vcovHC(lr1,method="white1"))
#fe model interpret but here we say thw whole years. 10% change in health exp, expect to have 0.2 years increase in life exp
#during this 16 years. white1: heteroscastecity se.

# Robustness checks
#dentified 8 extreme value countries rrom Africa; civil wars and/or epidemic
outliers <- c('ERI', 'GNB', 'RWA', 'LBR', 'GNQ', 'LSO', 'BWA', 'ZAF', 'SWZ')

panel_health_no_outlier <- pdata.frame(health[!(country_code %in% outliers)],
                                       index = c('country_name', 'year'))

#compare to: diff3_4
robust1 <- plm(
    diff(lifeexp) ~ lag(diff(log(hexppc)), 6) + 
        lag(diff(log(pop)), 6) +
        lag(diff(diff(log(hexppc))), 0:5) + 
        lag(diff(diff(log(pop))), 0:5) + year, 
    data = panel_health_no_outlier, model = 'pooling'
)
coeftest(robust1, vcov=vcovHC(robust1,type="HC0",cluster="group"))
#here we have significant coeff. to sum up the results are similar. so africa is not really distording at all.
#


#Heterogeneity poor vs rich countries
health[,`:=`(mean_pop = mean(pop), mean_gdp = mean(gdppc)), by=country_name ]
health[,`:=`(gdp_cat = ifelse(mean_gdp < median(mean_gdp), 'poor', 'rich')) ]
panel_health <- pdata.frame(health, index = c('country_name', 'year'))

robust2 <- plm(
    diff(lifeexp) ~ lag(diff(log(hexppc)), 6) + 
        lag(diff(log(pop)), 6) +
        lag(diff(diff(log(hexppc))), 0:5) + 
        lag(diff(diff(log(pop))), 0:5) + year, 
    data = panel_health, subset = (gdp_cat == 'poor'), model = 'pooling'
)
coeftest(robust2, vcov=vcovHC(robust2,type="HC0",cluster="group"))
#heath exp is not significant and coeff is lower. 

robust3 <- plm(
    diff(lifeexp) ~ lag(diff(log(hexppc)), 6) + 
        lag(diff(log(pop)), 6) +
        lag(diff(diff(log(hexppc))), 0:5) + 
        lag(diff(diff(log(pop))), 0:5) + year, 
    data = panel_health, subset = (gdp_cat == 'rich'), model = 'pooling'
)
coeftest(robust3, vcov=vcovHC(robust3,type="HC0",cluster="group"))
#the opposite. higher coeff higher significance. seems that effect is stronger for rich countries. after controlling for population.
#poor countries can be distorting because of bad data collection also. richer countries on average benefit more from increased health expenditures than poor countres.

#pooled OLS
pols <- lm(lifeexp ~ log(hexppc) + log(pop) + factor(year), data=health)
coeftest(pols, vcov=sandwich)
#just for exploration code. 

#use fd if you see unit root. this can be used only to balanced panels.
#if you reject the test, then it means it is statinary. argument for using fd model.
