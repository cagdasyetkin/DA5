# Clear the console
cat("\f")

# Clean up the memory
rm(list=ls())

library(haven)
library(data.table)
library(lspline)
library(ggthemes)
library(wbstats)
library(sqldf)

library(arm)
library(readr)
library(lmtest)
library(splines)
library(stargazer)
library(sandwich)


library(wbstats)
library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(plm)  
library(ggthemes)

#RnD_exp(% of GDP) is Y
#logNumberOfArticles is X

new_cache <- wbcache()

gdp_vars <- wbsearch(pattern = "GDP")
#NY.GDP.PCAP.PP.CD GDP per capita (current US$)

pop_vars <- wbsearch(pattern = "population")
#SP.POP.TOTL

rnd_vars <- wbsearch(pattern = "r&d")
# IP.JRN.ARTC.SC Scientific and technical journal articles 

research_vars <- wbsearch(pattern = "research")
#GB.XPD.RSDV.GD.ZS Research and development expenditure (% of GDP)


data <- wb(indicator = c("NY.GDP.PCAP.PP.CD","SP.POP.TOTL", "IP.JRN.ARTC.SC", "GB.XPD.RSDV.GD.ZS"),
           return_wide = TRUE,
           POSIXct = TRUE,
           startdate = 2005,
           enddate = 2014)

write.csv(data, file = "raw_data.csv", fileEncoding = "macroman")

data <- data.table(data)

countries <- data.table(wbcountries())
countries <- countries[region != "Aggregates"]
countries <- countries$iso2c

data <- data[iso2c %in% countries]

names(data)

data[, iso2c := NULL]
data[, granularity := NULL]
data[, date_ct := NULL]

data <- setNames(data, c("country_code", "year", "country", "RnD_exp", "Sci_art", "gdp", "pop"))


data <- data[(!is.na(RnD_exp)) & (!is.na(Sci_art)) 
                  & (!is.na(pop)) & (!is.na(gdp)), ]
data[,pop := pop/10^6]
data[,country_code := factor(country_code)]

data[, gdp_cat := cut(gdp, c(0, median(gdp), Inf), labels = c("low_gdp", "high_gdp"))]




# nonmissing years for each country
# we have 66 countries with missing data
data[, COUNT := .N , by = country_code]
unique(data[COUNT!=10,list(country_code,COUNT)])
unique(data[,list(country_code,COUNT)])

# statistics from 2010
#some highly skewed - take logs!
summary(data[year == 2010, list(RnD_exp, Sci_art, gdp, pop)])  
summary(data[ , list(RnD_exp, Sci_art, gdp, pop)])  

stargazer(data, out = 'data.txt')

#histograms
ggplot(data, aes(x=RnD_exp)) +
  geom_histogram()+
  xlab('Research and development expenditure (% of GDP)')

ggplot(data, aes(x=log(Sci_art))) +
  geom_histogram()+
  xlab('Scientific and technical journal articles')

# Cross country relationship between RnD expenditure and log scientific articles
ggplot(aes(log(Sci_art), RnD_exp), data=data) + 
  geom_point(size = 2, aes(col = log(gdp))) +
  geom_smooth(method = "loess", se=F) +
  scale_color_distiller("log gdp", palette = "Spectral") +
  ylab('RnD expenditure as % of GDP') +
  xlab('ln number of scientific articles')

#repeat above plot for population
ggplot(aes(log(Sci_art), RnD_exp), data=data) + 
  geom_point(size = 2, aes(col = log(pop))) +
  geom_smooth(method = "loess", se=F) +
  scale_color_distiller("log pop", palette = "Spectral") +
  ylab('RnD expenditure as % of GDP') +
  xlab('ln number of scientific articles')


ccols1 <- lm(formula= RnD_exp ~ log(Sci_art), data=data[year==2006])
coeftest(ccols1, vcov=sandwich)


ccols2 <- lm(formula= RnD_exp ~ log(Sci_art) + log(gdp) + log(pop), data=data[year==2006])
coeftest(ccols2, vcov=sandwich)

stargazer(ccols1, ccols2, type = "text")

world_trend <- data[,lapply(mget(c("Sci_art","RnD_exp","gdp")) ,weighted.mean,w=pop), 
                      by = year]

world_trend[order(year),`:=`(rellnRnD = log(Sci_art) - first(log(Sci_art)),
                             rellngdp = log(gdp) - first(log(gdp))) ]


p1 <- ggplot(world_trend, aes(year, Sci_art)) + 
  geom_line(size = 1, color = 'darkgreen', aes(group = 1)) + 
  ggtitle('Average number of scientific articles') 

p2 <- ggplot(aes(x = year), data = world_trend) +
  geom_line(aes(group = 1, y = rellnRnD), size = 1, linetype = 'dotted', color = 'blue', show.legend = T) +
  geom_line(aes(group = 1, y = rellngdp), size = 1, linetype = 'dashed', color = 'firebrick') +
  ylab('log change from 2005') +
  ggtitle('Log change in RnD and Gdp (red colored) from 2005')

grid.arrange(p1, p2)


# Analysis 
# data[data$Sci_art == 0, .N]
# data[data$Sci_art == 0, Sci_art := NA]
panel_data <- pdata.frame(data, index = c('country', 'year'))
head(panel_data)

stargazer(panel_data, type = "text", title="Descriptive statistics", digits=1)


# First differences
diff1 <- plm(
  diff(RnD_exp) ~ diff(log(Sci_art)) + year, 
  data = panel_data, model = 'pooling'
)
#obtain clustered standard error
coeftest(diff1, vcov=vcovHC(diff1,type="HC0",cluster="group"))


diff2 <- plm(
  diff(RnD_exp) ~ diff(log(Sci_art)) + lag(diff(log(Sci_art)), 1:2) + year,
  data = panel_data, model = 'pooling'
)
coeftest(diff2, vcov=vcovHC(diff2,type="HC0",cluster="group"))


diff3 <- plm(
  diff(RnD_exp) ~ diff(log(Sci_art)) + lag(diff(log(Sci_art)), 1:4) + year, 
  data = panel_data, model = 'pooling'
)

coeftest(diff3, vcov=vcovHC(diff3,type="HC0",cluster="group"))

stargazer(diff1, diff2, diff3, type = 'html', out = 'FD.html')

# FE > use within, country fixed effect and year fixed effect > twoways
fe1 <- plm( 
  RnD_exp ~ log(Sci_art), data = panel_data, 
  model = 'within', effect = 'twoways'
)
coeftest(fe1, vcov=vcovHC(fe1,type="HC0",cluster="group"))

fe2 <- plm( 
  RnD_exp ~ log(Sci_art) + log(pop), data = panel_data, 
  model = 'within', effect = 'twoways'
)
coeftest(fe2, vcov=vcovHC(fe2,type="HC0",cluster="group"))

stargazer(fe1, fe2, type = 'html', out = 'FE.html')



