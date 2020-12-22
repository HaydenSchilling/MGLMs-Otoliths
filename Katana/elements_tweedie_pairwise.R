# Tweedie GLM for elements Analysis

library(utils)
library(stringr)
library(mvabund)
library(tweedie)
library(statmod)


mydata <- read.csv("Otolith_data_mmol_mol_Ca.csv", header = T)
mydata$Site <-  as.factor(str_sub(as.character(mydata[,1]), start = -1)) # Get site from sample name
mydata2 <- mydata
mydata <- subset(mydata2, Site != "A")
E_data <- mydata[,c(1:14,20)]
names(E_data)[1] <- "ID"

elements <- mvabund(E_data[,c(2:8,10:14)]) # select only element data
boxplot(elements) # inspect data ranges
meanvar.plot(elements) # check mean variance relationship
plot(elements~E_data$Site) # quick rough plot of elements, x axis label wrong


#M1 <- gllvm(y = elements, formula = mydata$pop, family = "tweedie", Power=1.01, plot = T)


# Tweedie Function
fit3 <- manyany("glm", elements, data = E_data, elements ~ Site,
                family = tweedie(var.power = 1.75, link.power = 0), var.power = 1.75)
plot(fit3)
qqnorm(fit3$residuals)
# Null model for Tweedie
fitN <- manyany("glm", elements, data = E_data, elements ~ 1,
                family = tweedie(var.power = 1.75, link.power = 0), var.power = 1.75)
plot(fitN)
anova_results <- anova(fitN, fit3, p.uni = "unadjusted", nBoot = 999)
capture.output(anova_results,file="elements_anova_results_No_A.doc")

save(fit3, file = "Elements Tweedie Model_NO_A.rda")



##########################################


mydata <- read.csv("Otolith_data_mmol_mol_Ca.csv", header = T)
mydata$Site <-  as.factor(str_sub(as.character(mydata[,1]), start = -1)) # Get site from sample name
mydata2 <- mydata
mydata <- subset(mydata2, Site != "N")
E_data <- mydata[,c(1:14,20)]
names(E_data)[1] <- "ID"

elements <- mvabund(E_data[,c(2:8,10:14)]) # select only element data
boxplot(elements) # inspect data ranges
meanvar.plot(elements) # check mean variance relationship
plot(elements~E_data$Site) # quick rough plot of elements, x axis label wrong


#M1 <- gllvm(y = elements, formula = mydata$pop, family = "tweedie", Power=1.01, plot = T)


# Tweedie Function
fit3 <- manyany("glm", elements, data = E_data, elements ~ Site,
                family = tweedie(var.power = 1.75, link.power = 0), var.power = 1.75)
plot(fit3)
qqnorm(fit3$residuals)
# Null model for Tweedie
fitN <- manyany("glm", elements, data = E_data, elements ~ 1,
                family = tweedie(var.power = 1.75, link.power = 0), var.power = 1.75)
plot(fitN)
anova_results <- anova(fitN, fit3, p.uni = "unadjusted", nBoot = 999)
capture.output(anova_results,file="elements_anova_results_No_N.doc")

save(fit3, file = "Elements Tweedie Model_NO_N.rda")


######################################


mydata <- read.csv("Otolith_data_mmol_mol_Ca.csv", header = T)
mydata$Site <-  as.factor(str_sub(as.character(mydata[,1]), start = -1)) # Get site from sample name
mydata2 <- mydata
mydata <- subset(mydata2, Site != "L")
E_data <- mydata[,c(1:14,20)]
names(E_data)[1] <- "ID"

elements <- mvabund(E_data[,c(2:8,10:14)]) # select only element data
boxplot(elements) # inspect data ranges
meanvar.plot(elements) # check mean variance relationship
plot(elements~E_data$Site) # quick rough plot of elements, x axis label wrong


#M1 <- gllvm(y = elements, formula = mydata$pop, family = "tweedie", Power=1.01, plot = T)


# Tweedie Function
fit3 <- manyany("glm", elements, data = E_data, elements ~ Site, 
                family = tweedie(var.power = 1.75, link.power = 0), var.power = 1.75)
plot(fit3)
qqnorm(fit3$residuals)
# Null model for Tweedie
fitN <- manyany("glm", elements, data = E_data, elements ~ 1, 
                family = tweedie(var.power = 1.75, link.power = 0), var.power = 1.75)
plot(fitN)
anova_results <- anova(fitN, fit3, p.uni = "unadjusted", nBoot = 999)
capture.output(anova_results,file="elements_anova_results_No_L.doc")

save(fit3, file = "Elements Tweedie Model_NO_L.rda")


paste("THIS SCRIPT HAS FINISHED")
