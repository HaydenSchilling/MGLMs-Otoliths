# Tweedie MGLM for Otolith Elemental Analysis

library(utils)
library(stringr)
library(mvabund)
library(tweedie)
library(statmod)

# load the data
mydata <- read.csv("../Data/Otolith_data_mmol_mol_Ca.csv", header = T)
mydata$Site <-  as.factor(str_sub(as.character(mydata[,1]), start = -1)) # Get site from sample name
E_data <- mydata[,c(1:14,20)]
names(E_data)[1] <- "ID"

elements <- mvabund(E_data[,c(2:8,10:14)]) # select only element data
boxplot(elements) # inspect data ranges
meanvar.plot(elements) # check mean variance relationship
plot(elements~E_data$Site) # quick rough plot of elements, x axis label wrong
# fit1 <- manylm(elements ~ E_data$Site)
# summary(fit1)
# plot(fit1)

# Tweedie Function
fit3 <- manyany("glm", elements, data = E_data, elements ~ Site, 
                family = tweedie(var.power = 1.01), var.power = 1.01)
plot(fit3)
qqnorm(fit3$residuals)

# Null model for Tweedie
fitN <- manyany("glm", elements, data = E_data, elements ~ 1, 
                family = tweedie(var.power = 1.01), var.power = 1.01)
plot(fitN)
anova_results <- anova(fitN, fit3, p.uni = "unadjusted", nBoot = 9999) # This could be very slow
capture.output(anova_results,file="elements_anova_results.doc")

save(fit3, file = "../Data/Elements Tweedie Model.rda")

paste("THIS SCRIPT HAS FINISHED")
