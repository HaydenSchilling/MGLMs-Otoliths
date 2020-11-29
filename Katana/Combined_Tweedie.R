# Tweedie GLM for Shape and element combined Analysis

library(utils)
library(mvabund)
library(tweedie)
library(statmod)

mydata <- read.csv("Otolith_data_mmol_mol_Ca_and_shape.csv", header = T)

elements <- mvabund(mydata[,c(3:77)]) # select only shape data
boxplot(elements) # inspect data ranges
meanvar.plot(elements) # check mean variance relationship
plot(elements~mydata$pop) # quick rough plot of elements, x axis label wrong
fit1 <- manylm(elements ~ mydata$pop)
summary(fit1)
plot(fit1)

# Tweedie Function
fit3 <- manyany("glm", elements, data = mydata, elements ~ pop, 
                 family = tweedie(var.power = 1.9), var.power = 1.9)
# plot(fit3)
# qqnorm(fit3$residuals)
# # Null model for Tweedie
fitN <- manyany("glm", elements, data = mydata, elements ~ 1, 
               family = tweedie(var.power = 1.9), var.power = 1.9)
# plot(fitN)
anova_results <- anova(fitN, fit3, p.uni = "unadjusted", nBoot = 9999)
capture.output(anova_results,file="Combined_anova_results.doc")

save(fit3, file = "Combined_Tweedie_Model.rda")

paste("THIS SCRIPT HAS FINISHED")
