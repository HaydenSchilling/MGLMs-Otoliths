# Tweedie GLM for Otolith Shape and Element combined Analysis
# This was run on a HPC

library(utils)
library(mvabund)
library(tweedie)
library(statmod)
library(DHARMa)

# Load the data
mydata <- read.csv("../Data/Otolith_data_mmol_mol_Ca_and_shape.csv", header = T)

elements <- mvabund(mydata[,c(3:77)]) # select data
#elements_only <- mvabund(mydata[,c(3:14)]) # select data

boxplot(elements) # inspect data ranges
meanvar.plot(elements) # check mean variance relationship
plot(elements~mydata$pop) # quick rough plot of elements, x axis label wrong
#fit1 <- manylm(elements ~ mydata$pop)
#summary(fit1)
#plot(fit1)

# Tweedie Function
fit3 <- manyany("glm", elements, data = mydata, elements ~ pop, 
                 family = tweedie(var.power = 1.01), var.power = 1.01)
plot(fit3)
qqnorm(fit3$residuals)
qqline(fit3$residuals)
# # Null model for Tweedie
fitN <- manyany("glm", elements, data = mydata, elements ~ 1, 
                 family = tweedie(var.power = 1.01), var.power = 1.01)
plot(fitN)
anova_results <- anova(fitN, fit3, p.uni = "unadjusted", nBoot = 9999) # this will be very slow
capture.output(anova_results,file="Combined_anova_results.doc")

save(fit3, file = "../Data/Combined_Tweedie_Model.rda")

paste("THIS SCRIPT HAS FINISHED")


