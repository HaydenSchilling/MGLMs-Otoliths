# Tweedie GLM for Shape Analysis

library(utils)
library(mvabund)
library(tweedie)
library(statmod)

mydata <- read.csv("shape_data.csv", header = T)

elements <- mvabund(mydata[,c(13:75)]) # select only shape data
boxplot(elements) # inspect data ranges
meanvar.plot(elements) # check mean variance relationship
plot(elements~mydata$pop) # quick rough plot of elements, x axis label wrong
fit1 <- manylm(elements ~ mydata$pop)
summary(fit1)
plot(fit1)

#M1 <- gllvm(y = elements, formula = mydata$pop, family = "tweedie", Power=1.01, plot = T)

# Tweedie Function
fit3 <- manyany("glm", elements, data = mydata, elements ~ pop, 
                family = tweedie(var.power = 2, link.power = 0), var.power = 2)
# plot(fit3)
# qqnorm(fit3$residuals)
# # Null model for Tweedie
 fitN <- manyany("glm", elements, data = mydata, elements ~ 1, 
                 tweedie(var.power = 2, link.power = 0), var.power = 2)
# plot(fitN)
 anova_results <- anova(fitN, fit3, p.uni = "unadjusted", nBoot = 999)
 capture.output(anova_results,file="shape_anova_results.doc")

save(fit3, file = "Shape_Tweedie_Model.rda")

paste("THIS SCRIPT HAS FINISHED")

