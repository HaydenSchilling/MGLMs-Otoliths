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
families = list(tweedie(var.power = 1.75, link.power = 0))
families <- rep_len(families, length.out=12)

families2 <- list(tweedie(var.power = 2, link.power = 0))
families2 <- rep_len(families, length.out=63)

families3 <- c(families, families2)


fit3 <- manyany("glm", elements, data = mydata, elements ~ pop, 
                family = families3, var.power = c(rep.int(1.75,12), rep.int(2,63)))
# plot(fit3)
# qqnorm(fit3$residuals)
# # Null model for Tweedie
fitN <- manyany("glm", elements, data = mydata, elements ~ 1, 
                family = families3, var.power = c(rep.int(1.75,12), rep.int(2,63)))
# plot(fitN)
anova_results <- anova(fitN, fit3, p.uni = "unadjusted", nBoot = 999)
capture.output(anova_results,file="Combined_anova_results.doc")

save(fit3, file = "Combined_Tweedie_Model.rda")

paste("THIS SCRIPT HAS FINISHED")
