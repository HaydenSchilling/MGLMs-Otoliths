# make mean variance plots
par(mfrow=c(1,3))


library(utils)
library(mvabund)
library(tweedie)
library(statmod)
library(gllvm)


mydata <- read.csv("Otolith_data_mmol_mol_Ca.csv", header = T)
mydata$Site <-  as.factor(str_sub(as.character(mydata[,1]), start = -1)) # Get site from sample name
E_data <- mydata[,c(1:14,20)]
names(E_data)[1] <- "ID"

elements1 <- mvabund(E_data[,c(2:8,10:14)]) # select only element data
#boxplot(elements) # inspect data ranges
meanvar.plot(elements1) # check mean variance relationship

mydata <- read.csv("shape_data.csv", header = T)

elements2 <- mvabund(mydata[,c(13:75)]) # select only shape data
meanvar.plot(elements2) # check mean variance relationship


mydata <- read.csv("Otolith_data_mmol_mol_Ca_and_shape.csv", header = T)

elements3 <- mvabund(mydata[,c(3:77)]) # select only shape data
#boxplot(elements) # inspect data ranges
meanvar.plot(elements3) # check mean variance relationship

# exported to illustrator to put in a single figure panel
