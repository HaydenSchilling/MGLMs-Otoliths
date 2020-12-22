# Tweedie GLM for Otolith Shape Analysis
# This was run on a HPC

library(utils)
library(mvabund)
library(tweedie)
library(statmod)

# Load the data
mydata <- read.csv("Data/Otolith_data_mmol_mol_Ca_and_shape.csv", header = T)
mydata$Site <-  as.factor(str_sub(as.character(mydata[,1]), start = -1)) # Get site from sample name
E_data <- mydata[,c(15:77,83)]
names(E_data)[64] <- "ID"

elements <- mvabund(E_data[,c(1:63)]) # select only shape data
boxplot(elements) # inspect data ranges
meanvar.plot(elements) # check mean variance relationship
shape_only <- elements
#plot(elements~mydata$pop) # quick rough plot of elements, x axis label wrong
#fit1 <- manylm(elements ~ mydata$pop)
#summary(fit1)
#plot(fit1)

# Test PERMANOVA Appropriateness
library(vegan)

dist_euc <- vegdist(elements, method = "euclidean")

Disp <- betadisper(dist_euc, mydata$Site)
anova(Disp) #(P = 0.659) # Therefore assumption for PERMANOVA is met 

# Now run PERMANOVA
Per <- adonis(dist_euc ~ E_data$ID, permutations = 9999)
Per # Significant group differences (P <0.0001)

## NMDS

library(MASS)
fit2 <- isoMDS(dist_euc, k=2) # k is the number of dim
fit2 # view results

# plot solution
x <- fit2$points[,1]
y <- fit2$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Nonmetric MDS", type="n")
text(x, y, labels = E_data$ID, cex=.7)
# shows some separation but also large differencs in dispersion
# shows some separation but not a lot

# do ggplot
mydata$NMDS_x <- x
mydata$NMDS_y <- y



library(ggplot2)
ggplot(mydata, aes(x=x, y=y, col = Site, shape = Site)) + geom_point() +
  theme_classic() + xlim(c(-2.1,2.1)) +ylim(c(-2.1,2.1))+
  scale_color_manual(values = c("black", "red", "blue"), name = "Population",
                     breaks=c("A", "L", "N"),
                     labels=c("Agra", "Lucknow", "Narora")) +
  scale_shape_manual(values = c("circle", "triangle", "square"), name = "Population",
                     breaks=c("A", "L", "N"),
                     labels=c("Agra", "Lucknow", "Narora"))+
  ylab("NMDS 2") + xlab("NMDS 1") +
  theme(axis.title = element_text(face="bold", size = 14),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(face = "bold", size = 10),
        legend.position = c(0.8,0.2),
        legend.background = element_rect(colour="black"))


points_to_plot_shape_MDS <- mydata %>% dplyr::select(NMDS_x, NMDS_y, Site) %>% mutate(Data = "b) Shape Data")
head(points_to_plot_shape_MDS)
write_csv(points_to_plot_shape_MDS, "Data/Shape NMDS Ordination.csv")


#ggsave("plots/Shape NMDS.png", dpi = 600, width =12, height = 12, units = "cm")

# Tweedie Function
fit3 <- manyany("glm", elements, data = mydata, elements ~ pop, 
                family = tweedie(var.power = 2, link.power = 0), var.power = 2)

plot(fit3)
qqnorm(residuals.manyany(fit3))
qqline(residuals.manyany(fit3))

fit_s <- fit3

# # Null model for Tweedie
fitN <- manyany("glm", elements, data = mydata, elements ~ 1, 
                 family = tweedie(var.power = 2, link.power = 0), var.power = 2)

anova_results <- anova(fitN, fit3, p.uni = "unadjusted", nBoot = 999)
capture.output(anova_results,file="shape_anova_results.doc")

save(fit3, file = "../Data/Shape_Tweedie_Model.rda")

paste("THIS SCRIPT HAS FINISHED")


### Try pairwise by elimination

library(utils)
library(mvabund)
library(tweedie)
library(statmod)
library(stringr)

# Load the data
mydata <- read.csv("Data/Otolith_data_mmol_mol_Ca_and_shape.csv", header = T)

mydata <- subset(mydata, pop != "LL")
mydata$Site <-  as.factor(str_sub(as.character(mydata[,1]), start = -1)) # Get site from sample name
E_data <- mydata[,c(15:77,83)]
names(E_data)[64] <- "ID"

elements <- mvabund(E_data[,c(1:63)]) # select only shape data
boxplot(elements) # inspect data ranges
meanvar.plot(elements) # check mean variance relationship
shape_only <- elements
#plot(elements~mydata$pop) # quick rough plot of elements, x axis label wrong
#fit1 <- manylm(elements ~ mydata$pop)


# Tweedie Function
fit3 <- manyany("glm", elements, data = mydata, elements ~ pop, 
                family = tweedie(var.power = 2, link.power = 0), var.power = 2)

plot(fit3)
qqnorm(residuals.manyany(fit3))
qqline(residuals.manyany(fit3))

# # Null model for Tweedie
fitN <- manyany("glm", elements, data = mydata, elements ~ 1, 
                family = tweedie(var.power = 2, link.power=0), var.power = 2)

anova_results <- anova(fitN, fit3, p.uni = "unadjusted", nBoot = 999)
capture.output(anova_results,file="shape_anova_results NO L.doc")
anova_results
save(fit3, file = "Data/Shape_Tweedie_Model_NO_L.rda")

paste("THIS SCRIPT HAS FINISHED")


### Now NO A

mydata <- read.csv("Data/Otolith_data_mmol_mol_Ca_and_shape.csv", header = T)

mydata <- subset(mydata, pop != "AA")
mydata$Site <-  as.factor(str_sub(as.character(mydata[,1]), start = -1)) # Get site from sample name
E_data <- mydata[,c(15:77,83)]
names(E_data)[64] <- "ID"

elements <- mvabund(E_data[,c(1:63)]) # select only shape data
boxplot(elements) # inspect data ranges
meanvar.plot(elements) # check mean variance relationship
shape_only <- elements
#plot(elements~mydata$pop) # quick rough plot of elements, x axis label wrong
#fit1 <- manylm(elements ~ mydata$pop)


# Tweedie Function
fit3 <- manyany("glm", elements, data = mydata, elements ~ pop, 
                family = tweedie(var.power = 2, link.power = 0), var.power = 2)

plot(fit3)
qqnorm(residuals.manyany(fit3))
qqline(residuals.manyany(fit3))

# # Null model for Tweedie
fitN <- manyany("glm", elements, data = mydata, elements ~ 1, 
                family = tweedie(var.power = 2, link.power=0), var.power = 2)

anova_results <- anova(fitN, fit3, p.uni = "unadjusted", nBoot = 999)
capture.output(anova_results,file="shape_anova_results NO A.doc")

save(fit3, file = "Data/Shape_Tweedie_Model_NO_A.rda")

paste("THIS SCRIPT HAS FINISHED")

### NOW NO N

mydata <- read.csv("Data/Otolith_data_mmol_mol_Ca_and_shape.csv", header = T)

mydata <- subset(mydata, pop != "NN")
mydata$Site <-  as.factor(str_sub(as.character(mydata[,1]), start = -1)) # Get site from sample name
E_data <- mydata[,c(15:77,83)]
names(E_data)[64] <- "ID"

elements <- mvabund(E_data[,c(1:63)]) # select only shape data
boxplot(elements) # inspect data ranges
meanvar.plot(elements) # check mean variance relationship
shape_only <- elements
#plot(elements~mydata$pop) # quick rough plot of elements, x axis label wrong
#fit1 <- manylm(elements ~ mydata$pop)


# Tweedie Function
fit3 <- manyany("glm", elements, data = mydata, elements ~ pop, 
                family = tweedie(var.power = 2, link.power = 0), var.power = 2)

plot(fit3)
qqnorm(residuals.manyany(fit3))
qqline(residuals.manyany(fit3))

# # Null model for Tweedie
fitN <- manyany("glm", elements, data = mydata, elements ~ 1, 
                family = tweedie(var.power = 2, link.power=0), var.power = 2)

anova_results <- anova(fitN, fit3, p.uni = "unadjusted", nBoot = 999) 
capture.output(anova_results,file="shape_anova_results NO N.doc")

save(fit3, file = "Data/Shape_Tweedie_Model_NO_N.rda")

paste("THIS SCRIPT HAS FINISHED")