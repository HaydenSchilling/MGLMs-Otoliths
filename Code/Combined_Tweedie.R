# Tweedie GLM for Otolith Shape and Element combined Analysis
# This was run on a HPC

library(utils)
library(mvabund)
library(tweedie)
library(statmod)
library(DHARMa)

# Load the data
mydata <- read.csv("Data/Otolith_data_mmol_mol_Ca_and_shape.csv", header = T)
mydata$Site <-  as.factor(str_sub(as.character(mydata[,1]), start = -1)) # Get site from sample name
E_data <- mydata[,c(3:77,83)]
names(E_data)[76] <- "ID"
elements <- mvabund(mydata[,c(3:77)]) # select data
#elements_only <- mvabund(mydata[,c(3:14)]) # select data

boxplot(elements) # inspect data ranges
meanvar.plot(elements) # check mean variance relationship
elements_and_shape <- elements
plot(elements~E_data$ID) # quick rough plot of elements, x axis label wrong
#fit1 <- manylm(elements ~ mydata$pop)
#summary(fit1)
#plot(fit1)

# Test PERMANOVA Appropriateness
library(vegan)

dist_euc <- vegdist(elements, method = "euclidean")

Disp <- betadisper(dist_euc, mydata$Site)
anova(Disp) #(P = 0.0004) # Therefore assumption for PERMANOVA is NOT met 

# Now run PERMANOVA
Per <- adonis(dist_euc ~ E_data$ID, method = "euclidean", permutations = 9999)
Per # Significant group differences (P <0.0001)



# nmds
library(MASS)
fitO <- isoMDS(dist_euc, k=2) # k is the number of dim
fitO # view results

# plot solution
x <- fitO$points[,1]
y <- fitO$points[,2]
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
  theme_classic() +# xlim(c(-2.1,2.1)) +ylim(c(-2.1,2.1))+
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

points_to_plot_combined_MDS <- mydata %>% dplyr::select(NMDS_x, NMDS_y, Site) %>% mutate(Data = "c) Combined Data")
head(points_to_plot_combined_MDS)
write_csv(points_to_plot_combined_MDS, "Data/Combined NMDS Ordination.csv")




# Tweedie Function
families = list(tweedie(var.power = 1.75, link.power = 0))
families <- rep_len(families, length.out=12)

families2 <- list(tweedie(var.power = 2, link.power = 0))
families2 <- rep_len(families, length.out=63)

families3 <- c(families, families2)


fit3 <- manyany("glm", elements, data = mydata, elements ~ pop, 
                family = families3, var.power = c(rep.int(1.75,12), rep.int(2,63)))
plot(fit3)
#plot(fit3, log="x")

qqnorm(fit3$residuals)
qqline(fit3$residuals)
fit_c <- fit3

# # Null model for Tweedie
fitN <- manyany("glm", elements, data = mydata, elements ~ 1, 
                 family = families3, var.power = c(rep.int(1.75,12), rep.int(2,63)))
plot(fitN)
anova_results <- anova(fitN, fit3, p.uni = "unadjusted", nBoot = 999) # this will be very slow
capture.output(anova_results,file="Combined_anova_results.doc")

save(fit3, file = "../Data/Combined_Tweedie_Model.rda")

paste("THIS SCRIPT HAS FINISHED")


# Do Ordination

library(ecoCopula)
element_LV=cord(fitN, n.samp = 50)
plot(element_LV,biplot = TRUE,site.col = E_data$ID)
plot(element_LV,site.col = E_data$ID)
plot(fitN)
qqnorm(residuals(fitN))
qqline(residuals(fitN))

library(gllvm)
fitX <- gllvm(y = elements, family = "tweedie", Power = 1.9)
ordiplot.gllvm(fitX)
