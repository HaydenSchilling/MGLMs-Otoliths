# Tweedie MGLM for Otolith Elemental Analysis

library(utils)
library(stringr)
library(mvabund)
library(tweedie)
library(statmod)
library(tidyverse)

# load the data
mydata <- read.csv("Data/Otolith_data_mmol_mol_Ca.csv", header = T)
mydata$Site <-  as.factor(str_sub(as.character(mydata[,1]), start = -1)) # Get site from sample name
E_data <- mydata[,c(1:14,20)]
names(E_data)[1] <- "ID"

elements <- mvabund(E_data[,c(2:8,10:14)]) # select only element data
boxplot(elements) # inspect data ranges
meanvar.plot(elements) # check mean variance relationship
elements_only <- elements
plot(elements~E_data$Site) # quick rough plot of elements, x axis label wrong
# fit1 <- manylm(elements ~ E_data$Site)
# summary(fit1)
# plot(fit1)


# Test PERMANOVA Appropriateness
library(vegan)

dist_euc <- vegdist(elements, method = "euclidean")

Disp <- betadisper(dist_euc, E_data$Site)
anova(Disp) #(P = 0.0003) # Therefore PERMANOVA is not appropriate (P < 0.05)

# Do not do PERMANOVA as the above test should a violation of assumptions
Per <- adonis(elements ~ E_data$Site, method = "euclidean", permutations = 9999)
Per

# nmds
## NMDS

library(MASS)
fit <- isoMDS(dist_euc, k=2) # k is the number of dim
fit # view results

# plot solution
x <- fit$points[,1]
y <- fit$points[,2]
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


points_to_plot_elements_MDS <- mydata %>% dplyr::select(NMDS_x, NMDS_y, Site) %>% mutate(Data = "a) Chemistry Data")
head(points_to_plot_elements_MDS)
write_csv(points_to_plot_elements_MDS, "Data/Elements NMDS Ordination.csv")

# Tweedie Function
# do model selection to set var.power

fit3 <- manyany("glm", elements, data = E_data, elements ~ Site, 
                family = tweedie(var.power = 1.75, link.power = 0), var.power = 1.75)
plot(fit3)

qqnorm(residuals.manyany(fit3))
qqline(residuals.manyany(fit3))

fit_e <- fit3 # for plotting later

# Null model for Tweedie
fitN <- manyany("glm", elements, data = E_data, elements ~ 1, 
                family = tweedie(var.power = 1.75, link.power = 0), var.power = 1.75)
plot(fitN)
anova_results <- anova(fitN, fit3, p.uni = "unadjusted", nBoot = 999, nCores = 3) # This could be very slow
capture.output(anova_results,file="elements_anova_results.doc")

save(fit3, file = "../Data/Elements Tweedie Model.rda")

paste("THIS SCRIPT HAS FINISHED")

# Do Ordination

library(ecoCopula)
element_LV=cord(fitN)
plot(element_LV,biplot = TRUE,site.col = E_data$Site)
plot(element_LV,site.col = E_data$Site)
plot(fitN)
qqnorm(residuals(fitN))
qqline(residuals(fitN))


library(gllvm)
fitX <- gllvm(y = elements, family = "tweedie", Power = 1.9)
ordiplot.gllvm(fitX)


### plots of elements
mydata_long <- mydata %>% dplyr::select(-Ca) %>% pivot_longer(cols=c(Ni:Ba), values_to="Ratio", names_to="Element")
head(mydata_long)

mydata_long_sum <- mydata_long %>% group_by(Site, Element) %>% summarise(Mean_ratio = mean(Ratio), SD_Ratio = sd(Ratio), n=n(), SE_Ratio = SD_Ratio/sqrt(n))

ggplot(mydata_long_sum, aes(x=Site, y = Mean_ratio)) + geom_bar(stat="identity", fill="grey60") +
  facet_wrap(~Element, scales="free_y") + theme_classic()+
  ylab(bquote(bold("Element:Ca (mmol mol"^-1*")")))+
  geom_errorbar(aes(ymin=Mean_ratio-SE_Ratio, ymax=Mean_ratio + SE_Ratio),width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  scale_x_discrete(breaks=c("A", "L", "N"), labels=c("Agra", "Lucknow", "Narora"))+
  theme(axis.title = element_text(face="bold", size=12),
        axis.text.y = element_text(colour="black", size = 10),
        axis.text.x = element_text(colour="black", size = 10, angle=45, hjust = 1),
        strip.text = element_text(size= 12, face="bold"))

ggsave("Figures/Elements.pdf", height = 15, width = 18, units ="cm", dpi = 600)
ggsave("Figures/Elements.png", height = 15, width = 18, units ="cm", dpi = 600)
