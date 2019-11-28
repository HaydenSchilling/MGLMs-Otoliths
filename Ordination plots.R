# Tweedie ordination for elements Analysis

# Load packages and data
#install.packages("mvabund")
library(mvabund)
# 
# data(spider) # example data from mvabund
# 
#install.packages("boral")
library(boral)
# y <- spider$abund
# fit.lvmp <- boral(y = y, family = "poisson", num.lv = 2, row.eff = "fixed")
# summary(fit.lvmp)
# fit.lvmp$hpdintervals
# plot(fit.lvmp)
# 
# fit.lvmnb <- boral(y = y, family = "negative.binomial", num.lv = 2, row.eff = "fixed", offset = 1)
# summary(fit.lvmnb)
# fit.lvmnb$hpdintervals
# plot(fit.lvmnb)
# 
# lvsplot(fit.lvmnb, biplot = FALSE)


## COmbined shape and element
# load data
Fish = read.csv("Otolith_data_mmol_mol_Ca_and_shape.csv", header =T)

str(Fish)
# Convert fish to integers if needed # not sure if needed
#for(n in names(Fish)[1:38]){
#  Fish[[n]]<-as.integer(Fish[[n]])}


# Separate out the columns into important bits
Location = as.factor(Fish[,2])
Fish_community = as.mvabund(Fish[,c(3:77)]) #Select species (<3.4mm ESD)
#Flow = Fish[,3]
#DepthLocation <- as.factor(paste(Location, Depth, sep = "_"))
#LocationDepth <- as.character(Fish$LocationDepth)
#FLow_matrix <- matrix(Flow,nrow=length(Flow),ncol=19,byrow=FALSE)



X <- model.matrix(~ Location) 


fam = rep("tweedie", 75)
#fam[25] = "poisson"
#fam[24] = "poisson"
#fam[23] = "poisson"
#fam[22] = "poisson"
#fam[21] = "poisson"
#fam[1] = "lnormal"
#fam[3] = "lnormal"
#fam[4] = "lnormal"
#fam[2] = "lnormal"

# Attach the above together as a list
all = list(Location = Location, Fish_community = Fish_community, fam = fam)

fit.fishnb_no_mycts <- boral(y = Fish_community, family = fam,lv.control = list(num.lv = 2, type = "independent", distmat = NULL), #num.lv = 2, 
                             save.model = TRUE, data = all)
combined_boral <- fit.fishnb_no_mycts
save(combined_boral, file = "combined_boral.Rdata")
load("combined_boral.Rdata")
summ <- summary(combined_boral)
summ
combined_boral$hpdintervals
plot.boral(combined_boral) # to check assumptions
#variance <- calc.varpart(fit.fishnb_no_mycts, groupX = NULL) # should give variance explained

lvsplot(combined_boral, biplot = F, est="median", return.vals = TRUE)
lvs_data <- lvsplot(combined_boral, biplot = T, est="median", ind.spp = 10, return.vals = TRUE)


plot((-(combined_boral$lv.median)), col =c("blue", "red","green3") [Location], 
     xlab = "Latent Variable 1", ylab = "Latent Variable 2")


NMDS_C = data.frame(MDS1 = lvs_data$scaled.lvs[,2], MDS2 = lvs_data$scaled.lvs[,1], group = Location, data = "c) Combined Data")

library(ggplot2)
p1 <- ggplot(data = NMDS_C, aes(MDS1, MDS2))+
  xlab("Latent Variable 1") + ylab("Latent Variable 2") +
  geom_point(aes(color = Location, shape = Location), size = 2) +  theme_bw()+ # scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())#+
  #geom_path(data=df_ell, aes(x=MDS1, y=MDS2, group = group), show.legend = FALSE, size=1, linetype=2)+
  #annotate("text",x=(NMDS.mean$MDS1),y=(NMDS.mean$MDS2),label=NMDS.mean$group)

p1


#### Now do just shape:

## Larval fish analysis
# load data
Fish = read.csv("Otolith_data_mmol_mol_Ca_and_shape.csv", header =T)

str(Fish)
# Convert fish to integers if needed # not sure if needed
#for(n in names(Fish)[1:38]){
#  Fish[[n]]<-as.integer(Fish[[n]])}


# Separate out the columns into important bits
Location = as.factor(Fish[,2])
Fish_community = as.mvabund(Fish[,c(15:77)]) #Select species (<3.4mm ESD)
#Flow = Fish[,3]
#DepthLocation <- as.factor(paste(Location, Depth, sep = "_"))
#LocationDepth <- as.character(Fish$LocationDepth)
#FLow_matrix <- matrix(Flow,nrow=length(Flow),ncol=19,byrow=FALSE)



X <- model.matrix(~ Location) 


fam = rep("tweedie", 63)
#fam[25] = "poisson"
#fam[24] = "poisson"
#fam[23] = "poisson"
#fam[22] = "poisson"
#fam[21] = "poisson"
#fam[1] = "lnormal"
#fam[3] = "lnormal"
#fam[4] = "lnormal"
#fam[2] = "lnormal"

# Attach the above together as a list
all = list(Location = Location, Fish_community = Fish_community, fam = fam)

fit.fishnb_no_mycts <- boral(y = Fish_community, family = fam,lv.control = list(num.lv = 2, type = "independent", distmat = NULL), #num.lv = 2, 
                             save.model = TRUE, data = all)
shape_boral <- fit.fishnb_no_mycts
save(shape_boral, file = "shape_boral.Rdata")

load("shape_boral.Rdata")

summ <- summary(fit.fishnb_no_mycts)
summ
fit.fishnb_no_mycts$hpdintervals
plot.boral(shape_boral) # to check assumptions
#variance <- calc.varpart(fit.fishnb_no_mycts, groupX = NULL) # should give variance explained

lvsplot(shape_boral, biplot = F, est="median", return.vals = TRUE)
lvs_data <- lvsplot(shape_boral, biplot = F, est="median", ind.spp = 10, return.vals = TRUE)


plot((-(shape_boral$lv.median)), col =c("blue", "red","green3") [Location], 
     xlab = "Latent Variable 1", ylab = "Latent Variable 2")


NMDS_S = data.frame(MDS1 = lvs_data$scaled.lvs[,2], MDS2 = lvs_data$scaled.lvs[,1], group = Location, data = "b) Shape Data")

library(ggplot2)
p1 <- ggplot(data = NMDS_S, aes(MDS1, MDS2))+
  xlab("Latent Variable 1") + ylab("Latent Variable 2") +
  geom_point(aes(color = Location, shape = Location), size = 2) +  theme_bw()+ # scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())#+
#geom_path(data=df_ell, aes(x=MDS1, y=MDS2, group = group), show.legend = FALSE, size=1, linetype=2)+
#annotate("text",x=(NMDS.mean$MDS1),y=(NMDS.mean$MDS2),label=NMDS.mean$group)

p1


### Now just do elements:

## COmbined shape and element
# load data
Fish = read.csv("Otolith_data_mmol_mol_Ca_and_shape.csv", header =T)

str(Fish)
# Convert fish to integers if needed # not sure if needed
#for(n in names(Fish)[1:38]){
#  Fish[[n]]<-as.integer(Fish[[n]])}


# Separate out the columns into important bits
Location = as.factor(Fish[,2])
Fish_community = as.mvabund(Fish[,c(3:14)]) #Select species (<3.4mm ESD)
#Flow = Fish[,3]
#DepthLocation <- as.factor(paste(Location, Depth, sep = "_"))
#LocationDepth <- as.character(Fish$LocationDepth)
#FLow_matrix <- matrix(Flow,nrow=length(Flow),ncol=19,byrow=FALSE)



X <- model.matrix(~ Location) 


fam = rep("tweedie", 12)
#fam[25] = "poisson"
#fam[24] = "poisson"
#fam[23] = "poisson"
#fam[22] = "poisson"
#fam[21] = "poisson"
#fam[1] = "lnormal"
#fam[3] = "lnormal"
#fam[4] = "lnormal"
#fam[2] = "lnormal"

# Attach the above together as a list
all = list(Location = Location, Fish_community = Fish_community, fam = fam)

fit.fishnb_no_mycts <- boral(y = Fish_community, family = fam,lv.control = list(num.lv = 2, type = "independent", distmat = NULL), #num.lv = 2, 
                             save.model = TRUE, data = all)
element_boral <- fit.fishnb_no_mycts
save(element_boral, file = "element_boral.Rdata")
load("element_boral.Rdata")
summ <- summary(element_boral)
summ
element_boral$hpdintervals
plot.boral(element_boral) # to check assumptions
#variance <- calc.varpart(fit.fishnb_no_mycts, groupX = NULL) # should give variance explained

lvsplot(element_boral, biplot = F, est="median", return.vals = TRUE)
lvs_data <- lvsplot(element_boral, biplot = T, est="median", ind.spp = 10, return.vals = TRUE)


plot((-(element_boral$lv.median)), col =c("blue", "red","green3") [Location], 
     xlab = "Latent Variable 1", ylab = "Latent Variable 2")


NMDS_E = data.frame(MDS1 = lvs_data$scaled.lvs[,2], MDS2 = lvs_data$scaled.lvs[,1], group = Location, data = "a) Chemistry Data")

library(ggplot2)
p1 <- ggplot(data = NMDS_E, aes(MDS1, MDS2))+
  xlab("Latent Variable 1") + ylab("Latent Variable 2") +
  geom_point(aes(color = Location, shape = Location), size = 2) +  theme_bw()+ # scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())#+
#geom_path(data=df_ell, aes(x=MDS1, y=MDS2, group = group), show.legend = FALSE, size=1, linetype=2)+
#annotate("text",x=(NMDS.mean$MDS1),y=(NMDS.mean$MDS2),label=NMDS.mean$group)

p1


### Combine all data sets
library(dplyr)
combined_plot_data <- rbind_list(NMDS_E, NMDS_S, NMDS_C)

library(plyr)
combined_plot_data$group <- revalue(combined_plot_data$group, c("AA"="Agra", "NN"="Narora", "LL" = "Lucknow"))

str(combined_plot_data)

pALL <- ggplot(data = combined_plot_data, aes(MDS1, MDS2))+
  xlab("Latent Variable 1") + ylab("Latent Variable 2") +
  geom_point(aes(color = group, shape = group), size = 2, alpha = 0.7) +  theme_bw()+ # scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~data) + theme_classic() +
  scale_colour_manual(values = c("black", "red", "blue")) +
  theme(axis.title.x = element_text(face="bold", colour="black", size = 18),
        axis.text.x  = element_text(colour="black", size = 14), 
        axis.title.y = element_text(face="bold", colour="black", size = 18),
        axis.text.y  = element_text(colour="black", size = 14),
        axis.ticks = element_line(colour="black"),
        strip.text = element_text(colour="black", face = "bold", size = 12, hjust=0),
        strip.background = element_rect(colour = "white"),
        #legend.justification=c(1,0), legend.position="right",
        panel.border = element_rect(colour = "black", fill=NA, size = 1),
        #legend.key.size = unit(1, "cm"),
        legend.title=element_blank(),
        legend.text = element_text(size = 12, face = "bold"),
        legend.position = "bottom",
        panel.spacing = unit(1, "lines") # adjusts spacing of panels
  )

pALL

ggsave("Ordination plots.pdf", units = "cm", width = 21, height = 10)
ggsave("Ordination plots.png", units = "cm", width = 21, height = 10, dpi = 600)
