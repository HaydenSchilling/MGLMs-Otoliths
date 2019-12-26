# Distance based test to see if this data meets the assumptions

library(vegan)

# load data
Fish = read.csv("../Data/Otolith_data_mmol_mol_Ca_and_shape.csv", header =T)

str(Fish)

# Separate out the columns into important bits
Location = as.factor(Fish[,2])
Fish_elements = dist(Fish[,c(3:14)]) #Select elements only

# Run the test
fitX <- betadisper(Fish_elements, Location)
summary(fitX)
fitX
anova(fitX) ### There are significant differences in variance between these groups for element data
plot(fitX) ## Visualise differences (AA is very different)

######################################################################################

### Now do shape data

# Separate out the columns into important bits
Location = as.factor(Fish[,2])
Fish_shape = dist(Fish[,c(15:77)]) #Select shape data only

# Run the test
fitX <- betadisper(Fish_shape, Location)
summary(fitX)
fitX
anova(fitX) ### NO significant differences in variance between these groups for shape data
plot(fitX) ## Visualise differences (Lots of overlap)

###########################################################################################

### Now do combined element and shape data

# Separate out the columns into important bits
Location = as.factor(Fish[,2])
Fish_community = dist(Fish[,c(3:77)]) #Select elements and shape
#Flow = Fish[,3]

# Run the test
fitX <- betadisper(Fish_community, Location)
summary(fitX)
fitX
anova(fitX) ### significant differences in variance between these groups for shape data
plot(fitX) ## Visualise differences (AA is very different), similar to just elements


### NMDS visual (using element data)

Fish_community = Fish[,c(3:14)] # change this to do shape or elements and shape

example_NMDS=metaMDS(Fish_community,k=2,trymax=100, distance = "euclidean")
plot(example_NMDS)



MDS_xy <- data.frame(example_NMDS$points)
MDS_xy$Location <- Location

library(ggplot2)
ggplot(MDS_xy, aes(MDS1, MDS2, color = Location)) + geom_point() + theme_bw()
# gave similar pattern but differences in dispersal are very evident
