# NMDS plots

library(tidyverse)

elements <- read_csv("Data/Elements NMDS Ordination.csv")
shapes <- read_csv("Data/Shape NMDS Ordination.csv")
combined <- read_csv("Data/Combined NMDS Ordination.csv")

full_dat <- bind_rows(elements, shapes, combined)

library(plyr)
full_dat$Site <- revalue(full_dat$Site, c("A"="Agra", "N"="Narora", "L" = "Lucknow"))

str(full_dat)

pALL <- ggplot(data = full_dat, aes(NMDS_x, NMDS_y))+
  xlab("nMDS 1") + ylab("nMDS 2") +
  geom_point(aes(color = Site, shape = Site), size = 2, alpha = 0.7) +  theme_bw()+ # scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_wrap(~Data, scales = "free") + theme_classic() +  
  #scale_x_continuous(labels = c("-1", "-0.5", "0", "0.5", "1"), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  #scale_y_continuous(labels = c("-1", "-0.5", "0", "0.5", "1"), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  scale_colour_manual(values = c("black", "red", "blue")) +
  theme(axis.title.x = element_text(face="bold", colour="black", size = 16),
        axis.text.x  = element_text(colour="black", size = 12), 
        axis.title.y = element_text(face="bold", colour="black", size = 16),
        axis.text.y  = element_text(colour="black", size = 12),
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

ggsave("plots/NMDS Ordination plots.pdf", units = "cm", width = 21, height = 10)
ggsave("plots/NMDS Ordination plots.png", units = "cm", width = 21, height = 10, dpi = 600)
