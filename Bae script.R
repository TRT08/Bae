
####################################
####GETTING STARTED - RUN FIRST ####
####################################
require(reshape) || install.packages("reshape") 
require(plyr) || install.packages("plyr") 
require(RODBC) || install.packages("RODBC") 
require(splitstackshape) || install.packages("splitstackshape") 
require(data.table) || install.packages("data.table") 
require(car) || install.packages("car") 
require(doBy) || install.packages("doBy")
require(chron) || install.packages("chron")
require(SDMTools) || install.packages("SDMTools")

source("F:/DATA/SLBE/Manuscripts/Benthos-analysis/biomassfreqscript.R")
source("F:/DATA/SLBE/R scripts/General-Scripts/GobyThemes.R") 

options(scipen=999) #Keeps scientific notation away
names(Freq.All.benthos)[names(Freq.All.benthos)=="Sp.YearSER"] <- "YearSER"
setwd("F:/DATA/SLBE/Manuscripts/Bae/")
####################################

interestingSERS <- c("2013-154", "2013-305", "2013-462")

CompiledBenthosLog <- CompiledBenthosLog[CompiledBenthosLog$YearSER %in% interestingSERS,]
ben.all.L <- ben.all.L [ben.all.L $YearSER %in% interestingSERS,]
CompiledBenthosLogandFreq <- CompiledBenthosLogandFreq[CompiledBenthosLogandFreq$YearSER %in% interestingSERS,]
Freq.All.benthos <- Freq.All.benthos[Freq.All.benthos$YearSER %in% interestingSERS,]
NA.omit.biomass <- NA.omit.biomass[NA.omit.biomass$YearSER %in% interestingSERS,]

############MAKE PLOTS OF THE BT DATA##############
####### Get BT data from Database############
db <- "F:/DATA/SLBE/AVBOT Database.accdb"
con2 <- odbcConnectAccess2007(db)
BT <- sqlFetch(con2, "BT Environmental Data")
close(con2)
colnames(BT) <- make.names(colnames(BT), unique = TRUE)
DownCasts <- BT[which(BT$Cast.Direction=="Down"),]
names(DownCasts)[names(DownCasts)=="YearSer"] <- "BT.YearSer"

BT2 <- join(DownCasts, CompiledBenthosLog, by="BT.YearSer")
BT2 <- BT2[BT2$YearSER %in% interestingSERS,]

setwd("F:/DATA/SLBE/Manuscripts/Bae/Figs/")
for (i in names(BT2)[c(10,12:16, 19:20)]){
  ggplot(BT2, aes(x=BT2[[i]], y=Depth..m., group= as.factor(Date), colour=as.factor(Date))) + geom_point()+ geom_path(aes(colour=as.factor(Date))) + scale_y_reverse() + ggtitle(paste(i))+ xlab(paste(i))  + Goby_theme
  ggsave(filename = paste(i, ".png"), plot = last_plot() ,height = 4, width = 6)
}        

####################IMAGEJ DATA############################

setwd("F:/DATA/SLBE/Manuscripts/Bae/ImageJ/Rock analysis/Data/")
files <- list.files()
files <- as.factor(files)

filelist<- list()
for (i in levels(files)){
  tab <- read.table(i, header = TRUE, sep = "\t")
  tab$File <- paste(i)
  filelist[[i]]<- tab
}

rocks <- rbind.fill(filelist)
write.csv(rocks, "All_rocks.csv")
rm(files, i, filelist, tab)

####ANALYZE rock data######

library(doBy)
mean(rocks$Area)
sd(rocks$Area)
min(rocks$Area)
max(rocks$Area)
nrow(rocks)


###########MAPS THE BAES###############
library(maps)
library(mapdata)
library(ggplot2)
library(rgdal) 
library(maptools)

abun <- read.csv("F:/DATA/SLBE/Manuscripts/Bae/Doc/Abundances GL.csv")
#remove the ones we don't want to map

setwd("F:/DATA/SLBE/Manuscripts/Bae/Maps/") 

locs1 = readOGR(dsn=".", layer="BaeLocs")
locs1  <- spTransform(locs1 , CRS("+proj=longlat"))
locs1@data$Id = rownames(locs1@data)
DB <- data.frame(FID=locs1$Id)
DB <- cbind(DB, locs1@coords)
DB$NP <- rep("NP",nrow(DB))



#locs1 <- fortify(locs1)
ggplot(DB) + aes(coords.x1,coords.x2) + geom_point()


abun <- abun[abun$Include.on.map. == "Y", ]
abun$sqm.adj <- log10(abun$sqm.adj)
abun$sqm.adj[is.na(abun$sqm.adj)] <- 0

levels(abun$Present.on.map)[levels(abun$Present.on.map)=="N"] <- "Count"
levels(abun$Present.on.map)[levels(abun$Present.on.map)=="P"] <- "Present"


# get the names
MRmap <- fortify(map('worldHires', c('Canada'), fill=TRUE, plot=FALSE))
us = map_data("state")

#sub2=expression(atop(paste("Mean ", italic("Clad.")),paste("score")))

setwd("F:/DATA/SLBE/Manuscripts/Bae/Figs/")

MYcolors <-c("#999999", "#000000")

gg <- ggplot()
gg <- gg + geom_point(data=DB, aes(x=coords.x1, y=coords.x2, color=NP), color="light grey", shape=1) 
gg <- gg + geom_map(data=MRmap , map=MRmap, aes(x=long, y=lat, map_id=region), fill="white", color="black") + 
  geom_polygon(data=us,aes(x=long, y=lat, group=group), colour="black", fill="white", alpha=0)  +        
  coord_map(xlim = c(-93,-76),ylim = c(41, 49)) + Goby_theme
gg <- gg + geom_point(data=abun, aes(x=Long, y=Lat, shape = Present.on.map, size = sqm.adj, color=Age))  + 
  scale_colour_manual(values=MYcolors) +
  scale_size_continuous(range = c(2, 7)) + 
  xlab("Longitude") +  ylab("Latitude") + Goby_theme +
  labs(size = "log(Count)", shape = "Shape Key", color="Life stage") +
  theme(legend.position = c(.9, .6), legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="black"))+
  guides(color=guide_legend(), size = guide_legend())
gg

png("BaeLocs.png", width = 6, height = 7, units = 'in', pointsize = 12, res = 400)
gg
dev.off()


#p <- p + geom_point(data=df, aes(x=Long.1, y=Lat.1, size = meanClad, color=CV)) + 
#  scale_size_continuous(range = c(3, 7), name=sub2, guide=guide_legend(override.aes=aes(shape=1))) +
#  scale_color_gradient(low="gray87", high="gray0" ) + xlab("Longitude") +  ylab("Latitude") + Goby_theme 



