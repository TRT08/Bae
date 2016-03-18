library(plyr)
library(ggplot2)
library(lubridate)

setwd("C:/Users/trtucker/Desktop/weather")

#files <- list.files()
#bw <- do.call("rbind.fill",lapply(files,FUN=function(files){read.table(files,header=TRUE,fill = TRUE)}))
#bw$YY[bw$YY <1999] <- bw$YY[bw$YY <1999]+1900
#bw$date <- as.Date(with(bw, paste(YY, MM, DD, sep="-")), "%Y-%m-%d")
#write.csv(bw, "Buoy_all.csv")


bw <- read.csv("Buoy_all.csv", header=TRUE)
bw <- subset(bw, select=-c(X.1,X))

bw$DayNum <- yday(bw$date)

bw <- bw[!(bw$YY == 2010 & bw$DayNum < 100),] #remove messed up 2010 data

bw[bw == 9999] <- NA 
bw[bw == 999] <- NA 
bw[bw == 99] <- NA 

bw$YYF <- as.factor(bw$YY)
ggplot(data = bw, aes(x = DayNum, y = WTMP, colour = YYF)) + geom_line()

plot(bw$date, bw$WTMP, type='p', ylab='Mortality', main='lowess')
lines(lowess(bw$date ~ bw$WTMP, f=0.02), col='blue')# seasonal component (2% of data)
lines(lowess(bw$WTMP ~ bw$date, f=2/3))# trend (2/3 of data)
