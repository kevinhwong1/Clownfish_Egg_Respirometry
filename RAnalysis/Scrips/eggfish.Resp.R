#Title: Respiration Calculations
#Author: KH Wong
#Modified from: HM Putnam
#Date Last Modified: 20181109
#See Readme file for details

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 
library(devtools)
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') install_github('colin-olito/LoLinR') 
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate') 
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron') 
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("stringr" %in% rownames(installed.packages()) == 'FALSE') install.packages('stringr') 


#Read in required libraries
##### Include Versions of libraries
#install_github('colin-olito/LoLinR')
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library('plyr')
library('dplyr')
library('stringr')
#Required Data files

# Set Working Directory:
setwd("~/MyProjects/Clownfish_Egg_Respirometry/RAnalysis/") #set working

##### Respiration #####
path.p<-"Data/Respirometry/" #the location of all your respirometry files 

# bring in the respiration file names
file.names<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders
#basename above removes the subdirectory name from the file
#add file names that include the subdirectory name
#file.names.full<- list.files(path = path.p, pattern = "csv$", recursive = TRUE) #list all csv file names in the folder and subfolders

#generate a 3 column dataframe with specific column names
Resp.R <- data.frame(matrix(NA, ncol=6))
colnames(Resp.R) <- c("Date", "Run", "Sample.ID","Chamber.ID","Intercept", "umol.L.sec")

Resp.Rb <- data.frame(matrix(NA, ncol=6))
colnames(Resp.Rb) <- c("Date", "Run", "Sample.ID","Chamber.ID","Intercept", "umol.L.sec")


#Load Sample Info
Sample.Info <- read.csv(file="Data/Clown_Sample_Info.csv", header=T) #read sample.info data
#includes treatment, tank, chamber volume, animal size/number etc for normalization
rename <- Sample.Info$Chamber.ID
samp <- Sample.Info$Sample.ID
run <- str_sub(file.names, 10, 13)#file.names
date <- str_sub(file.names, 1, str_length(file.names)-9)#file.names

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Resp.Data <-read.table(file.path(path.p,file.names[i]), skip = 72, header=T, sep=",", na.string="NA", fill = TRUE, as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Resp.Data$Time.Min. <- seq.int(0.017, (nrow(Resp.Data))*0.25, by=0.25) #set time in min
  #Resp.Data[ Resp.Data[,] == "No Sensor" ] <- as.numeric(runif(nrow(Resp.Data), min=0, max=100)) #convert any vials with no data
  Resp.Data.N <- Resp.Data[,3:26]
  
  
  for(j in 1:(ncol(Resp.Data.N))){
    model <- rankLocReg(
      xall=Resp.Data$Time.Min., yall=as.numeric(Resp.Data.N[, j]),
      alpha=0.4, method="pc", verbose=TRUE)
  
  pdf(paste0("~/MyProjects/Clownfish_Egg_Respirometry/RAnalysis/Output/Resp_Plots/",date[i], "_",run[i],"_",rename[j],"_regression.pdf"))
  plot(model)
  dev.off()
  
  Resp.Rb[j,1] <- as.character(date[i]) #stores the date
  Resp.Rb[j,2] <- as.character(run[i]) #stores the run number
  Resp.Rb[j,3] <- as.character(samp[j+(i-1)*ncol(Resp.Data.N)]) #stores the chamber ID
  Resp.Rb[j,4] <- as.character(rename[j]) #stores the chamber ID
  Resp.Rb[j,5:6] <- model$allRegs[i,c(4,5)] #inserts slope and intercept in the dataframe

  
  }
  Resp.R <- rbind(Resp.R, Resp.Rb)
}

Resp.R <- Resp.R[-1,]

write.csv(Resp.R, paste0("~/MyProjects/Clownfish_Egg_Respirometry/RAnalysis/Output/Resp_rates.csv", sep=""))

Resp.R$Merge.ID <- paste0(Resp.R$Date,"_", Resp.R$Run,"_", Resp.R$Sample.ID)
Sample.Info$Merge.ID <- paste0(as.character(Sample.Info$Date), "_", as.character(Sample.Info$Run), "_", as.character(Sample.Info$Sample.ID))

A <- Resp.R$Merge.ID
B <- Sample.Info$Merge.ID

x <- A%in%B

#Merge data with sample info
Data <- merge(Resp.R,Sample.Info, by="Merge.ID" )

Data <- subset(Data, Exclude=!"low.conc")

#Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
Data$umol.sec <- Data$umol.L.sec * Data$Chamber.Vol.L


Data.BK <- subset(Data, Type == "BLANK") #subset to blank data only
plot(Data.BK$umol.sec , ylab="umol O2 sec-1")

#exclude outlier blanks


## BLANKS HAVE TO BE SPECIFIC TO RESPONSE VARIABLE (I.E., PHOTO OR RESP) AND Tank, Treatment, and Timepoint etc...
#Calculate Photo blank rate
blnk <- aggregate(umol.sec ~ Type, data=Data, FUN=mean) #calculate the mean blank rate
Blanks <- subset(blnk, Type == "BLANK") #subset by blanks only
Data$blnk <- Blanks[1,2] #assign blanks to dataframe

plot(blnk$umol.sec ~blnk$Type)

#Account for blank rate Subtract Blank by the temperature blank
Data$umol.sec.corr <- Data$umol.sec-Data$blnk

Data.Sperm <- subset(Data, Type == "SPERM") #subset to coral data only
Data.Egg <- subset(Data, Type == "EGG") #subset to coral data only
Data.Lar <- subset(Data, Type == "LARVAE") #subset to coral data only

#SPERM
#normalize to larval number and min
sperm.conc <- read.csv("Data/CASA_Output.csv")
sperm.conc$Merge.ID <- paste0(as.character(sperm.conc$Date), "_", as.character(sperm.conc$Run), "_", as.character(sperm.conc$Sample.ID))
Data.Sperm$Merge.ID <- paste0(as.character(Data.Sperm$Date.x), "_", as.character(Data.Sperm$Run.x), "_", as.character(Data.Sperm$Coral.ID))

A <- Data.Sperm$Merge.ID 
B <- sperm.conc$Merge.ID 

x <- A%in%B

Data.Sperm <- merge(Data.Sperm,sperm.conc, by="Merge.ID" )
  
Data.Sperm$umol.sperm.sec <- Data.Sperm$umol.sec.corr/Data.Sperm$Conc.sperm.ml
Data.Sperm$umol.sperm.min <- Data.Sperm$umol.sperm.sec*60
Data.Sperm$amol.sperm.min <- Data.Sperm$umol.sperm.min*1000*1000*1000*1000

plot(Data.Sperm$amol.sperm.min ~Data.Sperm$Coral.ID, ylab="amol O2 sperm-1 min-1", las=2)
plot(Data.Sperm$amol.sperm.min ~Data.Sperm$Condition, ylab="amol O2 sperm-1 min-1", las=2)

plot(Data.Sperm$Total.Motility ~Data.Sperm$amol.sperm.min,  las=2)
plot(Data.Sperm$Progressive ~Data.Sperm$amol.sperm.min,  las=2)
plot(Data.Sperm$Slow ~Data.Sperm$amol.sperm.min,  las=2)

pdf("Output/sperm.characteristics.pdf")
par(mfrow=c(2,2))
plot(Data.Sperm$amol.sperm.min ~Data.Sperm$Condition, ylab="amol O2 sperm-1 min-1", las=2)
plot(Data.Sperm$Total.Motility ~Data.Sperm$Condition,  las=2)
plot(Data.Sperm$Progressive ~Data.Sperm$Condition,  las=2)
plot(Data.Sperm$Slow ~Data.Sperm$Condition,  las=2)
dev.off()

#EGGS
#normalize to larval number and min
Data.Egg$umol.egg.sec <- Data.Egg$umol.sec.corr/Data.Egg$Org.Number
Data.Egg$umol.egg.min <- Data.Egg$umol.egg.sec*60
Data.Egg$nmol.egg.min <- Data.Egg$umol.egg.min*1000

plot(Data.Egg$nmol.egg.min ~as.factor(Data.Egg$Org.Number), ylab="nmol O2 egg-1 min-1")


#Larvae
#normalize to larval number and min
Data.Lar$umol.larva.sec <- Data.Lar$umol.sec.corr/Data.Lar$Org.Number
Data.Lar$umol.larva.min <- Data.Lar$umol.larva.sec*60
Data.Lar$nmol.larva.min <- Data.Lar$umol.larva.min*1000

plot(Data.Lar$nmol.larva.min ~Data.Lar$Org.Number, ylab="nmol O2 larva-1 min-1")

