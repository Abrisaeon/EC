#EDDY COVARIANCE POST-PROCESSING USING REDDYPROC:
#Script modified to include extra quality checks and summary stats related to gaps
#Includes: 
#1) What percentage of data are NA due to rainfall, instrument downtime, etc.
#2) Despiking according to plausibility thresholds (Brummer, pers.comm. 2019) - followed by summary stats on gaps
#3) Selection of quality criteria - Fundamental research grade or FLUXNET grade (Foken 2012) - followed by summary stats on gaps
#4) Preparation of a ReddyProc input file (calculations, dates & headers)
#5) Standard ReddyProc workflow (after summary stats on gaps)


#install.packages("REddyProc")

library(REddyProc)
library(dplyr)
#explore examples 
vignette('useCase')

#------------------------------------------------------------------------------------------------------------------------------------------------------------
#Prior to REddyproc postprocessing, data from Easyflux needs to be subset, plausibility checks run and quality filtered
setwd("C:/Users/Abri.SAEON/Documents/SAEON/Emergency data repository/Jonkershoek EC/rtest")

dat<-read.csv("CR3000 Jonkershoek EC_FluxCopy.csv",skip=0)

#Subset the columns that REddyProc needs
d<-subset(dat, select = c(TIMESTAMP, Fc_molar, Fc_qc_grade, LE, LE_qc_grade, H, H_qc_grade, Rn, u_star, Tc_Avg, RH_Avg, VPD_air,Tsoil_Avg1, Tsoil_Avg2))

#Evaluate time series completeness (% NAs due to rainfall, instrument downtime, etc.)
#%NAs in time series
(sum(is.na(d$Fc_molar))/nrow(d))*100 #CO2 flux
(sum(is.na(d$LE))/nrow(d))*100       #Evapotranspiration


#Plausibility thresholds (Brummer, 2019 - pers.comm.)
d$Fc_molar<- ifelse(d$Fc_molar < -50 | d$Fc_molar > 50, NA, d$Fc_molar)
d$LE<- ifelse(d$LE < -200 | d$LE > 1000, NA, d$LE)
d$H<- ifelse(d$H < -200 | d$H > 1000, NA, d$H)

#Summary stats for time series completeness (after implausible spikes are thrown out)
#%NAs in time series after plausibility filter
(sum(is.na(d$Fc_molar))/nrow(d))*100 #CO2 flux
(sum(is.na(d$LE))/nrow(d))*100       #Evapotranspiration


#DECIDE ON DATA QUALITY APPROPRIATE FOR YOUR PURPOSE
#Exclude data with Foken (2012) quality grading > 6 (general use, ie. for FLUXNET)
d$Fc_molar<- ifelse(d$Fc_qc_grade > "6", NA, d$Fc_molar )
d$LE<- ifelse(d$LE_qc_grade > "6", NA, d$LE )
d$H<- ifelse(d$H_qc_grade > "6", NA, d$H )

#Summary stats for time series completeness after a < grading 6 filter
#%NAs in time series after plausibility filter
(sum(is.na(d$Fc_molar))/nrow(d))*100 #CO2 flux
(sum(is.na(d$LE))/nrow(d))*100       #Evapotranspiration
#gaps of around 30% after this step are apparently not uncommon - ours is high

#Exclude data with Foken (2012) quality grading > 3 (fundamental research grade)
#d$Fc_molar<- ifelse(d$Fc_qc_grade > "3", NA, d$Fc_molar )
#d$LE<- ifelse(d$LE_qc_grade > "3", NA, d$LE )
#d$H<- ifelse(d$H_qc_grade > "3", NA, d$H )

#Summary stats for time series completeness after a < grading 3 filter
#%NAs in time series after plausibility filter
(sum(is.na(d$Fc_molar))/nrow(d))*100 #CO2 flux
(sum(is.na(d$LE))/nrow(d))*100       #Evapotranspiration

#
d$Rn<- ifelse(d$Rn < 0, 0, d$Rn )
d$Tsoil<- rowMeans(d[c('Tsoil_Avg1', 'Tsoil_Avg2')], na.rm=TRUE)
x<-subset(d, select = c(TIMESTAMP, Fc_molar, LE, H, Rn, u_star, Tc_Avg, RH_Avg, VPD_air, Tsoil))

#Prepare date format for REddyProc
x$datetime <- strftime(x$TIMESTAMP, format = "%Y-%m-%d %H:%M")  
x$year <- strftime(x$datetime, format = "%Y")  
x$DT<- as.POSIXlt(x$datetime)
x$doy<- strftime(x$datetime, format = "%j")
x$doy<-as.numeric(x$doy)
x$hour<-strftime(x$datetime, format = "%H")
x$min<-strftime(x$datetime, format = "%M")
x$HM<- as.numeric(x$hour) + (as.numeric(x$min)/60) 
z<-subset(x, select = c(year, doy, HM, Fc_molar, LE, H, Rn, Tc_Avg , Tsoil, RH_Avg, VPD_air,u_star))
colnames(z)<- c("Year", "DoY", "Hour", "NEE", "LE", "H", "Rg", "Tair" , "Tsoil", "rH", "VPD", "Ustar")

write.csv(z,"JHK.csv", sep="\t", col.names=T,quote=F, row.names=T, na = "NA")

#now ensure the saved file starts with the first half hour of data of the first day
#also ensure the saved file ends with the last half hour of data for the last day


dat<-read.csv("JHK.csv",skip=0)

#Summary stats about gaps
#%NAs in NEE timeseries
(sum(is.na(dat$NEE))/nrow(dat))*100
(sum(is.na(dat$LE))/nrow(dat))*100


dat$VPD<- ifelse(dat$VPD < 0, NA, dat$VPD )

#Replace long runs of equal NEE values by NA
dat<- filterLongRuns(dat, "NEE")

#Add time stamp in POSIX time format
EddyDataWithPosix <- fConvertTimeToPosix(
  dat, 'YDH', Year = 'Year',
  Day = 'DoY',Hour = 'Hour')%>% 
  filterLongRuns("NEE")

eddyC <- sEddyProc$new(
  'JHK', EddyDataWithPosix,
  c('NEE','Rg','Tair','VPD', 'Ustar'))

#Fingerprint plot of raw NEE data
eddyC$sPlotFingerprintY('NEE', Year = 2019)

#DETERMINING USTAR THRESHOLDS (USING REDDYPROC)
#--------------------------------------------------------
#Determine Ustar thresholds
eddyC$sEstimateUstarScenarios(
  nSample = 100L, probs = c(0.05, 0.5, 0.95))
eddyC$sGetEstimatedUstarThresholdDistribution()

# inspect the thresholds to be used by default (if you used REddyProc to determine thresholds)
eddyC$sGetUstarScenarios()

((sum(dat$Ustar > 0.13))/nrow(dat))*100 #what percentage of the data have Ustar>0.19
((sum(dat$Ustar > 0.19))/nrow(dat))*100 #what percentage of the data have Ustar>0.19
((sum(dat$Ustar > 0.44))/nrow(dat))*100 #what percentage of the data have Ustar>0.44

#uStarTh <- eddyC$sEstUstarThresholdDistribution(
#  nSample = 100L, probs = c(0.05, 0.5, 0.95))
#uStarTh %>%
#  filter( aggregationMode == "year") %>%
#  select( uStar, "5%", "50%", "95%")

#uStarThAgg <- eddyC$sGetEstimatedUstarThresholdDistribution()

#uStarSuffixes <- colnames(eddyC$sGetUstarScenarios())[-1]

#uStarThAnnual <-
#  usGetAnnualSeasonUStarMap(uStarTh)[-2]
#uStarSuffixes <- colnames(uStarThAnnual)[-1]
#print(uStarThAnnual)

#SETTING USTAR THRESHOLDS MANUALLY 
#--------------------------------------------------------
#0.1 to 0.13 are options. Test various choices and see what gets rid of rubbish data (Brummer, 2019 pers.com)
#If you want to set uStar thresholds manually, use the following:
#uStar <- 0.13
#eddyC$sMDSGapFillAfterUstar('NEE', uStarTh = uStar)


#GAP FILLING
#--------------------
# Fill a dataset that has been filtered With manually set uStar threshold
#eddyC$sMDSGapFillAfterUstar('NEE', uStarTh = uStar)

#grep("NEE_uStar_f$",names(eddyC$sExportResults()), value = TRUE)
#eddyC$sPlotFingerprintY('NEE_uStar_f', Year = 2019)

#-----------------------------------------------------------------
#GAP FILLING
# Fill a dataset that has been filtered With REddyProc determined uStar thresholds
eddyC$sMDSGapFillUStarScens('NEE')

#sEddyProc_sMDSGapFillAfterUstar('NEE', 
                                #uStarVar = "Ustar", uStarTh = uStarThAnnual[, 
                                #c("season", uStarSuffixes), drop = FALSE], 
                                #uStarSuffix = "uStar", isFlagEntryAfterLowTurbulence = FALSE, 
                                #isFilterDayTime = TRUE, swThr = 0, 
                                #RgColName = "Rg")

#eddyC$sMDSGapFillAfterUStarDistr('NEE',
                                      #UstarThres.df = uStarThAnnual,
                                      #UstarSuffix.V.s = uStarSuffixes,
                                      #FillAll = TRUE
#)


grep("NEE_U05_f$",names(eddyC$sExportResults()), value = TRUE) #Filled
grep("NEE_U05_fsd$",names(eddyC$sExportResults()), value = TRUE) #Estimated standard deviation of uncertainty

grep("NEE_U50_f$",names(eddyC$sExportResults()), value = TRUE)
grep("NEE_U50_fsd$",names(eddyC$sExportResults()), value = TRUE)

grep("NEE_U95_f$",names(eddyC$sExportResults()), value = TRUE)
grep("NEE_U95_fsd$",names(eddyC$sExportResults()), value = TRUE)

eddyC$sPlotFingerprintY('NEE_U50_f', Year = 2019)
eddyC$sPlotFingerprintY('NEE_U05_f', Year = 2019)
eddyC$sPlotFingerprintY('NEE_U95_f', Year = 2019)

#FLUX PARTITIONING
#-----------------------------------------------------------------------------------------------------------------------
#Enter site coordinates and time zone to use in calculation of sunrise & sunset
eddyC$sSetLocationInfo(LatDeg = -33.9, LongDeg = 18.9, TimeZoneHour = +2)  
eddyC$sMDSGapFill('Tair', FillAll = FALSE, minNWarnRunLength = NA)      
eddyC$sMDSGapFill('VPD', FillAll = FALSE, minNWarnRunLength = NA)


eddyC$sMRFluxPartitionUStarScens()

#variable uStarSuffixes was defined above at the end of uStar threshold estimation
#resPart <- lapply(uStarSuffixes, function(suffix){
#  eddyC$sMRFluxPartition(Suffix.s = suffix)
#})


grep("GPP_U05_f$|Reco",names(eddyC$sExportResults()), value = TRUE)
grep("GPP_U50_f$|Reco",names(eddyC$sExportResults()), value = TRUE)
grep("GPP_U95_f$|Reco",names(eddyC$sExportResults()), value = TRUE)
grep("GPP_uStar_f$|Reco",names(eddyC$sExportResults()), value = TRUE)

eddyC$sPlotFingerprintY('GPP_uStar_f', Year = 2019)



eddyC$sPlotFingerprintY('GPP_U05_f', Year = 2019)
eddyC$sPlotFingerprintY('GPP_U50_f', Year = 2019)
eddyC$sPlotFingerprintY('GPP_U95_f', Year = 2019)
eddyC$sPlotFingerprintY('Reco_U05', Year = 2019)
eddyC$sPlotFingerprintY('Reco_U05', Year = 2019)
eddyC$sPlotFingerprintY('Reco_U50', Year = 2019)
eddyC$sPlotFingerprintY('Reco_U95', Year = 2019)

FilledEddyData <- eddyC$sExportResults()
uStarSuffixes <- colnames(eddyC$sGetUstarScenarios())[-1]


#Calculate mean GPP across all years for each u* scenario
GPPAggCO2 <- sapply( uStarSuffixes, function(suffix) {
  GPPHalfHour <- FilledEddyData[[paste0("GPP_",suffix,"_f")]]
  mean(GPPHalfHour, na.rm = TRUE)
})

#Convert GPP from ??molCO2m???2s???1 to gCm???2yr???1.
molarMass <- 12.011
GPPAgg <- GPPAggCO2 * 1e-6 * molarMass * 3600*24*365.25

#??molCO2m???2s???1 to gCm???2yr???1.
print(GPPAgg) #The difference between these aggregated values is a first estimate of uncertainty range in GPP due to uncertainty of the u??? threshold.

(max(GPPAgg) - min(GPPAgg)) / median(GPPAgg) 

#Export post processed results to a .csv
FilledEddyData <- eddyC$sExportResults()
CombinedData <- cbind(EddyDataWithPosix, FilledEddyData)
fWriteDataframeToFile(CombinedData, 'JNK-Results.txt', Dir.s = "C:/Users/Abri.SAEON/Documents/SAEON/Emergency data repository/Jonkershoek EC/rtest")

########################

d<-read.csv("JNK-Results.csv",skip=0)

plot(NEE_U05_f~ Date.Time, data = d)
ggplot(d, aes(x=Date.Time, y=NEE_U05_f, group = 1)) + geom_line(colour = "red")   

ggplot(d, aes(x=Date.Time, y=Reco_U95)) + geom_point() 

ggplot(d, aes(x=Date.Time, y=GPP_U95_f)) + geom_point() 

#################################
#Openeddy
#-----------
dat<-read.csv("JHK.csv",skip=0)
EddyData.F <- read_eddy("JHK.csv")

head(EddyData.F)
str(EddyData.F)


#dat$VPD<- ifelse(dat$VPD < "0", "NA", dat$VPD )
#dat$VPD<- as.numeric(dat$VPD)

EddyDataWithPosix.F <- fConvertTimeToPosix(
  EddyData.F, 'YDH', Year.s = 'Year',
  Day.s = 'DoY',Hour.s = 'Hour')


variables <- c('NEE', 'LE', 'H', 'Rg','Tair', 'Tsoil', 'VPD', 'Ustar')
EddyProc.C <- sEddyProc$new(siteyear, EddyDataWithPosix.F, variables)
EddyProc.C$sSetLocationInfo(lat, long, tz)  # site location info

# See the content
str(EddyProc.C)
EddyProc.C$sPrintFrames(NumRows.i = 6L)

# Plot individual months/years to console (of current R graphics device)
if (plot_to_console) {
  for (Var in variables) {
    EddyProc.C$sPlotFingerprintY(Var, Year.i = year)
    title(main = Var, line = 0.2)
  }
  for (Var in variables) {
    EddyProc.C$sPlotHHFluxesY(Var, Year.i = year)
    title(ylab = Var, line = 2)
  }
}



eddyC <- sEddyProc$new(
  'JHK', EddyDataWithPosix.F,
  c('NEE','Rg','Tair','VPD', 'Ustar'))

colnames(dat)[colnames(dat)=="DateTime"] <- "timestamp"
colnames(dat)[colnames(dat)=="Rg"] <- "GR"
#***********************************************

dat$ne <- ifelse((dat$wnd_dir_compass > 22.51)&(dat$wnd_dir_compass < 67.5),1,0)
dat$e <- ifelse((dat$wnd_dir_compass > 67.51)&(dat$wnd_dir_compass < 112.5),1,0)
dat$se <- ifelse((dat$wnd_dir_compass > 112.5)&(dat$wnd_dir_compass < 157.5),1,0)
dat$s <- ifelse((dat$wnd_dir_compass > 157.51)&(dat$wnd_dir_compass < 202.5),1,0)
dat$sw <- ifelse((dat$wnd_dir_compass > 202.51)&(dat$wnd_dir_compass < 247.5),1,0)
dat$w <- ifelse((dat$wnd_dir_compass > 247.51)&(dat$wnd_dir_compass < 292.5),1,0)
dat$nw <- ifelse((dat$wnd_dir_compass > 292.51)&(dat$wnd_dir_compass < 337.5),1,0)
dat$dummy1<- ifelse((dat$wnd_dir_compass > 337.51)&(dat$wnd_dir_compass < 360),1,0)
dat$dummy <- ifelse((dat$wnd_dir_compass > 0)&(dat$wnd_dir_compass < 22.5),1,0)
dat$dummy1<-as.numeric(dat$dummy1)
dat$dummy<-as.numeric(dat$dummy)
dat$n<-dat$dummy1+dat$dummy

secount <- sum(na.omit(dat$se))		#Count days with rain as identified above
scount <- sum(na.omit(dat$s))
swcount <- sum(na.omit(dat$sw))
wcount <- sum(na.omit(dat$w))
nwcount <- sum(na.omit(dat$nw))
ncount <- sum(na.omit(dat$n))
necount <- sum(na.omit(dat$ne))
ecount <- sum(na.omit(dat$e))

#Bind all the count data together and select relevant columns for the rest of the analysis
wind<-cbind(ncount, necount, ecount, secount, scount, swcount, wcount, nwcount )

#Put the data in the long form for plotting and set column names
dw <- melt(wind ,  id.vars = 'count', variable.name = 'dir')					#Melt the dataframe into the long form in order to plot it
dw <- melt(wind)
dw<-subset(dw, select = c(X2, value))


#Create an object that contains a graph
ggplot(dw, aes(x=X2, y=value)) + geom_bar(aes(fill=X2), stat="identity", position="stack") + scale_fill_grey(start = .5, end = .5) + coord_polar(theta="x", direction = 10, start = 5)  + theme_bw() + ylab("Half hour count") + xlab("") + scale_y_continuous(limits = c(0, 2000))+ ggtitle("JEC - winddir") + theme(legend.position = "none") + theme(axis.text=element_text(size=13)) + theme(axis.text.x = element_text(angle = 0, hjust = 1))


dat$date <- as.POSIXct(dat$TIMESTAMP, tz = "Africa/Johannesburg")
dat$dates <- strftime(dat$date, format = "%Y-%m-%d %H:%M")  

summaryPlot(select(dat, date:wnd_dir_compass))
summaryPlot(select(dat, date:Fc_molar))

jecwind<- subset(dat, select = c(dates, wnd_spd, wnd_dir_compass ))
colnames(jecwind)<- c("datetime", "wind_speed", "wind_dir")

windRose(jecwind, key.footer = "knots")
