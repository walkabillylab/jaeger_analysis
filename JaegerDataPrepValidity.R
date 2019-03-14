# Created by Javad Rahimipour Anaraki on 26/02/18 updated 01/02/19
# Ph.D. Candidate
# Department of Computer Science
# Memorial University of Newfoundland
# jra066 [AT] mun [DOT] ca | www.cs.mun.ca/~jra066

#   input: Data from Jaeger Oxycon Pro
#  output: Interpolated data

rm(list=ls())
#========================Libraries=========================
list.of.packages <-
  c("stringr",
    "data.table",
    "kimisc",
    "imputeTS",
    "dplyr")
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages)

library(stringr)
library(data.table)
library(kimisc)
library(imputeTS)
library(dplyr)

#=========================Variables========================
OS <- Sys.info()
if (OS["sysname"] == "Windows") {
  path <-
    "Z:/Research/dfuller/Walkabilly/studies/smarphone_accel/data/Jaeger/"
  intrPath <-
    "Z:/Research/dfuller/Walkabilly/studies/smarphone_accel/data/"
} else {
  path <-
    "/Volumes/hkr-storage/Research/dfuller/Walkabilly/studies/smarphone_accel/data/Jaeger/"
  intrPath <-
    "/Volumes/hkr-storage/Research/dfuller/Walkabilly/studies/smarphone_accel/data/"
}
setwd(path)

#Required user id to be processed
uid <- "128"

#Imputation method
method <- "spline"

#Timezone
timeZone <- "America/St_Johns"
Sys.setenv(TZ = timeZone)


#===================Data files location====================
filenames <-
  list.files(paste(path, uid, sep = ""),
             pattern = "[0-9].csv",
             full.names = TRUE)

intervals <-
  fread(paste(intrPath, "intervals.csv", sep = ""),
        sep = ",",
        data.table = FALSE)

timing <-
  fread(paste0(path, uid, "/timing.csv"),
        sep = ",",
        data.table = FALSE)


#================Loop over wear locations==================
for (i in 1:length(filenames)) {
  #Reading data
  inData <-
    read.table(
      file = filenames[i],
      header = FALSE,
      skip = 11,
      skipNul = TRUE,
      blank.lines.skip = TRUE,
      sep = ",",
      col.names = c(
        "Time",
        "PHR",
        "VO2",
        "VO2kg",
        "VCO2",
        "VE",
        "BF",
        "RER",
        "EqO2",
        "Speed",
        "Elev.",
        "MET"
      )
    )
  
  #Filtering data out based on uid and start and end date
  usrInfo <- intervals[which(intervals$userid == uid),]
  startDate <- usrInfo[, "start"]
  endDate <- usrInfo[, "end"]
  
  #Fix Jaeger time
  colnames(inData)[1] <- "record_time"
  inData[, "record_time"] <- trimws(as.character(inData[, "record_time"]))
  exctOne <- union(which(inData$record_time == "1:00:00"), which(inData$record_time == "01:00:00"))
  likeOne <- which(inData$record_time %like% "01:00:0")
  
  if(length(exctOne) >= 2) {
    exctOne <- exctOne[length(exctOne)] - 1
    inData[1:exctOne, "record_time"] <- paste0("00:", substring(inData[1:exctOne, "record_time"], 1, 5))
    
  } else if (length(likeOne) > 0) {
    likeOne <- likeOne[1] - 1
    inData[1:likeOne, "record_time"] <- paste0("00:", substring(inData[1:likeOne, "record_time"], 1, 5))
  }
  
  #Convert Jaeger Time to desired record_time
  for (j in 1:nrow(inData)) {
    pos <-
      unlist(strsplit(inData[j, "record_time"], "[:]"))
    if (length(pos) > 2) {
      inData[j, "record_time"] <-
        as.character(
          as.POSIXlt(startDate, TZ = timeZone) + as.numeric(pos[3]) + as.numeric(pos[2]) *
            60 + as.numeric(pos[1]) * 60 * 60
        )
    } else {
      inData[j, "record_time"] <-
        as.character(as.POSIXlt(startDate, TZ = timeZone) + as.numeric(pos[2]) + as.numeric(pos[1]) *
                       60)
    }
  }
  
  #Fill in intervals in second level
  cData <-
    data.frame(seq.POSIXt(
      as.POSIXlt(startDate, TZ = timeZone),
      as.POSIXlt(endDate, TZ = timeZone),
      by = "sec"
    ))
  
  #Merging in data with second level data
  colnames(cData) <- "record_time"
  cData$record_time <- as.character(cData$record_time)
  lData <- merge(inData, cData, by = "record_time", all.y = TRUE)
  
  #Convert zeros and "-" to NAs and remove all NA columns | convert columns to numeric
  for (j in 2:ncol(lData)) {
    lData[which(lData[, j] == 0), j] <- NA
    lData[which(lData[, j] == "-"), j] <- NA
    lData[, j] <- as.numeric(as.character(lData[, j]))
  }
  
  #Remove columns with no values / save data for calculating MSE
  lData[, which(colSums(is.na(lData)) == nrow(lData))] <- NULL
  
  #Apply bandpass filter to RER column
  lData[which(lData$RER < 0.7), "RER"] <- 0.7
  lData[which(lData$RER > 1.2), "RER"] <- 1.2
  
  #Column for Energy Expenditure
  lData$EE <- 0
  
  #Set labels
  lData$activity <- "-"
  for (l in 1:nrow(timing)) {
    rows <-
      which((
        format(lData$record_time, format = '%H:%M:%S') >= as.POSIXlt(timing[l, "Start time"])
      ) &
        (
          format(lData$record_time, format = '%H:%M:%S') <= as.POSIXlt(timing[l, "End time"])
        ))
    lData[rows, "activity"] <- timing[l, "Task name"]
    
    #Impute data
    lData[rows, ] <-
      na.interpolation(lData[rows, ], option = method)
  }
  
  #Set label for intervals without label and apply imputation
  rows <- which(lData[, "activity"] == "-")
  lData[rows, "activity"] <- "transit"
  lData[rows,] <- na.interpolation(lData[rows,], option = method)
  
  #Keep original data for plotting
  plotData <- lData
  
  #Generating data using Sline curve fitting model
  lData$VO2 <- smooth.spline(1:nrow(lData), lData$VO2, df = 50)$y
  lData$VCO2 <- smooth.spline(1:nrow(lData), lData$VCO2, df = 50)$y
  lData$RER <- lData$VCO2 / lData$VO2
  
  #Calculate Energy Expenditure (Kcal/min)
  lData$EE <- -(0.550 * (lData$VCO2 / 1000) - 4.47 * (lData$VO2 / 1000))

  #Plot and save
  jpeg(paste0(path, uid, '/',unlist(strsplit(basename(filenames), "[.]"))[1], "_VO2.jpeg"))
  org <- cbind(plotData$VO2, lData$VO2)
  plot(
    org[,1],
    col = "green",
    ylab = "Original VO2",
    xlab = "Time",
    type = "l"
  )
  par(new = TRUE)
  lines(
    org[,2],
    xlab = NA,
    ylab = NA,
    col = "red"
  )
  axis(side = 4)
  mtext('Smoothed VO2', side = 4, line = 2)
  dev.off()
  jpeg(paste0(path, uid, '/',unlist(strsplit(basename(filenames), "[.]"))[1], "_VO2_BlandAltman.jpeg"))
  plot(org[,2] - org[,1])
  dev.off()
  
  jpeg(paste0(path, uid, '/',unlist(strsplit(basename(filenames), "[.]"))[1], "_VCO2.jpeg"))
  org <- cbind(plotData$VCO2, lData$VCO2)
  plot(
    org[,1],
    col = "green",
    ylab = "Original VCO2",
    xlab = "Time",
    type = "l"
  )
  par(new = TRUE)
  lines(
    org[,2],
    xlab = NA,
    ylab = NA,
    col = "red"
  )
  axis(side = 4)
  mtext('Smoothed VCO2', side = 4, line = 2)
  dev.off()
  jpeg(paste0(path, uid, '/',unlist(strsplit(basename(filenames), "[.]"))[1], "_VCO2_BlandAltman.jpeg"))
  plot(org[,2] - org[,1])
  dev.off()
  
  jpeg(paste0(path, uid, '/',unlist(strsplit(basename(filenames), "[.]"))[1], "_RER.jpeg"))
  org <- cbind(plotData$RER, lData$RER)
  plot(
    org[,1],
    col = "green",
    ylab = "Original RER",
    xlab = "Time",
    type = "l"
  )
  par(new = TRUE)
  lines(
    org[,2],
    xlab = NA,
    ylab = NA,
    col = "red"
  )
  axis(side = 4)
  mtext('Smoothed RER', side = 4, line = 2)
  dev.off()
  
  jpeg(paste0(path, uid, '/',unlist(strsplit(basename(filenames), "[.]"))[1], "_RER_BlandAltman.jpeg"))
  plot(org[,2] - org[,1])
  dev.off()
  
  #Validate data
  avglData <- lData %>%
    group_by(activity) %>%
    summarise(
      avgRER = mean(RER),
      avgsV.O2 = mean(VO2),
      avgsV.CO2 = mean(VCO2)
    )
  
  #Save the results as a CSV file
  fileName <-
    paste(unlist(strsplit(basename(filenames), "[.]"))[1], "_labeled.csv", sep = "")
  write.csv(lData, paste(path, uid, "/", fileName, sep = ""), row.names = FALSE)
  
  #Save validation as a XLSX file
  fileName <-
    paste0(unlist(strsplit(basename(filenames), "[.]"))[1], "_Validity.csv")
  write.csv(as.data.frame(avglData), paste0(path, uid, "/", fileName), row.names = FALSE)
}
