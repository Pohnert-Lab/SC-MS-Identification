library(ggplot2)
library(RCurl)
set.seed(123)

# Genus HR positive -------------------------------------------------------
# Read spectra csv files
Algae <- read.csv(text=getURL("https://raw.github.com/Pohnert-Lab/SC-MS-Identification/master/Metadata/S10_metadata_md_positive.csv",.opts=curlOptions(followlocation = TRUE))) # S10

peak.list <- readRDS(gzcon(url("https://github.com/Pohnert-Lab/SC-MS-Identification/blob/master/RData_spectra/Spectra_S5_md_gen_HR_pos.RData?raw=true"))) # S5

npeaks<-sapply(Algae$ID, function(x){   
  if(as.character(x) %in% names(peak.list)){   # Check if ID is part of used peak list
    length(peak.list[[as.character(x)]]$m.z)}    # get number of peaks in spectrum for sp. ID  
  else NA                                      # if ID is not part of peak.list set NA
})             

peaks.df<-data.frame(ID = Algae$ID, Genus=Algae$Genus, npeaks = npeaks)

peaks.df<-subset(peaks.df, !is.na(peaks.df$npeaks)) #Remove spectra that are not part of the peak list

medians<-aggregate(npeaks ~ Genus, peaks.df, median)

# postscript(file = "npeaks_pos_gen.eps", family = "Courier") #Create Adobe-readable postscript. Courier can be used without a problem by my Adobe

## ggplot2 
p<- ggplot(peaks.df, aes(Genus, log10(npeaks), fill = Genus)) # left axis fluorescence

p<- p + geom_boxplot()

p<- p + geom_text(data = medians, aes(label = npeaks, y = log10(npeaks)))

# p<- p + geom_jitter(shape = 21, position = position_jitter(width = 0.2, height = 0))

p<- p + labs(x = "Microalgal genus",
             y = "Nr. of peaks (log10)")

p<- p + guides(fill = FALSE)

p<- p + theme_classic(base_size = 13)

plot(p)

# dev.off()




# Genus integer positive --------------------------------------------------

# Read spectra csv files
Algae <- read.csv(text=getURL("https://raw.github.com/Pohnert-Lab/SC-MS-Identification/master/Metadata/S10_metadata_md_positive.csv",.opts=curlOptions(followlocation = TRUE))) # S10

peak.list <- readRDS(gzcon(url("https://github.com/Pohnert-Lab/SC-MS-Identification/blob/master/RData_spectra/Spectra_S7_md_gen_int_pos.RData?raw=true"))) # S7

npeaks<-sapply(Algae$ID, function(x){   
  if(as.character(x) %in% names(peak.list)){   # Check if ID is part of used peak list
    length(peak.list[[as.character(x)]]$m.z)}    # get number of peaks in spectrum for sp. ID  
  else NA                                      # if ID is not part of peak.list set NA
})             

peaks.df<-data.frame(ID = Algae$ID, Genus=Algae$Genus, npeaks = npeaks)

peaks.df<-subset(peaks.df, !is.na(peaks.df$npeaks)) #Remove spectra that are not part of the peak list

medians<-aggregate(npeaks ~ Genus, peaks.df, median)

# postscript(file = "npeaks_pos_int_gen.eps", family = "Courier") #Create Adobe-readable postscript. Courier can be used without a problem by my Adobe

## ggplot2 
p<- ggplot(peaks.df, aes(Genus, log10(npeaks), fill = Genus)) # left axis fluorescence

p<- p + geom_boxplot()

p<- p + geom_text(data = medians, aes(label = npeaks, y = log10(npeaks)))

# p<- p + geom_jitter(shape = 21, position = position_jitter(width = 0.2, height = 0))

p<- p + labs(x = "Microalgal genus",
             y = "Nr. of peaks (log10)")

p<- p + guides(fill = FALSE)

p<- p + theme_classic(base_size = 13)

plot(p)

# dev.off()




# Species HR positive -----------------------------------------------------

## Read spectra csv files

Algae <- read.csv(text=getURL("https://raw.github.com/Pohnert-Lab/SC-MS-Identification/master/Metadata/S12_metadata_csd_positive.csv",.opts=curlOptions(followlocation = TRUE))) # S12

peak.list <- readRDS(gzcon(url("https://github.com/Pohnert-Lab/SC-MS-Identification/blob/master/RData_spectra/Spectra_S3_csd_sp_HR_pos.RData?raw=true"))) # S3

npeaks<-sapply(Algae$ID, function(x){   
  if(as.character(x) %in% names(peak.list)){   # Check if ID is part of used peak list
    length(peak.list[[as.character(x)]]$m.z)}    # get number of peaks in spectrum for sp. ID  
  else NA                                      # if ID is not part of peak.list set NA
})             

peaks.df<-data.frame(ID = Algae$ID, Species=paste(Algae$Genus, Algae$Species), npeaks = npeaks)

peaks.df<-subset(peaks.df, !is.na(peaks.df$npeaks)) #Remove spectra that are not part of the peak list

medians<-aggregate(npeaks ~ Species, peaks.df, median)

# postscript(file = "npeaks_pos_sp.eps", family = "Courier") #Create Adobe-readable postscript. Courier can be used without a problem by my Adobe

## ggplot2 
p<- ggplot(peaks.df, aes(Species, log10(npeaks), fill = Species)) # left axis fluorescence

p<- p + geom_boxplot()

p<- p + geom_text(data = medians, aes(label = npeaks, y = log10(npeaks)))

# p<- p + geom_jitter(shape = 21, position = position_jitter(width = 0.2, height = 0))

p<- p + labs(x = "Microalgal species",
             y = "Nr. of peaks (log10)")

p<- p + guides(fill = FALSE)

p<- p + theme_classic(base_size = 13)

plot(p)

# dev.off()





# Species integer positive ------------------------------------------------

## Read spectra csv files

Algae <- read.csv(text=getURL("https://raw.github.com/Pohnert-Lab/SC-MS-Identification/master/Metadata/S12_metadata_csd_positive.csv",.opts=curlOptions(followlocation = TRUE))) # S12

peak.list <- readRDS(gzcon(url("https://github.com/Pohnert-Lab/SC-MS-Identification/blob/master/RData_spectra/Spectra_S9_csd_sp_int_pos.RData?raw=true"))) # S9

npeaks<-sapply(Algae$ID, function(x){   
  if(as.character(x) %in% names(peak.list)){   # Check if ID is part of used peak list
    length(peak.list[[as.character(x)]]$m.z)}    # get number of peaks in spectrum for sp. ID  
  else NA                                      # if ID is not part of peak.list set NA
})             

peaks.df<-data.frame(ID = Algae$ID, Species=paste(Algae$Genus, Algae$Species), npeaks = npeaks)

peaks.df<-subset(peaks.df, !is.na(peaks.df$npeaks)) #Remove spectra that are not part of the peak list

medians<-aggregate(npeaks ~ Species, peaks.df, median)

# postscript(file = "npeaks_pos_int_sp.eps", family = "Courier") #Create Adobe-readable postscript. Courier can be used without a problem by my Adobe

## ggplot2 
p<- ggplot(peaks.df, aes(Species, log10(npeaks), fill = Species)) # left axis fluorescence

p<- p + geom_boxplot()

p<- p + geom_text(data = medians, aes(label = npeaks, y = log10(npeaks)))

# p<- p + geom_jitter(shape = 21, position = position_jitter(width = 0.2, height = 0))

p<- p + labs(x = "Microalgal species",
             y = "Nr. of peaks (log10)")

p<- p + guides(fill = FALSE)

p<- p + theme_classic(base_size = 13)

plot(p)

# dev.off()



# Genus HR negative -------------------------------------------------------

# Read spectra csv files
Algae <- read.csv(text=getURL("https://raw.github.com/Pohnert-Lab/SC-MS-Identification/master/Metadata/S11_metadata_md_negative.csv",.opts=curlOptions(followlocation = TRUE))) # URL S11

peak.list <- readRDS(gzcon(url("https://github.com/Pohnert-Lab/SC-MS-Identification/blob/master/RData_spectra/Spectra_S6_md_gen_HR_neg.RData?raw=true"))) # S6

npeaks<-sapply(Algae$ID, function(x){   
  if(as.character(x) %in% names(peak.list)){   # Check if ID is part of used peak list
  length(peak.list[[as.character(x)]]$m.z)}    # get number of peaks in spectrum for sp. ID  
  else NA                                      # if ID is not part of peak.list set NA
  })             

peaks.df<-data.frame(ID = Algae$ID, Genus=Algae$Genus, npeaks = npeaks)

peaks.df<-subset(peaks.df, !is.na(peaks.df$npeaks)) #Remove spectra that are not part of the peak list

medians<-aggregate(npeaks ~ Genus, peaks.df, median)

# postscript(file = "npeaks_neg_gen.eps", family = "Courier") #Create Adobe-readable postscript. Courier can be used without a problem by my Adobe

## ggplot2 
p<- ggplot(peaks.df, aes(Genus, log10(npeaks), fill = Genus)) # left axis fluorescence

p<- p + geom_boxplot()

p<- p + geom_text(data = medians, aes(label = npeaks, y = log10(npeaks)))

# p<- p + geom_jitter(shape = 21, position = position_jitter(width = 0.2, height = 0))

p<- p + labs(x = "Microalgal genus",
             y = "Nr. of peaks (log10)")

p<- p + guides(fill = FALSE)

p<- p + theme_classic(base_size = 13)

plot(p)

# dev.off()



# Genus integer negative --------------------------------------------------

# Read spectra csv files
Algae <- read.csv(text=getURL("https://raw.github.com/Pohnert-Lab/SC-MS-Identification/master/Metadata/S11_metadata_md_negative.csv",.opts=curlOptions(followlocation = TRUE))) # URL S11

peak.list <- readRDS(gzcon(url("https://github.com/Pohnert-Lab/SC-MS-Identification/blob/master/RData_spectra/Spectra_S8_md_gen_int_neg.RData?raw=true"))) # S8

npeaks<-sapply(Algae$ID, function(x){   
  if(as.character(x) %in% names(peak.list)){   # Check if ID is part of used peak list
    length(peak.list[[as.character(x)]]$m.z)}    # get number of peaks in spectrum for sp. ID  
  else NA                                      # if ID is not part of peak.list set NA
})             

peaks.df<-data.frame(ID = Algae$ID, Genus=Algae$Genus, npeaks = npeaks)

peaks.df<-subset(peaks.df, !is.na(peaks.df$npeaks)) #Remove spectra that are not part of the peak list

medians<-aggregate(npeaks ~ Genus, peaks.df, median)

# postscript(file = "npeaks_neg_int_gen.eps", family = "Courier") #Create Adobe-readable postscript. Courier can be used without a problem by my Adobe

## ggplot2 
p<- ggplot(peaks.df, aes(Genus, log10(npeaks), fill = Genus)) # left axis fluorescence

p<- p + geom_boxplot()

p<- p + geom_text(data = medians, aes(label = npeaks, y = log10(npeaks)))

# p<- p + geom_jitter(shape = 21, position = position_jitter(width = 0.2, height = 0))

p<- p + labs(x = "Microalgal genus",
             y = "Nr. of peaks (log10)")

p<- p + guides(fill = FALSE)

p<- p + theme_classic(base_size = 13)

plot(p)

# dev.off()






# Species HR negative -----------------------------------------------------

## Read spectra csv files

Algae <- read.csv(text=getURL("https://raw.github.com/Pohnert-Lab/SC-MS-Identification/master/Metadata/S13_metadata_csd_negative.csv",.opts=curlOptions(followlocation = TRUE))) # S13

peak.list <- readRDS(gzcon(url("https://github.com/Pohnert-Lab/SC-MS-Identification/blob/master/RData_spectra/Spectra_S2aS4_csd_sp_HR_neg.RData?raw=true"))) # S2aS4

npeaks<-sapply(Algae$ID, function(x){   
  if(as.character(x) %in% names(peak.list)){   # Check if ID is part of used peak list
    length(peak.list[[as.character(x)]]$m.z)}    # get number of peaks in spectrum for sp. ID  
  else NA                                      # if ID is not part of peak.list set NA
})             

peaks.df<-data.frame(ID = Algae$ID, Species=paste(Algae$Genus, Algae$Species), npeaks = npeaks)

peaks.df<-subset(peaks.df, !is.na(peaks.df$npeaks)) #Remove spectra that are not part of the peak list

medians<-aggregate(npeaks ~ Species, peaks.df, median)

# postscript(file = "npeaks_neg_sp.eps", family = "Courier") #Create Adobe-readable postscript. Courier can be used without a problem by my Adobe

## ggplot2 
p<- ggplot(peaks.df, aes(Species, log10(npeaks), fill = Species)) # left axis fluorescence

p<- p + geom_boxplot()

p<- p + geom_text(data = medians, aes(label = npeaks, y = log10(npeaks)))

# p<- p + geom_jitter(shape = 21, position = position_jitter(width = 0.2, height = 0))

p<- p + labs(x = "Microalgal species",
             y = "Nr. of peaks (log10)")

p<- p + guides(fill = FALSE)

p<- p + theme_classic(base_size = 13)

plot(p)

# dev.off()







# Species integer negative ------------------------------------------------

## Read spectra csv files

Algae <- read.csv(text=getURL("https://raw.github.com/Pohnert-Lab/SC-MS-Identification/master/Metadata/S12_metadata_csd_negative.csv",.opts=curlOptions(followlocation = TRUE))) # S12

peak.list <- readRDS(gzcon(url("https://github.com/Pohnert-Lab/SC-MS-Identification/blob/master/RData_spectra/Spectra_S8_sp_int_neg.RData?raw=true"))) # S8

npeaks<-sapply(Algae$ID, function(x){   
  if(as.character(x) %in% names(peak.list)){   # Check if ID is part of used peak list
    length(peak.list[[as.character(x)]]$m.z)}    # get number of peaks in spectrum for sp. ID  
  else NA                                      # if ID is not part of peak.list set NA
})             

peaks.df<-data.frame(ID = Algae$ID, Species=paste(Algae$Genus, Algae$Species), npeaks = npeaks)

peaks.df<-subset(peaks.df, !is.na(peaks.df$npeaks)) #Remove spectra that are not part of the peak list

medians<-aggregate(npeaks ~ Species, peaks.df, median)

# postscript(file = "npeaks_neg_int_species.eps", family = "Courier") #Create Adobe-readable postscript. Courier can be used without a problem by my Adobe

## ggplot2 
p<- ggplot(peaks.df, aes(Species, log10(npeaks), fill = Species)) # left axis fluorescence

p<- p + geom_boxplot()

p<- p + geom_text(data = medians, aes(label = npeaks, y = log10(npeaks)))

# p<- p + geom_jitter(shape = 21, position = position_jitter(width = 0.2, height = 0))

p<- p + labs(x = "Microalgal species",
             y = "Nr. of peaks (log10)")

p<- p + guides(fill = FALSE)

p<- p + theme_classic(base_size = 13)

plot(p)

# dev.off()


