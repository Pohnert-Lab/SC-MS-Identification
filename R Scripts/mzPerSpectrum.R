set.seed(123)
library(ggplot2)
library(reshape2)
library(RCurl)
# Genus positive ----------------------------------------------------------

## Read spectra csv files

Algae <- read.csv(text=getURL("https://raw.github.com/Pohnert-Lab/SC-MS-Identification/master/Metadata/S11_metadata_md_positive.csv",.opts=curlOptions(followlocation = TRUE))) # S11

peak.list <- readRDS(gzcon(url("https://github.com/Pohnert-Lab/SC-MS-Identification/blob/master/RData_spectra/Spectra_S5_md_gen_HR_pos.RData?raw=true"))) # S5


brk <- seq(100,1000,by=10)

freq.peaks<-lapply(Algae$ID, function(x){   
  if(as.character(x) %in% names(peak.list)){   # Check if ID is part of used peak list
    
    mz<-peak.list[[as.character(x)]]$m.z
    
    histo<-hist(mz, breaks = brk, include.lowest = TRUE, plot = FALSE)
    
    return(frequency = histo$counts)
    
    }  # get number of peaks in spectrum for sp. ID  
  else NA                                      # if ID is not part of peak.list set NA
})       

names(freq.peaks)<-as.character(Algae$ID)

df.freq.peaks<-as.data.frame(freq.peaks)

ranges <- paste(head(brk,-1), brk[-1], sep=" - ")

df.freq.peaks<- as.data.frame(t(df.freq.peaks), stringsAsFactors = F)
 
df.freq.peaks<- data.frame(Genus = as.character(sapply(rownames(df.freq.peaks), function(x) Algae$Genus[which(as.character(Algae$ID) == x)])), df.freq.peaks)

df.freq.peaks[, 2]<- as.numeric(as.character(df.freq.peaks[,2]))
  
NrSpectraPerGroup<-aggregate(df.freq.peaks[1], by = list(df.freq.peaks$Genus), length)
    


agg<-aggregate(df.freq.peaks[2:91], by = list(df.freq.peaks$Genus), FUN = sum)  
  
colnames(agg)<-c("Genus", ranges)

agg[-1]<-agg[-1]/NrSpectraPerGroup$Genus

agg.melt<-melt(agg)

k<-rep(NULL, 90)

k<-ranges[rep(c(T,rep(F,10)), 15)] # 

p <- ggplot(agg.melt, aes(x = variable, y = value, fill = Genus))

p <- p + geom_bar(stat = "identity")

p <- p + ylab("Counts per spectrum") + scale_x_discrete(name="m/z window", breaks = k)

p<- p + facet_grid(Genus ~.)

p <- p + theme_bw()

plot(p)

# ggsave(filename = paste0("CountsOfMZperGroup_gen_pos.eps"),device = "eps", width = 210, height = 297, units = "mm")




# Genus negative ----------------------------------------------------------
## Read spectra csv files

Algae <- read.csv(text=getURL("https://raw.github.com/Pohnert-Lab/SC-MS-Identification/master/Metadata/S12_metadata_md_negative.csv",.opts=curlOptions(followlocation = TRUE))) # S12

peak.list <- readRDS(gzcon(url("https://github.com/Pohnert-Lab/SC-MS-Identification/blob/master/RData_spectra/Spectra_S6_md_gen_HR_neg.RData?raw=true"))) # S6

brk <- seq(100,1000,by=10)

freq.peaks<-lapply(Algae$ID, function(x){   
  if(as.character(x) %in% names(peak.list)){   # Check if ID is part of used peak list
    
    mz<-peak.list[[as.character(x)]]$m.z
    
    histo<-hist(mz, breaks = brk, include.lowest = TRUE, plot = FALSE)
    
    return(frequency = histo$counts)
    
  }  # get number of peaks in spectrum for sp. ID  
  else NA                                      # if ID is not part of peak.list set NA
})       

names(freq.peaks)<-as.character(Algae$ID)

df.freq.peaks<-as.data.frame(freq.peaks)

ranges <- paste(head(brk,-1), brk[-1], sep=" - ")

df.freq.peaks<- as.data.frame(t(df.freq.peaks), stringsAsFactors = F)

df.freq.peaks<- data.frame(Genus = as.character(sapply(rownames(df.freq.peaks), function(x) Algae$Genus[which(as.character(Algae$ID) == x)])), df.freq.peaks)

df.freq.peaks[, 2]<- as.numeric(as.character(df.freq.peaks[,2]))

NrSpectraPerGroup<-aggregate(df.freq.peaks[1], by = list(df.freq.peaks$Genus), length)



agg<-aggregate(df.freq.peaks[2:91], by = list(df.freq.peaks$Genus), FUN = sum)  

colnames(agg)<-c("Genus", ranges)

agg[-1]<-agg[-1]/NrSpectraPerGroup$Genus

agg.melt<-melt(agg)

k<-rep(NULL, 90)

k<-ranges[rep(c(T,rep(F,10)), 15)] # 

p <- ggplot(agg.melt, aes(x = variable, y = value, fill = Genus))

p <- p + geom_bar(stat = "identity")

p <- p + ylab("Counts per spectrum") + scale_x_discrete(name="m/z window", breaks = k)

p<- p + facet_grid(Genus ~.)

p <- p + theme_bw()

plot(p)

# ggsave(filename = paste0("CountsOfMZperGroup_gen_neg.eps"),device = "eps", width = 210, height = 297, units = "mm")





# Species positive --------------------------------------------------------
## Read spectra csv files

Algae <- read.csv(text=getURL("https://raw.github.com/Pohnert-Lab/SC-MS-Identification/master/Metadata/S13_metadata_csd_positive.csv",.opts=curlOptions(followlocation = TRUE))) # S13

peak.list <- readRDS(gzcon(url("https://github.com/Pohnert-Lab/SC-MS-Identification/blob/master/RData_spectra/Spectra_S1aS3_sp_HR_pos.RData?raw=true"))) # S1-S3


brk <- seq(100,1000,by=10)

freq.peaks<-lapply(Algae$ID, function(x){   
  if(as.character(x) %in% names(peak.list)){   # Check if ID is part of used peak list
    
    mz<-peak.list[[as.character(x)]]$m.z
    
    histo<-hist(mz, breaks = brk, include.lowest = TRUE, plot = FALSE)
    
    return(frequency = histo$counts)
    
  }  
  else NA                                    
})       

names(freq.peaks)<-as.character(Algae$ID)

df.freq.peaks<-as.data.frame(freq.peaks)

ranges <- paste(head(brk,-1), brk[-1], sep=" - ")

df.freq.peaks<- as.data.frame(t(df.freq.peaks), stringsAsFactors = F)

df.freq.peaks<- data.frame(Species = as.character(sapply(rownames(df.freq.peaks), function(x) Algae$Species[which(as.character(Algae$ID) == x)])), df.freq.peaks)

df.freq.peaks[, 2]<- as.numeric(as.character(df.freq.peaks[,2]))

NrSpectraPerGroup<-aggregate(df.freq.peaks[1], by = list(df.freq.peaks$Species), length)

agg<-aggregate(df.freq.peaks[2:91], by = list(df.freq.peaks$Species), FUN = sum)  

colnames(agg)<-c("Species", ranges)
agg.norm<-agg[-1]/NrSpectraPerGroup$Species

agg[-1]<-agg[-1]/NrSpectraPerGroup$Species

agg.melt<-melt(agg)

k<-rep(NULL, 90)

k<-ranges[rep(c(T,rep(F,10)), 15)] # 

p <- ggplot(agg.melt, aes(x = variable, y = value, fill = Species))

p <- p + geom_bar(stat = "identity")

p <- p + ylab("Counts per spectrum") + scale_x_discrete(name="m/z window", breaks = k)

p<- p + facet_grid(Species ~.)

p <- p + theme_bw()

p

# ggsave(filename = paste0("CountsOfMZperGroup_sp_pos.eps"),device = "eps", width = 210, height = 297, units = "mm")





# Species negative --------------------------------------------------------
## Read spectra csv files

Algae <- read.csv(text=getURL("https://raw.github.com/Pohnert-Lab/SC-MS-Identification/master/Metadata/S14_metadata_csd_negative.csv",.opts=curlOptions(followlocation = TRUE))) # S14

peak.list <- readRDS(gzcon(url("https://github.com/Pohnert-Lab/SC-MS-Identification/blob/master/RData_spectra/Spectra_S2aS4_sp_HR_neg.RData?raw=true"))) # S2aS4

brk <- seq(100,1000,by=10)

freq.peaks<-lapply(Algae$ID, function(x){   
  if(as.character(x) %in% names(peak.list)){   # Check if ID is part of used peak list
    
    mz<-peak.list[[as.character(x)]]$m.z
    
    histo<-hist(mz, breaks = brk, include.lowest = TRUE, plot = FALSE)
    
    return(frequency = histo$counts)
    
  }  
  else NA                                    
})       

names(freq.peaks)<-as.character(Algae$ID)

df.freq.peaks<-as.data.frame(freq.peaks)

ranges <- paste(head(brk,-1), brk[-1], sep=" - ")

df.freq.peaks<- as.data.frame(t(df.freq.peaks), stringsAsFactors = F)

df.freq.peaks<- data.frame(Species = as.character(sapply(rownames(df.freq.peaks), function(x) Algae$Species[which(as.character(Algae$ID) == x)])), df.freq.peaks)

df.freq.peaks[, 2]<- as.numeric(as.character(df.freq.peaks[,2]))

NrSpectraPerGroup<-aggregate(df.freq.peaks[1], by = list(df.freq.peaks$Species), length)

agg<-aggregate(df.freq.peaks[2:91], by = list(df.freq.peaks$Species), FUN = sum)  

colnames(agg)<-c("Species", ranges)
agg.norm<-agg[-1]/NrSpectraPerGroup$Species

agg[-1]<-agg[-1]/NrSpectraPerGroup$Species

agg.melt<-melt(agg)

k<-rep(NULL, 90)

k<-ranges[rep(c(T,rep(F,10)), 15)] # 

p <- ggplot(agg.melt, aes(x = variable, y = value, fill = Species))

p <- p + geom_bar(stat = "identity")

p <- p + ylab("Counts per spectrum") + scale_x_discrete(name="m/z window", breaks = k)

p<- p + facet_grid(Species ~.)

p <- p + theme_bw()

p

# ggsave(filename = paste0("CountsOfMZperGroup_sp_neg.eps"),device = "eps", width = 210, height = 297, units = "mm")

