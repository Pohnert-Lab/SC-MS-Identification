set.seed(123)

##############
## Packages ##
##############

library(caret)
library(reshape2)
library(pROC)
require(devtools)
library(RCurl)

###############
## Functions ##
###############


## BacteriaMS functions from Yang et al. 2017 DOI:10.1021/acs.analchem.7b03820
devtools::source_url("https://github.com/lmsac/BacteriaMS/blob/master/functions/bootstrapping.R?raw=TRUE")
devtools::source_url("https://github.com/lmsac/BacteriaMS/blob/master/functions/preprocess.R?raw=TRUE")
devtools::source_url("https://github.com/lmsac/BacteriaMS/blob/master/functions/score.R?raw=TRUE")
devtools::source_url("https://github.com/lmsac/BacteriaMS/blob/master/functions/search_database.R?raw=TRUE")


## Altered BacteriaMS functions
get.bootstrap.score.TB = function(bootstrap.result, database = NULL, IDs = NULL, level = "Genus") {
  bootstrap.names = sapply(bootstrap.result, function(x) x[1]) # get the top hit (as ID) of every bt scoring (length is identical to number of bt generated)
  bootstrap.classes = unlist(sapply(bootstrap.names, function(x) database[[level]][which(IDs == x)], simplify = FALSE),use.names = FALSE) # Get Genus/Species name out of IDs
  scores = table(bootstrap.classes) / length(bootstrap.classes)
  c(names(scores),scores[[1]]) # slightly different output
}
get.bootstrap.confidence.score.TB = function(search.database.result, database = NULL, IDs = NULL, level = "Genus") {
  bootstrap.names = sapply(search.database.result$bootstrap, function(x) x[1]) # get the top hit (as ID) of every bt scoring (length is identical to number of bt generated)
  bootstrap.classes = unlist(sapply(bootstrap.names, function(x) database[[level]][which(IDs == x)], simplify = FALSE),use.names = FALSE) # Get Genus/Species name out of IDs
  besthit.scoring = database[[level]][which(IDs == search.database.result$sample[1])] #Get the top hit as name from top hit from cos/eu/ieu scoring
  sum(bootstrap.classes == besthit.scoring) / length(bootstrap.classes) # Calculate confidence score
}

######################
## Data preparation ##
######################

## File location to Microalgae_metadata.csv and corresponding spectra
Algae <- read.csv(text=getURL("https://raw.github.com/Pohnert-Lab/SC-MS-Identification/master/Metadata/S10_metadata_genus_negative.csv",.opts=curlOptions(followlocation = TRUE))) # URL S10

# FolderToCSVFiles<- "Folder to CSV files"
# csvFiles<- list.files(FolderToCSVFiles, pattern = ".csv", recursive = T, full.names = T)

## Get spectra from GitHub as RData
peak.list <- readRDS(gzcon(url("https://github.com/Pohnert-Lab/SC-MS-Identification/blob/master/RData_spectra/Spectra_S2_gen_HR_neg.RData?raw=true"))) # S2

## Normalize intensities of peak.list
ref.peaklist <- lapply(peak.list, normalize.intensity)

######################
## Database search  ##
######################

## Perform database search using ref.peaklist, bootstrap.times = 500, use of methods: Cos/Eu/iEu
top1result.list<-list()
start<-Sys.time()

for(i in 1:length(ref.peaklist)){
  
  
  intermediate.result<-  search.datebase.bootstrap(
    as.data.frame(ref.peaklist[[i]]),
    reference.spectra = ref.peaklist[-i], 
    tolerance = 5,#5ppm
    bootstrap.times = 500,
    method = c('eu', 'ieu', 'cosine'))
  # List of the three scores: eu, ieu and cos for sample and the x bootstrap (based on bootstrap.times) results
  
  
  result <- lapply(intermediate.result, function(res) {
    list(
      res$sample,
      genus = get.bootstrap.score.TB(
        res$bootstrap, database = Algae, IDs = Algae$ID, level = "Genus"
      ),
      species = get.bootstrap.score.TB(
        res$bootstrap, database = Algae, IDs = Algae$ID, level = "Species"
      ),
      genus.score = get.bootstrap.confidence.score.TB(
        res, database = Algae, IDs = Algae$ID, level = "Genus"
      ),
      species.score = get.bootstrap.confidence.score.TB(
        res, database = Algae, IDs = Algae$ID, level = "Species"
      )
    )
  })
  # Create Overview of results from databasesearch
  
  
  ID <- names(ref.peaklist[i])
  
  top1result.part<-data.frame(real.identity  = paste(Algae$Genus[which(Algae$ID==ID)], Algae$Species[which(Algae$ID==ID)],sep = "_"),
                              ID = ID,
                              sample.eu = result$eu[[1]][[1,1]],
                              sample.eu.val = as.numeric(result$eu[[1]][[1,2]]),
                              sample.ieu = result$ieu[[1]][[1,1]],
                              sample.ieu.val = as.numeric(result$ieu[[1]][[1,2]]),
                              sample.cos = result$cosine[[1]][[1,1]],
                              sample.cos.val = as.numeric(result$cosine[[1]][[1,2]]),
                              bt.genus.eu.score = result$eu$genus.score[1],
                              bt.genus.ieu.score = result$ieu$genus.score[1],
                              bt.genus.cos.score = result$cosine$genus.score[1],
                              stringsAsFactors = FALSE
  )
  
  top1result.list[[i]]<-top1result.part
  
  
  ## Progress bar
  cat("\014") 
  print(paste0(round(i/length(ref.peaklist)*100,1),"%"))
  print(Sys.time()-start)
}

## Get dataframe with top hits from all similarity measures

top1result.df <- do.call(rbind, top1result.list)

## Reference genus names

reference<-unlist(sapply(as.character(top1result.df$ID), function(x) Algae$Genus[which(Algae$ID == x)], simplify = FALSE),use.names = FALSE) 


## Format dataframe for ROC analysis
# Get genus names from metadata for top hit
tophits.genus <- unlist(sapply(as.character(top1result.df$ID), function(x) Algae$Genus[which(Algae$ID == x)], simplify = FALSE),use.names = FALSE) 

# Create vector with booleans of genus hits
hits.eu <-  tophits.genus == unlist(sapply(as.character(top1result.df$sample.eu), function(x) Algae$Genus[which(Algae$ID == x)], simplify = FALSE),use.names = FALSE)
hits.ieu <- tophits.genus == unlist(sapply(as.character(top1result.df$sample.ieu), function(x) Algae$Genus[which(Algae$ID == x)], simplify = FALSE),use.names = FALSE)
hits.cos <- tophits.genus == unlist(sapply(as.character(top1result.df$sample.cos), function(x) Algae$Genus[which(Algae$ID == x)], simplify = FALSE),use.names = FALSE)

# Create vector with names of genus hits
genus.eu  <-unlist(sapply(as.character(top1result.df$sample.eu), function(x) Algae$Genus[which(Algae$ID == x)], simplify = FALSE),use.names = FALSE)
genus.ieu <-unlist(sapply(as.character(top1result.df$sample.ieu), function(x) Algae$Genus[which(Algae$ID == x)], simplify = FALSE),use.names = FALSE)
genus.cos <-unlist(sapply(as.character(top1result.df$sample.cos), function(x) Algae$Genus[which(Algae$ID == x)], simplify = FALSE),use.names = FALSE)

# Combine booleans/genus names with top1result.df
top1result.df<- cbind(top1result.df, reference, genus.cos, hits.cos, genus.eu,hits.eu, genus.ieu, hits.ieu)

# write.csv(top1result.df, file = "top1result_genus_highres_negative_BT500.csv")

## If only data analysis from top1result dataframe shall be done, active lower command line.

# top1result.df<- read.csv(text=getURL("https://raw.github.com/Pohnert-Lab/SC-MS-Identification/master/Identification_results/Dataset_S2.csv", .opts=curlOptions(followlocation = TRUE))) # URL S2

##################
## ROC analysis ##
##################

# setEPS()
## ROC for identification of single microalgal cells with Cos/Eu/iEu on genus level
roc.cos<-roc(top1result.df$hits.cos, as.numeric(top1result.df$sample.cos.val))
roc.eu<-roc(top1result.df$hits.eu, as.numeric(top1result.df$sample.eu.val))
roc.ieu<-roc(top1result.df$hits.ieu, as.numeric(top1result.df$sample.ieu.val))

# postscript("ROC_Scores_genus_negative.eps", width = 8, height = 8)
plotr<-plot(roc.cos, col = "blue", main = "ROC - Score (Cos/Eu/iEu)", asp = 1)
plotr<-plot(roc.eu, add = TRUE, col = "red", asp = 1)
plotr<-plot(roc.ieu, add = TRUE, col = "green",asp = 1)
legend(x = 0.2, y = 0.2, c("Cos","Eu", "iEu"), col = c("blue", "red", "green"), pch = 19)
# dev.off()

ci.cos<-ci.auc(roc.cos)
ci.eu<-ci.auc(roc.eu)
ci.ieu<-ci.auc(roc.ieu)

## ROC for identification of single microalgal cells with Cos/Eu/iEu with bootstrapping assessment on genus level
roc.cos.bt<-roc(top1result.df$hits.cos, top1result.df$bt.genus.cos.score)
roc.eu.bt <-roc(top1result.df$hits.eu, top1result.df$bt.genus.eu.score)
roc.ieu.bt<-roc(top1result.df$hits.ieu, top1result.df$bt.genus.ieu.score)


# postscript("ROC_BT_genus_negative.eps", width = 8, height = 8)
plotr.bt<-plot(roc.cos.bt, col = "blue", main = "ROC - BT confidence score (cos/eu/ieu)")
plotr.bt<-plot(roc.eu.bt, add = TRUE, col = "red")
plotr.bt<-plot(roc.ieu.bt, add = TRUE, col = "green")
legend(x = 0.2, y = 0.2, c("Cos","Eu", "iEu"), col = c("blue", "red", "green"), pch = 19)
# dev.off()

ci.cos.bt<-ci.auc(roc.cos.bt)
ci.eu.bt<-ci.auc(roc.eu.bt)
ci.ieu.bt<-ci.auc(roc.ieu.bt)

############################################
## Score vs. sensitivity/error rate plots ##
############################################


sc.seq<-seq(0,1,by = 0.01) # Scores to test

## Cos

# postscript("SensErr_Cos_genus_negative.eps", width = 8, height = 8)
df.cos<-data.frame(sc = top1result.df$sample.cos.val, hits = top1result.df$hits.cos)

sensitivities <- sapply(sc.seq,function(x) length(df.cos$hits[df.cos$hits == T & df.cos$sc>=x]) / length(df.cos$hits[df.cos$hits == T]))
errorrates <- sapply(sc.seq,function(x) 1-length(df.cos$hits[df.cos$hits == T & df.cos$sc>=x]) / length(df.cos$hits[df.cos$sc >= x]))

sens.cos<-data.frame(scores = sc.seq, sensitivity = sensitivities)
sens.cos<-sens.cos[order(sens.cos$scores,decreasing = F),]

err.cos<-data.frame(scores = sc.seq, errorrates = errorrates) 
err.cos <- err.cos[order(err.cos$scores, decreasing = F),]

plot(sens.cos, type = "l", ylim=c(0,1),xlim = c(0,1), col = "green", ylab = "Value", main = "Sensitivity - Errorrate vs. Score - Cos") 
lines(err.cos, col = "red")


# add vertical line with errorrate less equal 10%

abline( v = err.cos$scores[min(which(err.cos$errorrates<=0.1))])
err.cos$scores[min(which(err.cos$errorrates<=0.1))] # Threshold
err.cos$errorrates[min(which(err.cos$errorrates<=0.1))] # Errorrate
sens.cos$sensitivity[which(sens.cos$scores == err.cos$scores[min(which(err.cos$errorrates<=0.1))])] # Sensitivity 
# dev.off()

## Score vs. Sensitivity/Errorrate plots
#Eu

# postscript("SensErr_Eu_genus_negative.eps", width = 8, height = 8)
df.eu<-data.frame(sc = top1result.df$sample.eu.val, hits = top1result.df$hits.eu)

sensitivities <- sapply(sc.seq,function(x) length(df.eu$hits[df.eu$hits == T & df.eu$sc>=x]) / length(df.eu$hits[df.eu$hits == T]))
errorrates <- sapply(sc.seq,function(x) 1-length(df.eu$hits[df.eu$hits == T & df.eu$sc>=x]) / length(df.eu$hits[df.eu$sc >= x]))

sens.eu<-data.frame(scores = sc.seq, sensitivity = sensitivities)
sens.eu<-sens.eu[order(sens.eu$scores,decreasing = F),]

err.eu<-data.frame(scores = sc.seq, errorrates = errorrates) 
err.eu <- err.eu[order(err.eu$scores, decreasing = F),]

plot(sens.eu, type = "l", ylim=c(0,1),xlim = c(0,1), col = "green", ylab = "Value", main = "Sensitivity - Errorrate vs. Score - Eu")
lines(err.eu, col = "red")

# add vertical line with errorrate less equal 10%

abline( v = err.eu$scores[min(which(err.eu$errorrates<=0.1))])
err.eu$scores[min(which(err.eu$errorrates<=0.1))] # Threshold
err.eu$errorrates[min(which(err.eu$errorrates<=0.1))] # Errorrate
sens.eu$sensitivity[which(sens.eu$scores == err.eu$scores[min(which(err.eu$errorrates<=0.1))])] # Sensitivity 
# dev.off()

## Score vs. Sensitivity/Errorrate plots
#iEu

# postscript("SensErr_iEu_genus_negative.eps", width = 8, height = 8)
df.ieu<-data.frame(sc = top1result.df$sample.ieu.val, hits = top1result.df$hits.ieu)

sensitivities <- sapply(sc.seq,function(x) length(df.ieu$hits[df.ieu$hits == T & df.ieu$sc>=x]) / length(df.ieu$hits[df.ieu$hits == T]))
errorrates <- sapply(sc.seq,function(x) 1-length(df.ieu$hits[df.ieu$hits == T & df.ieu$sc>=x]) / length(df.ieu$hits[df.ieu$sc >= x]))

sens.ieu<-data.frame(scores = sc.seq, sensitivity = sensitivities)
sens.ieu<-sens.ieu[order(sens.ieu$scores,decreasing = F),]

err.ieu<-data.frame(scores = sc.seq, errorrates = errorrates) 
err.ieu <- err.ieu[order(err.ieu$scores, decreasing = F),]

plot(sens.ieu, type = "l", ylim=c(0,1),xlim = c(0,1), col = "green", ylab = "Value", main = "Sensitivity - Errorrate vs. Score - iEu") 
lines(err.ieu, col = "red")

# add vertical line with errorrate less equal 10%

abline( v = err.ieu$scores[min(which(err.ieu$errorrates<=0.1))])
err.ieu$scores[min(which(err.ieu$errorrates<=0.1))] # Threshold
err.ieu$errorrates[min(which(err.ieu$errorrates<=0.1))] # Errorrate
sens.ieu$sensitivity[which(sens.ieu$scores == err.ieu$scores[min(which(err.ieu$errorrates<=0.1))])] # Sensitivity 
# dev.off()


## Score
###

## Score vs. Sensitivity/Errorrate plots 
# Cos with bootstrap assessment

# postscript("SensErr_BT_Cos_genus_negative.eps", width = 8, height = 8)
df.cos.bt<-data.frame(sc = top1result.df$bt.genus.cos.score, hits = top1result.df$hits.cos)

sensitivities <- sapply(sc.seq,function(x) length(df.cos.bt$hits[df.cos.bt$hits == T & df.cos.bt$sc>=x]) / length(df.cos.bt$hits[df.cos.bt$hits == T]))
errorrates <- sapply(sc.seq,function(x) 1-length(df.cos.bt$hits[df.cos.bt$hits == T & df.cos.bt$sc>=x]) / length(df.cos.bt$hits[df.cos.bt$sc >= x]))

sens.cos.bt<-data.frame(scores = sc.seq, sensitivity = sensitivities)
sens.cos.bt<-sens.cos.bt[order(sens.cos.bt$scores,decreasing = F),]

err.cos.bt<-data.frame(scores = sc.seq, errorrates = errorrates) 
err.cos.bt <- err.cos.bt[order(err.cos.bt$scores, decreasing = F),]

plot(sens.cos.bt, type = "l", ylim=c(0,1),xlim = c(0,1), col = "green", ylab = "Value", main = "Sensitivity - Errorrate vs. Score + BT - Cos ") 
lines(err.cos.bt, col = "red")

# add vertical line with errorrate less equal 10%

abline( v = err.cos.bt$scores[min(which(err.cos.bt$errorrates<=0.1))])
err.cos.bt$scores[min(which(err.cos.bt$errorrates<=0.1))] # Threshold
err.cos.bt$errorrates[min(which(err.cos.bt$errorrates<=0.1))] # Errorrate
sens.cos.bt$sensitivity[which(sens.cos.bt$scores == err.cos.bt$scores[min(which(err.cos.bt$errorrates<=0.1))])] # Sensitivity 
# dev.off()

## Score

## Score vs. Sensitivity/Errorrate plots
#Eu with bootstrap assessment

# postscript("SensErr_BT_Eu_genus_negative.eps", width = 8, height = 8)
df.eu.bt<-data.frame(sc = top1result.df$bt.genus.eu.score, hits = top1result.df$hits.eu)

sensitivities <- sapply(sc.seq,function(x) length(df.eu.bt$hits[df.eu.bt$hits == T & df.eu.bt$sc>=x]) / length(df.eu.bt$hits[df.eu.bt$hits == T]))
errorrates <- sapply(sc.seq,function(x) 1-length(df.eu.bt$hits[df.eu.bt$hits == T & df.eu.bt$sc>=x]) / length(df.eu.bt$hits[df.eu.bt$sc >= x]))

sens.eu.bt<-data.frame(scores = sc.seq, sensitivity = sensitivities)
sens.eu.bt<-sens.eu.bt[order(sens.eu.bt$scores,decreasing = F),]

err.eu.bt<-data.frame(scores = sc.seq, errorrates = errorrates) 
err.eu.bt <- err.eu.bt[order(err.eu.bt$scores, decreasing = F),]

plot(sens.eu.bt, type = "l", ylim=c(0,1),xlim = c(0,1), col = "green", ylab = "Value", main = "Sensitivity - Errorrate vs. Score + BT - Eu ")  
lines(err.eu.bt, col = "red")

# add vertical line with errorrate less equal 10%

abline( v = err.eu.bt$scores[min(which(err.eu.bt$errorrates<=0.1))])
err.eu.bt$scores[min(which(err.eu.bt$errorrates<=0.1))] # Threshold
err.eu.bt$errorrates[min(which(err.eu.bt$errorrates<=0.1))] # Errorrate
sens.eu.bt$sensitivity[which(sens.eu.bt$scores == err.eu.bt$scores[min(which(err.eu.bt$errorrates<=0.1))])] # Sensitivity 
# dev.off()

## Score vs. Sensitivity/Errorrate plots
#iEu with bootstrap assessment

# postscript("SensErr_BT_iEu_genus_negative.eps", width = 8, height = 8)
df.ieu.bt<-data.frame(sc = top1result.df$bt.genus.ieu.score, hits = top1result.df$hits.ieu)

sensitivities <- sapply(sc.seq,function(x) length(df.ieu.bt$hits[df.ieu.bt$hits == T & df.ieu.bt$sc>=x]) / length(df.ieu.bt$hits[df.ieu.bt$hits == T]))
errorrates <- sapply(sc.seq,function(x) 1-length(df.ieu.bt$hits[df.ieu.bt$hits == T & df.ieu.bt$sc>=x]) / length(df.ieu.bt$hits[df.ieu.bt$sc >= x]))

sens.ieu.bt<-data.frame(scores = sc.seq, sensitivity = sensitivities)
sens.ieu.bt<-sens.ieu.bt[order(sens.ieu.bt$scores,decreasing = F),]

err.ieu.bt<-data.frame(scores = sc.seq, errorrates = errorrates) 
err.ieu.bt <- err.ieu.bt[order(err.ieu.bt$scores, decreasing = F),]

plot(sens.ieu.bt, type = "l", ylim=c(0,1), col = "green", ylab = "Value", main = "Sensitivity - Errorrate vs. Score + BT - iEu ")  
lines(err.ieu.bt, col = "red")

# add vertical line with errorrate less equal 10%

abline( v = err.ieu.bt$scores[min(which(err.ieu.bt$errorrates<=0.1))])
err.ieu.bt$scores[min(which(err.ieu.bt$errorrates<=0.1))] # Threshold
err.ieu.bt$errorrates[min(which(err.ieu.bt$errorrates<=0.1))] # Errorrate
sens.ieu.bt$sensitivity[which(sens.ieu.bt$scores == err.ieu.bt$scores[min(which(err.ieu.bt$errorrates<=0.1))])] # Sensitivity 
# dev.off()

#######################################
## Confusion matrix / Accuracy plots ##
#######################################

confMatrices<-list()
for (i in c("cos","eu","ieu")){
  
  confmat<-confusionMatrix(data = top1result.df[[paste0("genus.",i)]], reference = top1result.df$reference) # Create confusion matrix based on input genera and top predictions from cos/eu/ieu scoring
  
  confmat.percent<-t(t(as.matrix(confmat$table))/summary(top1result.df$reference)) # Get accuracy matrix (in %). accuracy = hits/total amount specta
  
  longtable<-melt(confmat.percent)
  
  longtable$value<-longtable$value*100
  
  frac<-melt(matrix(paste(t(as.matrix(confmat$table)),summary(top1result.df$reference), sep = "/"),dim(as.matrix(confmat$table)),byrow = FALSE),na.rm = T) # Fraction matrix 
  
  frac<-unname(sapply(as.character(frac$value), function(x) {
    if(grepl(x, pattern = "^0/")){ 
      ""}
    else{x}
  }))
  
  
  p<-ggplot(data = longtable, aes (Reference, Prediction, fill = value))
  
  p <- p +
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "white", high = "darkorchid2", mid = "darkturquoise", 
                         midpoint = 50, limit = c(0,100), space = "Lab", 
                         name="Accuracy:") +
    theme(panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
    coord_fixed()+
    geom_text(aes(x = Prediction, y = Reference, label = frac))+
    ggtitle(paste("Similarity measure:",i))
  
  plot (p)
  
  # ggsave(filename = paste0("confmat_negative_",i,".eps"),device = "eps", width = 20, height = 20, units = "cm")

  confMatrices[[i]]<-confmat
  }


