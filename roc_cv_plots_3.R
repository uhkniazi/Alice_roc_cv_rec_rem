# Name: roc_cv_plots_3.R
# Auth: Umar Niazi u.niazi@imperial.ac.uk
# Date: 22/06/2016
# Desc: roc curves with cv for recent remote data.


if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/master/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

## redefine plotting function
setGeneric('plot.cv.performance', def = function(ob, legend.pos='bottomright', ...) standardGeneric('plot.cv.performance'))
setMethod('plot.cv.performance', signature='CCrossValidation.LDA', definition = function(ob, legend.pos='bottomright', ...){
  # plot the ROC curves for cross validation and validation set 
  # cv error
  pred = ob@oPred.cv
  perf = ob@oPerf.cv
  auc = ob@oAuc.cv
  plot(perf, main=paste('ROC Prediction of for', ob@cPred),
       spread.estimate='stddev', avg='vertical', spread.scale=2)
  auc.cv = paste('CV AUC=', round(mean(as.numeric(auc@y.values)), digits = 2))
  cv.err = paste('CV Error=', round(mean(ob@ivCV.error.rate), 2))
  abline(0, 1, lty=2)
  
#   # validation error
#   pred = ob@oPred.val
#   perf = ob@oPerf.val
#   auc = ob@oAuc.val
#   plot(perf, add=T, lty=3, lwd=2, col=2)
#   auc.t = paste('Val AUC=', round(mean(as.numeric(auc@y.values)), digits = 2))
#   err.t = paste('Val Error=', round(ob@iTest.error, 2))
  legend(legend.pos, legend = c(auc.cv, cv.err))
})

### import the data
dfData = read.csv('Data_external/import.csv', header=T, stringsAsFactors=F)

## format data for variable selection
fGroups = dfData$Group
dfData = data.frame(TNFaTEFF= dfData$TNFaTEFF)
fGroups = factor(fGroups, levels=c('Recent', 'Remote'))

## select test set
test = sample(1:nrow(dfData), nrow(dfData) * 0.2, replace = F)
table(fGroups[test]); table(fGroups[-test])

## 10 fold nested cross validation with various variable combinations
# try models of various sizes with CV
dfData.train = as.data.frame(dfData)
colnames(dfData.train) = colnames(dfData)
dfData.test = data.frame(dfData[test, ])
colnames(dfData.test) = colnames(dfData)

oCV = CCrossValidation.LDA(test.dat = (dfData.test), train.dat = (dfData.train), test.groups = fGroups[test],
                           train.groups = fGroups, level.predict = 'Recent', boot.num = 1000, k.fold = 10)

plot.cv.performance(oCV)
# print variable names and 95% confidence interval for AUC
x = getAUCVector(oCV)
sig = (signif(quantile(x, probs = c(0.025, 0.975)), 2))
print(sig)
summary(x)

# fit ROC on original values instead of using a model
pred.2 = prediction(dfData.train, fGroups == 'Recent')
perf.2 = performance(pred.2, 'tpr', 'fpr')
plot(perf.2, add=T, lty=4, lwd=2)

############# second data set
dfData = read.csv('Data_external/import_all.csv', header=T, stringsAsFactors=F)

## look at TB and Remote
dfData = dfData[dfData$Group != 'Recent',]
## format data for variable selection
fGroups = dfData$Group
dfData = data.frame(TNFaTEFF= dfData$TNFaTEFF)
fGroups = factor(fGroups, levels=c('ATB', 'Remote'))

## select test set
test = sample(1:nrow(dfData), nrow(dfData) * 0.2, replace = F)
table(fGroups[test]); table(fGroups[-test])

## 10 fold nested cross validation with various variable combinations
# try models of various sizes with CV
dfData.train = as.data.frame(dfData)
colnames(dfData.train) = colnames(dfData)
dfData.test = data.frame(dfData[test, ])
colnames(dfData.test) = colnames(dfData)

oCV = CCrossValidation.LDA(test.dat = (dfData.test), train.dat = (dfData.train), test.groups = fGroups[test],
                           train.groups = fGroups, level.predict = 'ATB', boot.num = 1000, k.fold = 10)

plot.cv.performance(oCV)
# print variable names and 95% confidence interval for AUC
x = getAUCVector(oCV)
sig = (signif(quantile(x, probs = c(0.025, 0.975)), 2))
print(sig)
summary(x)

## plot raw data
pred.2 = prediction(dfData.train, fGroups == 'ATB')
perf.2 = performance(pred.2, 'tpr', 'fpr')
plot(perf.2, add=T, lty=4, lwd=2)

########## third data set
dfData = read.csv('Data_external/import_all.csv', header=T, stringsAsFactors=F)

## look at Recent and Remote
dfData = dfData[dfData$Group != 'ATB',]
## format data for variable selection
fGroups = dfData$Group
dfData = data.frame(TNFaTEFF= dfData$TNFaTEFF)
fGroups = factor(fGroups, levels=c('Recent', 'Remote'))

## select test set
test = sample(1:nrow(dfData), nrow(dfData) * 0.2, replace = F)
table(fGroups[test]); table(fGroups[-test])

## 10 fold nested cross validation with various variable combinations
# try models of various sizes with CV
dfData.train = as.data.frame(dfData)
colnames(dfData.train) = colnames(dfData)
dfData.test = data.frame(dfData[test, ])
colnames(dfData.test) = colnames(dfData)

oCV = CCrossValidation.LDA(test.dat = (dfData.test), train.dat = (dfData.train), test.groups = fGroups[test],
                           train.groups = fGroups, level.predict = 'Recent', boot.num = 1000, k.fold = 10)

plot.cv.performance(oCV)
# print variable names and 95% confidence interval for AUC
x = getAUCVector(oCV)
sig = (signif(quantile(x, probs = c(0.025, 0.975)), 2))
print(sig)
summary(x)

## plot raw data
pred.2 = prediction(dfData.train, fGroups == 'Recent')
perf.2 = performance(pred.2, 'tpr', 'fpr')
plot(perf.2, add=T, lty=4, lwd=2)


##############################################################################
library(lattice)
dfData = read.csv('Data_external/import_all.csv', header=T, stringsAsFactors=T)

xyplot(TNFaTEFF ~ Group, data=dfData, type='o')

library(nlme)


