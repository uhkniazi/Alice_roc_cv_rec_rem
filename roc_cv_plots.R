# File: roc_cv_plots.R
# Auth: u.niazi@imperial.ac.uk
# Date: 27/03/2015
# Desc: plots of roc and cross validation of the data


### data loading
dfData = read.csv(file.choose(), header = TRUE)

### data cleaning and formatting
head(dfData)
# remove the column on ifn.sfc
dfData = dfData[,-5]
# create new variables
ivIl2 = dfData$IL2....
ivIfn = dfData$IFNy....
ivDual = dfData$Dual....
ivTnf = dfData$TNFa.Teff
ivCombined = dfData$Combined.score
fStatus = dfData$Status
ivAge = dfData$Age
ivTime = dfData$Time
fGen = dfData$Gender
# relevel recent/remote status as we need to predict only recents 
fStatus = relevel(fStatus, 'Remote')
dfData$Status = fStatus

####### for ROC of 10 points of combined score and others using lda
## TAG_1 
dfData.sub = dfData[,c(4, 5, 6, 8, 10, 7)]
colnames(dfData.sub) = c('status', 'il2', 'ifn', 'tnf', 'time', 'dual')
## change variable of choice to build model on
## NOTE: choose one variable i.e. one line of code below
# fit the model on each variable separately
dfData.sub = dfData.sub[,c(1, 2)] # il2
dfData.sub = dfData.sub[,c(1, 3)] # ifn
dfData.sub = dfData.sub[,c(1, 4)] # tnf
dfData.sub = dfData.sub[,c(1, 6)] # dual
# for combined scroe
dfData.sub = dfData[,c(4, 9)]
colnames(dfData.sub) = c('status', 'combined')

## continue from here after choosing one variable (see TAG_1)
dfData.sub = na.omit(dfData.sub)
library(MASS)
library(ROCR)
set.seed(1)
# set these variables for outer and inner loop respectively
iBoot.1 = 10
iBoot.2 = 10

# list variables to store results for ROC
lPred = vector(mode = 'list', length = iBoot.1)
lLab = vector(mode = 'list', length=iBoot.1)
iCv.error = NULL
## outer loop for multiple predictions, used for ROC curve confidence intervals
for (oo in 1:iBoot.1){
  t.lPred = NULL
  t.lLab = NULL
  # inner loop to perform nested cross validation
  # the data points are less so it artifically produced more points by resampling
  for (o in 1:iBoot.2){
    # perform 4 fold cross validation
    k = 4
    folds = sample(1:k, nrow(dfData.sub), replace = T, prob = rep(1/k, times=k))
    for (i in 1:k){
      # check if selected training and test folds have a 0 
      # observation for a particular class  
      if ((length(unique(dfData.sub$status[folds != i])) < 2) || (length(unique(dfData.sub$status[folds == i])) < 2)) next
      # check if number of observations in the selected fold too small to fit model
      if (nrow(dfData.sub[folds != i,]) < 3) next
      # if both previous checks have passed then continue or else 
      # resample the data
      # if all fine continue
      # fit model on data NOT in fold
      fit = lda(status ~ ., data=dfData.sub[folds != i,])
      # generate name for list
      name = paste('pred',oo, o, i,sep='.' )
      # predict on data IN fold for ROC curve
      pred = predict(fit, newdata = dfData.sub[folds == i,])$posterior[,'Recent']
      t.lPred[[name]] = pred
      # generate name for list
      name = paste('label',oo,o, i,sep='.' )
      # actual class label for the fold, use TRUE or FALSE for ROC
      t.lLab[[name]] = dfData.sub$status[folds == i] == 'Recent'
      # predict for data IN fold to calculate cross validation error
      pred = predict(fit, newdata = dfData.sub[folds == i,])$class
      iCv.error = append(iCv.error, mean(pred != dfData.sub$status[folds == i]))
    }
  }
  t.lPred = unlist(t.lPred)
  t.lLab = unlist(t.lLab)
  # save data in ROC list
  lPred[[oo]] = t.lPred
  lLab[[oo]] = t.lLab
}

## fit the ROC object
pred = prediction(lPred, lLab)
perf = performance(pred, 'tpr', 'fpr')
auc = performance(pred, 'auc')

## plot ROC, with confidence intervals
plot(perf, main=paste('ROC of Prediction by LDA of Recent based on', colnames(dfData.sub)[2]),
     spread.estimate='stddev', avg='vertical', spread.scale=2)
auc = paste('auc=', signif(mean(as.numeric(auc@y.values)), digits = 3))
cv = paste('CV Error=', signif(mean(iCv.error), 3))
legend('bottomright', legend = c(auc, cv))
abline(0, 1, lty=2)

## continue from here next time to get the actual values and cutoffs
# the cutoff values for prediction and performance objects
dfCutoffs = data.frame(tp=unlist(pred@tp), fp=unlist(pred@fp),  
                cutoff=unlist(perf@alpha.values), tpr=unlist(perf@y.values), fpr=unlist(perf@x.values))

