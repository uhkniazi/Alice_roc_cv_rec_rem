# File: roc_cv_plots_2.R
# Auth: u.niazi@imperial.ac.uk
# Date: 27/03/2015
# Desc: plots of roc and cross validation of the data on new data


### data loading
dfData = read.csv(file.choose(), header = TRUE)

### data cleaning and formatting
head(dfData)
# create new variables
ivIfn = dfData$IFN
ivTnf = dfData$TNFa
ivCombined = dfData$Combined
fStatus = dfData$Status
# relevel recent/remote status as we need to predict only recents 
fStatus = relevel(fStatus, 'Remote')
dfData$Status = fStatus

####### for ROC of data using LDA
## TAG_1 
dfData.sub = dfData[,2:5]
colnames(dfData.sub) = c('status', 'ifn', 'tnf', 'combined')
## change variable of choice to build model on
## NOTE: choose one variable i.e. one line of code below
# fit the model on each variable separately
dfData.sub = dfData.sub[,c(1, 2)] # ifn
dfData.sub = dfData.sub[,c(1, 3)] # tnf
dfData.sub = dfData.sub[,c(1, 4)] # combined

## continue from here after choosing one variable (see TAG_1)
dfData.sub = na.omit(dfData.sub)
library(MASS)
library(ROCR)
set.seed(1)
# set these variables for outer and inner loop respectively
iBoot.1 = 20
iBoot.2 = 100

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
     spread.estimate='stddev', avg='vertical', spread.scale=2, lwd=2)
auc = paste('auc=', signif(mean(as.numeric(auc@y.values)), digits = 3))
cv = paste('CV Error=', signif(mean(iCv.error), 3))
legend('bottomright', legend = c(auc, cv))
abline(0, 1, lty=2)

# fit ROC on original values instead of using a model
pred.2 = prediction(dfData.sub[,2], dfData.sub[,1] == 'Recent')
perf.2 = performance(pred.2, 'tpr', 'fpr')
plot(perf.2, add=T, lty=4, lwd=2)


# the cutoff values for prediction and performance objects 
dfCutoffs = data.frame(tp=unlist(pred.2@tp), fp=unlist(pred.2@fp),  
                       cutoff=unlist(perf.2@alpha.values), tpr=unlist(perf.2@y.values), fpr=unlist(perf.2@x.values))
