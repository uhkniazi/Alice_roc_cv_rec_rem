# Alice_roc_cv_rec_rem
Cross validation and ROC curve analysis for cytokine levels of Recent vs Remote groups

# roc_cv_plots.R
The script takes the csv file including the data, and creates factors and variable names (data cleaning). Look at TAG_1
for selecting the variable of choice on which we need to calculate the cross validation and roc curve error rates. The 
nested cross validation has 2 nested loops, the inner loop performs k fold cross validation, on the dataset. But it also
makes sure that there are at least 2 observations for each class in the fold and out of fold (this check is necessary as 
we do not have a lot of data points). The inner loop repeats this cross validation iBoot.2 times, using the LDA model,
and puts the test set results in the list lPred and lLab (predicted probabilities for Recent and original labels for the
data). The outer loop iBoot.1 will hold the results of each repetition of the cross validation run in lPred and lLab lists.
This is used to calculate confidence intervals for the ROC curves. The cross validation error, auc, tpr and fpr are calculated
and plotted. An additional line is added for the original values without using cross validation (as it doesn't make sense, 
because there are no predicted values). The cutoff values are also printed in the data frame. 

# roc_cv_plots_2.R
Very similar to previous script, but for new data sets