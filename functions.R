# helper functions and libraries for heart_disease code
library(plyr)
library(knitr)
library(kableExtra)
library(pROC)
library(caTools)
library(pracma)
library(stargazer)
library(ggplot2)
library(rlist)

# get some prediction metrics about model
pred.power = function(model, testdata, printTable = FALSE, th=.5, optTH=FALSE) {
  predictTest = predict(model, type="response", newdata=testdata)
  t = table(predictTest > th, testdata$TenYearCHD)
  if (dim(t)[1] == 1) {t = rbind(t,c(0,0))}
  if(printTable){print(t);print(sum(t))}
  TP = t[2,2]
  TN = t[1,1]
  FP = t[2,1]
  FN = t[1,2]
  accuracy = (TP+TN)/sum(t)
  tpr = TP/(TP+FN) # sensitivity
  fnr = FN/(FN+TP)
  tnr = TN/(TN+FP)
  fpr = FP/(FP+TN)
  roc = roc(response = testdata$TenYearCHD, predictor = predictTest, plot=TRUE, quiet = TRUE)
  if (optTH == TRUE) {
    idx = which.max(roc$sensitivities + roc$specificities)
    cat("optimal threshold:",roc$thresholds[idx],"\n")
    cat("sensitivity:",roc$sensitivities[idx],"\n")
    cat("specificity:",roc$specificities[idx],"\n")
  }
  par(pty="s")
  plot(roc, grid=TRUE, print.auc=TRUE)
  metrics = data.frame("aic"=model$aic,"accuracy"=accuracy,"roc"=roc$auc,"tpr"=tpr,"tnr"=tnr,"fpr"=fpr,"fnr"=fnr)
  return(metrics)
}

roc_analysis = function(model, testdata, plotROC = FALSE, optTH = FALSE) {
  predictTest = predict(model, type="response", newdata=testdata)
  roc = roc(response = testdata$TenYearCHD, predictor = predictTest, plot=plotROC, quiet = TRUE)
  prevalence = count(testdata$TenYearCHD)$freq[2]/sum(count(testdata$TenYearCHD)$freq)
  roc$accuracy = (roc$sensitivities*prevalence) + (roc$specificities*(1-prevalence))
  if (optTH == TRUE) {
    idx = which.max(roc$sensitivities + roc$specificities)
    cat("optimal threshold:",roc$thresholds[idx],"\n")
    cat("sensitivity:",roc$sensitivities[idx],"\n")
    cat("specificity:",roc$specificities[idx],"\n")
    cat("sum:", roc$sensitivities[idx]+roc$specificities[idx],"\n")
    roc$optTHidx = idx
    
    sens75 = tail(which(abs(roc$sensitivities-.75)==min(abs(roc$sensitivities-.75))),n=1)
    cat("optimal threshold (sensitivity=.75):",roc$thresholds[sens75],"\nsensitivity:",
        roc$sensitivities[sens75],"\nspecificity:",roc$specificities[sens75],"\nsum:",
        roc$sensitivities[sens75]+roc$specificities[sens75])
    roc$sens75TH = roc$thresholds[sens75]
  }
  if (plotROC) {
    par(pty="s")
    plot(roc, grid=TRUE, print.auc=TRUE)
  }
  return(roc)
}

# custom proportion table
props = function(tbl){
  tbl[,2]/rowSums(tbl)}

# drop 1 test including BH multiple testing
d1 = function(model){
  d1.test = drop1(model, test="Chisq")
  d1.test$`Pr(>Chi)` = p.adjust(d1.test$`Pr(>Chi)`, method = "BH")
  d1.test = d1.test[order(d1.test$AIC),]
  print(d1.test)
}

# add 1 test incuding BH multiple testing
add1.test = function(model) {
  t = add1(model, scope = .~. + .^2, test="Chisq")
  t$`Pr(>Chi)` = p.adjust(t$`Pr(>Chi)`, method = "BH")
  t = t[order(t$`Pr(>Chi)`),]
  print(t)
}

# threshold analysis (sensitivity and specificity)
threshold.test = function(model, test, pltlab="", th.lims = c(.05,.95), npts = 50){
  th = linspace(th.lims[1],th.lims[2],npts)
  TP = FP = TN = FN = vector(mode = "numeric", length = npts)
  for (i in 1:npts) {
    m = pred.power(model,test,th = th[i])
    TP[i] = m$tpr; FP[i] = m$fpr; FN[i] = m$fnr; TN[i] = m$tnr
  }
  cols=c("violet","deepskyblue3","coral1","aquamarine3")
  plot(th,TP,lwd=2,col=cols[1],xlab = "Prediction Threshold",ylab = "",type = "l",ylim = c(0,1),main = pltlab)
  lines(th,TN,col=cols[2],lwd=2)
  legend("right",inset=.05,legend=c("Sensitivity","Specificity"),col=cols[1:2],lty = "solid",lwd=2)
  
}
