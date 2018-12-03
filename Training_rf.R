#read feature matrix from RegDB
train_0 = read.delim("training_all_RegDB.txt",stringsAsFactors = FALSE)
train =  train_0[,-(1:8)]    #RegDB features

train[train$IC_change == ".","IC_change"] = 0
train[train$IC_matched_change == ".","IC_matched_change"] = 0
train$IC_change = as.numeric(train$IC_change)
train$IC_matched_change = as.numeric(train$IC_matched_change)
train$label = "no"
train[train_0$confidence >=0.1 & train_0$effect_size<0,"label"] = "neg"
train[train_0$confidence >=0.1 & train_0$effect_size>0,"label"] = "pos"
train$label = as.factor(train$label)
#read DeepSEA 919 diff scores & functiontal sig score
ds.train = read.delim("DeepSEA_features/jobs-training/infile.vcf.out.funsig.txt",header =FALSE, stringsAsFactors = FALSE)
colnames(ds.train) = c("chrom","pos","target","ref","alt","funsig")
ds.fun = data.frame(funsig = ds.train$funsig)
#impute NA in funsig with mean
ds.fun[is.na(ds.fun)] <- mean(na.omit(ds.fun[,1]))
ds.diff = read.delim("DeepSEA_features/jobs-training/infile.vcf.out.diff.txt",stringsAsFactors = FALSE)
ds.diff = ds.diff[,-(1:5)]
#imput NA in ds.diff with 0
ds.diff[is.na(ds.diff)] <- 0
train.multi = cbind.data.frame(train[,-10], ds.diff, ds.fun, sapply(levels(train$label), function(x) (x == train$label)))

train.multi$neg = as.factor(as.numeric(train.multi$neg))
train.multi$pos = as.factor(as.numeric(train.multi$pos))
train.multi$no = as.factor(as.numeric(train.multi$no))

#Train random forest models
library(mlr)
rf = makeLearner(cl = "classif.randomForest",par.vals = list(ntree=500)) #set parameter in rf
rf = setPredictType(rf, "prob")

trainTask.1 = makeClassifTask(data = train.multi, target = "neg",positive = 1)
trainTask.1 = dropFeatures(trainTask.1,c("pos","no"))
trainTask.2 = makeClassifTask(data = train.multi, target = "pos",positive = 1)
trainTask.2 = dropFeatures(trainTask.2,c("neg","no"))
trainTask.3 = makeClassifTask(data = train.multi, target = "no",positive = 1)
trainTask.3 = dropFeatures(trainTask.3,c("neg","pos"))
set.seed(1000)
rforest1 = train(rf,trainTask.1)
set.seed(1000)
rforest2 = train(rf,trainTask.2)
set.seed(1000)
rforest3 = train(rf,trainTask.3)

#Make predictions
test_0 = read.delim("test_all_RegDB.txt",stringsAsFactors = FALSE)
test = test_0[,-(1:6)]
test[test$IC_change == ".","IC_change"] = 0
test[test$IC_matched_change == ".","IC_matched_change"] = 0
test$IC_change = as.numeric(test$IC_change)
test$IC_matched_change = as.numeric(test$IC_matched_change)

ds.test_0 = read.delim("DeepSEA_features/jobs-Dataset/infile.vcf.out.diff.txt",header = TRUE,stringsAsFactors = FALSE)
ds.test = ds.test_0[,-(1:5)]
ds.test[is.na(ds.test)] = 0
ds.fun.test = read.delim("DeepSEA_features/jobs-Dataset/infile.vcf.out.funsig.txt",header = FALSE, stringsAsFactors = FALSE)
ds.fun.test = data.frame(funsig = ds.fun.test[,6])
ds.fun.test[is.na(ds.fun.test)] <- mean(na.omit(ds.fun.test[,1]))

test.multi = cbind.data.frame(test,ds.test,ds.fun.test)

pred1 = predict(rforest1, newdata = test.multi)
pred2 = predict(rforest2, newdata = test.multi)
pred3 = predict(rforest3, newdata = test.multi)

pred = data.frame(neg = getPredictionProbabilities(pred1),pos = getPredictionProbabilities(pred2), no = getPredictionProbabilities(pred3))

a = colnames(pred)[apply(pred,1,which.max)]
pred$label = a
pred$prob = apply(pred[,1:3],1,max)
pred[pred$label=="no","label"] = 0
pred[pred$label=="pos","label"] = 1
pred[pred$label=="neg","label"] = -1

#Train random forest regression models on confidence scores
rf = makeLearner(cl = "regr.randomForest",par.vals = list(ntree=500))
rf = setPredictType(rf,"se")

trainTask = makeRegrTask(data = cbind.data.frame(train.multi[,-(930:932)],confidence = train_0$confidence), target = "confidence")
set.seed(1000)
rforest = train(rf, trainTask)
pred.reg = predict(rforest,newdata = test.multi)
#add predicted confidence scores
pred = cbind.data.frame(pred, confidence = as.data.frame(getPredictionResponse(pred.reg)), se = as.data.frame(getPredictionSE(pred.reg)))
write.table(cbind.data.frame(test_0[,1:6], pred[,-(1:3)]),"final_prediction/submission.txt",quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)

#make continuous value prediction based on prob.
pred$effect = 0
pred[pred$label == 1,"effect"] = pred[pred$label == 1,"pos"]
pred[pred$label == -1,"effect"] = -pred[pred$label == -1,"neg"]
pred[pred$label==0 & pred$neg >= pred$pos,"effect"] = pred[pred$label==0 & pred$neg >= pred$pos,"no"] -1
pred[pred$label==0 & pred$neg < pred$pos,"effect"] = 1- pred[pred$label==0 & pred$neg < pred$pos,"no"]
pred$effect = round(pred$effect,3)
#write.table(cbind.data.frame(test_0[,1:6],pred),"../final_prediction/make_continuous_value/prob/submission_4_prob.txt",quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)





