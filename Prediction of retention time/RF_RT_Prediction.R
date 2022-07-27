library(randomForest)
set.seed(0)
#Training
RT_train <- read.csv("D:\\Training_set.csv")
RT_descriptor<- RT_train[1:nrow(RT_train),2:ncol(RT_train)]
RT_descriptor<-data.frame(RT_descriptor)
RT_time<-RT_train[1:nrow(RT_train),1]

#Feature screening
important_feature<-vector()
for(i in 1:100){  
  RT.rf<- randomForest(RT_time~., data=RT_descriptor, ntree=100)    
  score<-importance(RT.rf)[,1]
  rank<-order(score,decreasing = TRUE)[1:100]
  important_feature<-c(important_feature,rank)
  print(i)
}  

important_feature<-table(important_feature)
important_feature <- as.data.frame(important_feature)

freq<-important_feature[,2]
seleted_freq<-which(freq>50)
selected_feature <- as.numeric(as.character(important_feature[seleted_freq,1]))
print(length(selected_feature))
print(selected_feature)

RT_descriptor<- RT_descriptor[,selected_feature]
RT_descriptor<-data.frame(RT_descriptor)

#mtry and ntree optimization
featnum<-c(50)
tree_num<-c(1000)
err<-vector()
for(i in featnum){  
  for (j in tree_num) {
    mtry_test <- randomForest(RT_time~., data=RT_descriptor, mtry=i, ntree=j)    
    err<- c(err,mean( mtry_test$mse ) )
  }
}  
min_err_index<-which.min(err)
out_time<-floor((min_err_index-1)/length(tree_num))
opti_mtry<-featnum[out_time+1]
opti_ntree<-tree_num[min_err_index-length(tree_num)*out_time]
print(opti_mtry)
print(opti_ntree)

RT.rf<-randomForest(RT_time~., data=RT_descriptor, mtry=opti_mtry, ntree=opti_ntree)
pred<-predict(RT.rf,RT_descriptor)
for (i in pred) {
  write.table(i,"D:\\Training_set_prediction.csv",sep="",row.names = FALSE,col.names = FALSE,append = TRUE)
}
print(RT.rf)

#Prediction
RT_test <- read.csv("D:\\Test.csv")
RT_descriptor<- RT_test[1:nrow(RT_test),2:ncol(RT_test)]
RT_descriptor<- RT_descriptor[,selected_feature]
RT_descriptor<-data.frame(RT_descriptor)
RT_time<-RT_test[1:nrow(RT_test),1]
pred<-predict(RT.rf,RT_descriptor)
for (i in pred) {
  write.table(i,"D:\\Test_Prediction.csv",sep="",row.names = FALSE,col.names = FALSE,append = TRUE)
}
print('Finished.')