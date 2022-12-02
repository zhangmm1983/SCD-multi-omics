# Author: YuchuanLiu
# Email: 517549918@qq.com
# 2020.12.28


suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-i","--input"), type = "character", default = NULL,
    help = "gene/protein expression matrix file,e.g. fpkm.xls(force). "),
  make_option(c("-m","--mapping"), type = "character", default = NULL,
    help = "Group file of gene/protein expression-matrix. Colnames of group info should be 'Group'." ),
  make_option(c("-d","--method"), type = "character", default = "glm",
    help = "method choose from: glm/svm/rf. [default: glm]" ),
  make_option(c("-r","--repeatN"), type = "numeric", default = 500,
    help = "Repeat number of Training subset and ValidationData.[default: 500]"),
  make_option(c("-n","--topN"), type = "numeric", default = 4,
    help = "After the 500 times repeated analysis, the final consensus model was comprised of the combination of proteins which were selected the most(topN) in the 500 repeats .[default: 4]")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "/data/software/R/R-v3.5.1/bin/Rscript glm_ROC.R -i fpkm.xls -m mapping.txt -d glm -r 500 -n 4");
opt = parse_args(opt_parser);

suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(e1071))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(sampling))
suppressPackageStartupMessages(library(randomForest))

input <- opt$input
mapping <- opt$mapping
repeatN <- opt$repeatN   ### numbers of repeat
topN <- opt$topN   ###top N gene/protein
method <- "glm"
seed <- "20210730"

### subsample ####
FUN_sampling <- function(merged_file){
  sub_test <- strata(merged_file,stratanames=("Group"),size=round(table(merged_file$Group)*0.3),method="srswor")
  trainData <- merged_file[-sub_test$ID_unit,]
  rownames(trainData) <- trainData[,1]
  trainData <- trainData[,-1]
  validationData <- merged_file[sub_test$ID_unit,]
  rownames(validationData) <- validationData[,1]
  validationData <- validationData[,-1]
  newlist <- list(trainData,validationData)
  return(newlist)
}



FUN_glm <- function(trainData,validationData){
  LR_fit <- glm(Group~.,data=trainData,family = binomial(),control=list(maxit=100))
  #LR_fit = randomForest(Group ~., data=trainData_single)
  fit_AIC <- LR_fit[[11]]
  summ <- summary(LR_fit)
  fit_pvalue <- summ$coefficients[1,][[4]]
  LR_predict <- predict(LR_fit,newdata = validationData,type="response")
  ROC <- roc(validationData$Group,LR_predict, quiet = T)
  AUC <- round(as.numeric(gsub("Area under the curve: ","",ROC$auc)),4)
  info <- cbind(colnames(trainData)[2],fit_AIC,fit_pvalue,AUC)
  newlist <- list(info,ROC,LR_predict)
  return(newlist)
}



########  test  glm #########

file1 <- read.csv(input,sep="\t",quote="",header=T,row.names=1)
file1_t <- as.data.frame(t(file1))

mapping <- read.csv(mapping,sep="\t",quote="",header=T,row.names=1)
merged_file <- merge(mapping[,'Group',drop=F],file1_t,  by.x='row.names',by.y='row.names')


##### repeat ######


dir.create("repeat_modules")

if(method=="glm"){
  cnt <- 1
  module_repeat <- NULL
  while (cnt <= repeatN) {

    set.seed(cnt+10)
    newlist <- FUN_sampling(merged_file)
    trainData <- newlist[[1]]
    validationData <- newlist[[2]]

    AUC_matrix <- NULL
    for(i in 2:ncol(trainData)){
      trainData_single <- trainData[,c(1,i)]
      validationData_single <- validationData[,c(1,i)]
      res_info <- FUN_glm(trainData_single,validationData_single)[[1]]
      AUC_matrix <- rbind(AUC_matrix,res_info)
    }
    AUC_matrix <- as.data.frame(AUC_matrix)
    rownames(AUC_matrix) <- AUC_matrix[,1]
    colnames(AUC_matrix) <- c("Index","fit_AIC","fit_pvalue","AUC")

    AUC_matrix_sort <- AUC_matrix[order(AUC_matrix$AUC,decreasing = T),]   #### 
    AUC_matrix_selected <- AUC_matrix_sort[which(as.matrix(AUC_matrix_sort$AUC)>0.6),] #### 
    AUC_matrix_1 <- AUC_matrix_selected[which(as.matrix(AUC_matrix_selected$AUC)==1),]

    if(nrow(AUC_matrix_1)>0){
      print(paste("there were ",length(nrow(AUC_matrix_1)),"AUC 1 terms",sep=""))
    }else if(length(AUC_matrix_selected)==0){
      print("AUC of all index are <= 0.6")
    }else{
      classifer <- rownames(AUC_matrix_selected)[1]
      classifer_auc <- as.numeric(as.matrix(AUC_matrix_selected$AUC[1]))
      Number <- length(rownames(AUC_matrix_selected))
      if(Number <2){
        print("there were only one index selected(AUC > 0.6)")
      }else{
        AUC_matrix_repeatedly <- NULL
        for(j in 2:Number){
          rownames_list <- c("Group",rownames(AUC_matrix_selected)[1:j])
          trainData_repeatedly <- trainData[,rownames_list]
          validationData_repeatedly <- validationData[,rownames_list]
          res_info <- FUN_glm(trainData_repeatedly,validationData_repeatedly)[[1]]
          res_info <- as.data.frame(res_info)[,-1]
          rownames(res_info) <- str_c(rownames_list[2:(j+1)],collapse='@@')
          AUC_matrix_repeatedly <- rbind(AUC_matrix_repeatedly,res_info)
        }
        AUC_matrix_repeatedly <- as.data.frame(AUC_matrix_repeatedly)
        AUC_matrix_repeatedly_sort <- AUC_matrix_repeatedly[order(AUC_matrix_repeatedly$AUC,decreasing = T),]
        value <- classifer_auc + 0.01
  
        AUC_0.01more <- AUC_matrix_repeatedly_sort[which(as.matrix(AUC_matrix_repeatedly_sort$AUC)>= value),]
        Index <- rownames(AUC_0.01more) 
        AUC_0.01more_res <- cbind(Index,AUC_0.01more)
        write.table(AUC_0.01more_res,paste("repeat_modules/repeat",cnt,"_AUC_",value,"more.txt",sep=""),sep="\t",quote=F,row.names=F)
    }  
  }
  print(cnt)
  cnt <- cnt +1
  module_repeat <- c(module_repeat,rownames(AUC_0.01more))
  }  ## while循环终止

  all_list <- str_split(str_c(module_repeat,collapse='@@'),"@@")

  counts <- as.data.frame(table(all_list))
  counts_sort <- counts[order(counts$Freq,decreasing = T),]
  write.table(counts_sort,paste("repeat_",repeatN,"times_sumList.txt",sep=""),sep="\t",quote=F,row.names=F)
  print(head(counts_sort,10))


  
  
  AUClist <- NULL
  for(i in 1:2021){
  seed <- i

  FUN_sampling <- function(merged_file){
    sub_test <- strata(merged_file,stratanames=("Group"),size=round(table(merged_file$Group)*0.3),method="srswor")
    trainData <- merged_file[-sub_test$ID_unit,]
    rownames(trainData) <- trainData[,1]
    trainData <- trainData[,-1]
    validationData <- merged_file[sub_test$ID_unit,]
    rownames(validationData) <- validationData[,1]
    validationData <- validationData[,-1]
    newlist <- list(trainData,validationData)
    return(newlist)
   }


  set.seed(seed)
  newlist <- FUN_sampling(merged_file)
  trainData <- newlist[[1]]
  validationData <- newlist[[2]]

  comb_trainData <- trainData[,c("Group",as.matrix(counts_sort$all_list[1:topN]))]
  comb_validationData <- validationData[,c("Group",as.matrix(counts_sort$all_list[1:topN]))]
  

  LR_fit <-randomForest(Group ~., data=comb_trainData,importance = TRUE,ntree=1000,nodesize=5)

  LR_predict <- predict(LR_fit,newdata = comb_trainData,type="prob")
  ROC_train <- roc(comb_trainData$Group,LR_predict[,1], quiet = T)
  AUC_train <- round(as.numeric(gsub("Area under the curve: ","",ROC_train$auc)),4)

  LR_predict <- predict(LR_fit,newdata = comb_validationData,type="prob")
  ROC_valid <- roc(comb_validationData$Group,LR_predict[,1], quiet = T)
  AUC_valid <- round(as.numeric(gsub("Area under the curve: ","",ROC_valid$auc)),4)

  print(paste(AUC_train,AUC_valid,sep="\t"))

  seedAUC  <- cbind(seed,AUC_train,AUC_valid)

  AUClist <- rbind(seedAUC,AUClist)

  
}
write.table(AUClist,"AUClist.xls",sep="\t")
  


  all_list <- str_split(str_c(module_repeat,collapse='@@'),"@@")

  counts <- as.data.frame(table(all_list))
  counts_sort <- counts[order(counts$Freq,decreasing = T),]
  write.table(counts_sort,paste("repeat_",repeatN,"times_sumList.txt",sep=""),sep="\t",quote=F,row.names=F)
  print(head(counts_sort,10))



#### 
  set.seed(seed)
  newlist <- FUN_sampling(merged_file)
  trainData <- newlist[[1]]
  validationData <- newlist[[2]]
  test_predict <- newlist[[3]]

  comb_trainData <- trainData[,c("Group",as.matrix(counts_sort$all_list[1:topN]))]
  comb_validationData <- validationData[,c("Group",as.matrix(counts_sort$all_list[1:topN]))]

  res_info <- FUN_rf(comb_trainData,comb_validationData)[[1]]
  ROC <- FUN_rf(comb_trainData,comb_validationData)[[2]]
  AUC <- as.numeric(res_info[2])
}


