suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-i","--input"), type = "character", default = NULL,
    help = "gene/protein expression matrix file,e.g. fpkm.xls(force). "),
  make_option(c("-m","--mapping"), type = "character", default = NULL,
    help = "Group file of gene/protein expression-matrix. Colnames of group info should be 'Group'." ),
  make_option(c("-s","--seed"), type = "numeric", default = NULL,
    help = "Repeat number of Training subset and ValidationData.[default: 500]")
);
opt_parser = OptionParser(option_list=option_list,epilogue = "/data/software/R/R-v3.5.1/bin/Rscript glm_ROC.R -i fpkm.xls -m mapping.txt -s 100");
opt = parse_args(opt_parser);

suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(e1071))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(sampling))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

n <- opt$topN
seed <- opt$seed
input <- opt$input
mapping <- opt$mapping
method <- "rf"
#topNlist <- "repeat_500times_sumList.txt"


file1 <- read.csv(input,sep="\t",quote="",header=T,row.names=1)
file1_t <- as.data.frame(t(file1))

mapping <- read.csv(mapping,sep="\t",quote="",header=T,row.names=1)


merged_file <- merge(mapping[,'Group',drop=F],file1_t,  by.x='row.names',by.y='row.names')


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

merged_file2 <- merged_file
rownames(merged_file2) <- merged_file2[,1]
merged_file2 <- merged_file2[,-1]

LR_fit_all <-randomForest(Group ~., data=merged_file2,importance = TRUE,ntree=1000,nodesize=5)

LR_predict_all <- predict(LR_fit_all,newdata = merged_file2,type="prob")
ROC_all <- roc(merged_file2$Group,LR_predict_all[,1], quiet = T)
AUC_all <- round(as.numeric(gsub("Area under the curve: ","",ROC_all$auc)),4)
print(ci(ROC_all))
#print(t.test(ci(ROC_all)))
pdf(paste(method,"_ROC_all.pdf",sep=""))
plot.roc(ROC_all,main="AUC",col="#1c61b6")
text(0.7,0.8, paste("AUC is ",round(AUC_all,2),sep=""),cex=1)
legend("bottomright", legend=c("classifer"), col=c("#1c61b6"), lwd=2)
dev.off()


AUClist <- NULL
set.seed(seed)
newlist <- FUN_sampling(merged_file)
trainData <- newlist[[1]]
validationData <- newlist[[2]]

merged_file2 <- merged_file
rownames(merged_file2) <- merged_file2[,1]
merged_file2 <- merged_file2[,-1]

LR_fit <-randomForest(Group ~., data=trainData,importance = TRUE,ntree=1000,nodesize=5)

LR_predict_train <- predict(LR_fit,newdata = trainData,type="prob")
ROC_train <- roc(trainData$Group,LR_predict_train[,1], quiet = T)
AUC_train <- round(as.numeric(gsub("Area under the curve: ","",ROC_train$auc)),4)

LR_predict_valid <- predict(LR_fit,newdata = validationData,type="prob")
ROC_valid <- roc(validationData$Group,LR_predict_valid[,1], quiet = T)
AUC_valid <- round(as.numeric(gsub("Area under the curve: ","",ROC_valid$auc)),4)


XX <- coords(ROC_valid, "best", best.method = "youden",ret=c("threshold", "sensitivity","specificity","npv","ppv","tpr","fpr","tnr","fnr","fdr","accuracy","precision","youden"),transpose = F)
index <- rownames(XX)
res <- cbind(index,XX)
colnames(res)[1]=colnames(XX[1,])[1]
write.table(res,"YoudenIndex_all.xls",sep="\t",quote=F,row.names=F)

print(ci(ROC_valid))
print(t.test(ci(ROC_valid)))
print(ROC_valid$p.value)

print(paste(AUC_train,AUC_valid,sep="\t"))

seedAUC  <- cbind(seed,AUC_train,AUC_valid)

AUClist <- rbind(seedAUC,AUClist)
write.csv(AUClist,paste(seed,"_AUClist.xls",sep=""))

pdf(paste(method,"_ROC_train.pdf",sep=""))
plot.roc(ROC_train,main="AUC",col="#1c61b6")
text(0.7,0.8, paste("AUC is ",round(AUC_train,2),sep=""),cex=1)
legend("bottomright", legend=c("classifer"), col=c("#1c61b6"), lwd=2)
dev.off()

pdf(paste("best",method,"_ROC_train.pdf",sep=""))
plot.roc(ROC_train,print.thres = "best")
text(0.7,0.8, paste("AUC is ",round(AUC_train,2),sep=""),cex=1)
legend("bottomright", legend=c("classifer"), col=c("#1c61b6"), lwd=2)
dev.off()

pdf(paste(method,"_ROC_valid.pdf",sep=""))
plot.roc(ROC_valid,main="AUC",col="#1c61b6")
text(0.7,0.8, paste("AUC is ",round(AUC_valid,2),sep=""),cex=1)
legend("bottomright", legend=c("classifer"), col=c("#1c61b6"), lwd=2)
dev.off()

pdf(paste("best",method,"_ROC_valid.pdf",sep=""))
plot.roc(ROC_valid,print.thres = "best")
text(0.7,0.8, paste("AUC is ",round(AUC_valid,2),sep=""),cex=1)
legend("bottomright", legend=c("classifer"), col=c("#1c61b6"), lwd=2)
dev.off()

n <- ncol(validationData)
z_score <- scale(validationData[,2:n])
z_score_mean <- apply(z_score,1,mean)
dataxx <- data.frame(rownames(validationData),validationData$Group)
#####################################
#data <-data.frame(z_score_mean,LR_predict,dataxx)
data <-data.frame(z_score_mean,LR_predict_valid[,1],dataxx)
colnames(data) <- c("mean","predicts","sample2","type2")

pdf(paste("seed_",seed,"_105pair_validationData_plot.pdf",sep=""))
ggplot(data,aes(x=predicts,y=mean,group=type2,color=type2))+
geom_point(size=4)+
geom_vline(xintercept=0.5,linetype="dotted")+
xlab("predict value")+
ylab("Average molecular intersity")+
ggtitle("testData pointplot n = 64")+
xlim(0,1)+
theme(legend.text = element_text(size=15,color="black"),legend.position="right",
      legend.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour="black"))+
theme(panel.grid=element_blank())+
theme(axis.text = element_text(size=10,color="black"))+
theme(axis.text.x = element_text(hjust =1, angle =45))+
theme(plot.subtitle = element_text(size=30,hjust =0,color="black"))+
theme(axis.title.x = element_text(size=17,hjust =0.5,color="black"))+
theme(axis.title.y = element_text(size=17,hjust =0.5,color="black"))+
theme(plot.title = element_text(hjust = 0.5)) +
scale_color_manual(values=c("#E29827","#922927"))+
#scale_color_manual(limits=c("1","2"),values=c("#E29827","#922927"))+
geom_text_repel(aes(label=sample2,vjust=-0.8,hjust=0.5),color="black",show.legend=FALSE)
dev.off()


