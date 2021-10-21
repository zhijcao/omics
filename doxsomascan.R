#
#Somascan Human Plasma data analysis
#
#library(outliers)
#library(beeswarm)
library(ggplot2)
library(magrittr)
library(mixOmics)
library(ggrepel)
library(survminer)
library(survival)
library(purrr)
library(tidyr)
library(dplyr)
library(openxlsx)
library(tibble)


library(plotROC)
library(OptimalCutpoints)
library(pROC)




rm(list=ls())

removeprotein <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/11-29-2016/soma 5 protein and calibrators remove.csv", sep=",", header=TRUE)
somaid <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/11-29-2016/somaid1.csv", sep=",", header=TRUE)
sample <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/sample_replace_duplicate_with_mean_04142018_subject33labchanged_time.csv")
#BCsample_ratio <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/BCsample_ratio.csv")
BCsample_ratio <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/BCsample_ratio_time.csv")
write.csv(BCsample_ratio, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/BCsample_ratio_time.csv")
BCsample_ratio <- inner_join(sample[,c(1:7)], BCsample_ratio[,-c(2:5)])

####PCA
length(names(sample)[-c(1:7)])
samplelog <- log(sample[,-c(1:7)])

BCsample_ratiolog <- log(BCsample_ratio[,-c(1:7)])
proteinname <- as.character(somaid$Target[somaid$SomaId %in% names(sample)[-c(1:7)]])

# tem <- data.frame(SomaId=names(samplelog), x=proteinname)
# tem1 <- inner_join(somaid, tem)  
# write.csv(tem1, "idmatch.csv")

names(samplelog) <- proteinname
names(BCsample_ratiolog) <- proteinname


pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/11082019/dox_somalog_pca_score_loading.pdf",
    width = 6, height = 5)
sample_pca <- pca(samplelog, scale=T)
plotIndiv(sample_pca, group=sample$TimePoint,  pch=sample$Status1, 
          legend = T, title = "pca")

#abnormal and normal
sample_pca <- pca(samplelog[sample$Status1=="Normal", ], scale=T)
plotIndiv(sample_pca,group=sample$TimePoint[sample$Status1=="Normal"],
          pch=sample$Status1[sample$Status1=="Normal"], 
          legend = T, title = "Normal.pca")

sample_pca <- pca(samplelog[sample$Status1=="Abnormal", ], scale=T)
plotIndiv(sample_pca,group=sample$TimePoint[sample$Status1=="Abnormal"],
          pch=sample$Status1[sample$Status1=="Abnormal"], 
                 legend = T, title = "Abnormal.pca")
#time point
sample_pca <- pca(samplelog[sample$TimePoint=="A", ], scale=T)
plotIndiv(sample_pca, ind.names = F, group=sample$Status1[sample$TimePoint=="A"], 
          legend = T, title = "A.pca")

sample_pca <- pca(samplelog[sample$TimePoint=="B", ], scale=T)
plotIndiv(sample_pca, ind.names = F, group=sample$Status1[sample$TimePoint=="B"], 
          legend = T, title = "B.pca")

sample_pca <- pca(samplelog[sample$TimePoint=="C", ], scale=T)
plotIndiv(sample_pca, ind.names = F, group=sample$Status1[sample$TimePoint=="C"], 
          legend = T, title = "C.pca")

dev.off()

###plsda
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/11082019/dox_somalog_plsda.pdf",
    width = 6, height = 5)
sample_plsda <- plsda(samplelog,Y=sample$TimePoint,  scale=T)
plotIndiv(sample_plsda, group=sample$TimePoint,  pch=sample$Status1, 
          legend = T, title = "~time posint plsda")
sample_plsda <- plsda(samplelog,Y=sample$Status1 ,  scale=T)
plotIndiv(sample_plsda, group=sample$TimePoint,  pch=sample$Status1, 
          legend = T, title = "~Status plsda")

###abnormal and normal
sample_plsda <- plsda(samplelog[sample$Status1=="Normal", ],Y=sample$TimePoint[sample$Status1=="Normal"],  scale=T)
plot(sample_plsda$loadings$X)
plot(sample_plsda$loadings.star[[1]])

plotIndiv(sample_plsda, group=sample$TimePoint[sample$Status1=="Normal"],  ind.names = F, 
          legend = T, title = "Normal plsda")
plotVar(sample_plsda, cutoff = 0.6, title = "Normal plsda")

sample_plsda <- plsda(samplelog[sample$Status1=="Abnormal", ],Y=sample$TimePoint[sample$Status1=="Abnormal"],  scale=T)
plotIndiv(sample_plsda, group=sample$TimePoint[sample$Status1=="Abnormal"],  pch=sample$Status1[sample$Status1=="Abnormal"], 
          legend = T, title = "Abnormal plsda")
plotVar(sample_plsda, cutoff = 0.6, title = "Abnormal plsda")

###timepoint
sample_plsda <- plsda(samplelog[sample$TimePoint=="A", ],Y=sample$Status1[sample$TimePoint=="A"],  scale=T)
plotIndiv(sample_plsda, ind.names = F,group=sample$Status1[sample$TimePoint=="A"],  
          legend = T, title = "A plsda")
plotVar(sample_plsda, cutoff = 0.6, title = "A plsda")
plot(sample_plsda$loadings$X)
plot(sample_plsda$loadings.star[[1]])

sample_plsda <- plsda(samplelog[sample$TimePoint=="B", ],Y=sample$Status1[sample$TimePoint=="B"],  scale=T)
plotIndiv(sample_plsda, ind.names = F,group=sample$Status1[sample$TimePoint=="B"],  
          legend = T, title = "B plsda")
plotVar(sample_plsda, cutoff = 0.6, title = "B plsda")
plot(sample_plsda$loadings$X)
plot(sample_plsda$loadings.star[[1]])

sample_plsda <- plsda(samplelog[sample$TimePoint=="C", ],Y=sample$Status1[sample$TimePoint=="C"],  scale=T)
plotIndiv(sample_plsda, ind.names = F,group=sample$Status1[sample$TimePoint=="C"],  
          legend = T, title = "C plsda")
plotVar(sample_plsda, cutoff = 0.6, title = "C plsda")

dev.off()

###sample Ratio timepoint
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/11082019/dox_somalog_ratio_plsda.pdf",
    width = 6, height = 5)
BCsample_ratio_plsda <- plsda(BCsample_ratiolog[BCsample_ratio$TimePoint=="B", ],Y=BCsample_ratio$Status1[BCsample_ratio$TimePoint=="B"],  scale=T)
plotIndiv(BCsample_ratio_plsda, ind.names = F,group=BCsample_ratio$Status1[BCsample_ratio$TimePoint=="B"],  
          legend = T, title = "B ratio plsda")

BCsample_ratio_plsda <- plsda(BCsample_ratiolog[BCsample_ratio$TimePoint=="C", ],Y=BCsample_ratio$Status1[BCsample_ratio$TimePoint=="C"],  scale=T)
plotIndiv(BCsample_ratio_plsda, ind.names = F,group=BCsample_ratio$Status1[BCsample_ratio$TimePoint=="C"],  
          legend = T, title = "C ratio plsda")

dev.off()

pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/11082019/dox_somalog_ratio_plsda_loading.pdf",
    width = 9, height = 11)
BCsample_ratio_plsda <- plsda(BCsample_ratiolog[BCsample_ratio$TimePoint=="B", ],Y=BCsample_ratio$Status1[BCsample_ratio$TimePoint=="B"],  scale=T)
plotLoadings(BCsample_ratio_plsda,comp=1, ndisplay=ndisplay, size.name=size.name, method="median",contrib="max", title="B ratio plsda comp1")
plotLoadings(BCsample_ratio_plsda,comp=2, ndisplay=ndisplay, size.name=size.name, method="median",contrib="max", title="B ratio plsda comp2")

BCsample_ratio_plsda <- plsda(BCsample_ratiolog[BCsample_ratio$TimePoint=="C", ],Y=BCsample_ratio$Status1[BCsample_ratio$TimePoint=="C"],  scale=T)
plotLoadings(BCsample_ratio_plsda,comp=1, ndisplay=ndisplay, size.name=size.name, method="median",contrib="max", title="C ratio plsda comp1")
plotLoadings(BCsample_ratio_plsda,comp=2, ndisplay=ndisplay, size.name=size.name, method="median",contrib="max", title="C ratio plsda comp2")

dev.off()

###plsda loading
ndisplay=100
size.name=0.5
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/11082019/dox_somalog_plsda loading1.pdf",
    width = 9, height = 11)

###abnormal and normal
sample_plsda <- plsda(samplelog[sample$Status1=="Normal", ],Y=sample$TimePoint[sample$Status1=="Normal"],  scale=T)

plotLoadings(sample_plsda,comp=1, ndisplay=ndisplay, size.name=size.name, method="median",contrib="max", title="Normal plsda comp1")
plotLoadings(sample_plsda,comp=2, ndisplay=ndisplay, size.name=size.name, method="median",contrib="max", title="Normal plsda comp2")

sample_plsda <- plsda(samplelog[sample$Status1=="Abnormal", ],Y=sample$TimePoint[sample$Status1=="Abnormal"],  scale=T)
plotLoadings(sample_plsda,comp=1, ndisplay=ndisplay, size.name=size.name, method="median",contrib="max", title="Abnormal plsda comp1")
plotLoadings(sample_plsda,comp=2, ndisplay=ndisplay, size.name=size.name, method="median",contrib="max", title="Abnormal plsda comp2")

###timepoint
sample_plsda <- plsda(samplelog[sample$TimePoint=="A", ],Y=sample$Status1[sample$TimePoint=="A"],  scale=T)
plotLoadings(sample_plsda,comp=1, ndisplay=ndisplay, size.name=size.name, method="median",contrib="max", title="A plsda comp1")
plotLoadings(sample_plsda,comp=2, ndisplay=ndisplay, size.name=size.name, method="median",contrib="max", title="A plsda comp2")

sample_plsda <- plsda(samplelog[sample$TimePoint=="B", ],Y=sample$Status1[sample$TimePoint=="B"],  scale=T)
?plotLoadings(sample_plsda,comp=1, ndisplay=ndisplay, size.name=size.name, method="median",contrib="max", title="B plsda comp1")
plotLoadings(sample_plsda,comp=2, ndisplay=ndisplay, size.name=size.name, method="median",contrib="max", title="B plsda comp2")

sample_plsda <- plsda(samplelog[sample$TimePoint=="C", ],Y=sample$Status1[sample$TimePoint=="C"],  scale=T)
plotLoadings(sample_plsda,comp=1, ndisplay=ndisplay, size.name=size.name, method="median",contrib="max", title="C plsda comp1")
plotLoadings(sample_plsda,comp=2, ndisplay=ndisplay, size.name=size.name, method="median",contrib="max", title="C plsda comp2")
selectVar(sample_plsda, comp=1)
sample_plsda$loadings$X
dev.off()





###linear regression
x <- names(sample)[1:20][8]
sample$Status1 <- factor(sample$Status1, levels = c("Normal", "Abnormal"))
str(sample)
contrast <- results$Contrast
re <- list()
for (x in names(sample)[-c(1:7)]){
  
  fml <- paste(x,"~ Status1+TimePoint + (1|SubjectID)")
  fml <-  as.formula(fml)
  fit <- lmer(fml, data=sample)
  summary(fit)
  results <- get_contrasts(fit, "Status1+TimePoint", adjust = "none")
  results$Contrast <- as.character(results$Contrast)
  rownames(results) <- results$Contrast
  results[[x]] <- p.adjust(results$p, method="BH")
  re[[x]] <- results[,x]
  
}

re1 <- cbind(results[, "Contrast"], data.frame(re))
re2 <- t(re1)
names(re1)
write.xlsx(re1, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/dox_mixed mode.xlsx")
str(re1)
names(re1) <- rownames(results1)

fdr <- function(x)p.adjust(x, method="BH")
re_fdr <- re1 %>% mutate_all(fdr)

names(re_fdr) <- paste(names(re_fdr), "FDR", sep = ": ")
names(re1) <- paste(names(re1), "p", sep = ": ")
stool_pfdr <- cbind(re1,re_fdr)

#####ROC analysis

#time point A intensity
data <- sample %>% filter(TimePoint=="A") %>% droplevels()
proteins <- names(data)[-c(1:7)]
roc_result <- data.frame()
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/dox_soma__A_intensity_roc plot_optimalcutoff_11202019.pdf", height = 4, width = 4)
i <- 1
for (protein in proteins){
  cat(i, ", ")
  dm <- data.frame(d=data[,"status"], m=data[,protein] )
  proteinname <- somaid$Target[somaid$SomaId==protein]
  
  rocoptimalcutpoint(dm=dm,title=paste(proteinname, ", A, Abnormal:"),n.cutoff=10,labelsize=3,method="Youden")
  
  roc_result<- rbind(roc_result,rocoptimalsesp(dm=dm,protein=proteinname, method="Youden"))
  
  i <- i+1
}
dev.off()
write.csv(roc_result, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/dox_soma_A_intensity_roc_aucsesp_10012019.csv")
roc_A_intensity_result <- roc_result


#time point B intensity
data <- sample %>% filter(TimePoint=="B") %>% droplevels()
proteins <- names(data)[-c(1:7)]
roc_result <- data.frame()
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/dox_soma__B_intensity_roc plot_optimalcutoff_11202019.pdf", height = 4, width = 4)
i <- 1
for (protein in proteins){
  cat(i, ", ")
  dm <- data.frame(d=data[,"status"], m=data[,protein] )
  proteinname <- somaid$Target[somaid$SomaId==protein]
  
  rocoptimalcutpoint(dm=dm,title=paste(proteinname, ", B, Abnormal:"),n.cutoff=10,labelsize=3,method="Youden")
  
  roc_result<- rbind(roc_result,rocoptimalsesp(dm=dm,protein=proteinname, method="Youden"))
  
  i <- i+1
}
dev.off()
write.csv(roc_result, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/dox_soma_B_intensity_roc_aucsesp_11202019.csv")
roc_B_intensity_result <- roc_result


#time point C intensity
data <- sample %>% filter(TimePoint=="C") %>% droplevels()
proteins <- names(data)[-c(1:7)]
roc_result <- data.frame()
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/dox_soma__C_intensity_roc plot_optimalcutoff_11202019.pdf", height = 4, width = 4)
i <- 1
for (protein in proteins){
  cat(i, ", ")
  dm <- data.frame(d=data[,"status"], m=data[,protein] )
  proteinname <- somaid$Target[somaid$SomaId==protein]
  
  rocoptimalcutpoint(dm=dm,title=paste(proteinname, ", C, Abnormal:"),n.cutoff=10,labelsize=3,method="Youden")
  
  roc_result<- rbind(roc_result,rocoptimalsesp(dm=dm,protein=proteinname, method="Youden"))
  
  i <- i+1
}
dev.off()
write.csv(roc_result, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/dox_soma_C_intensity_roc_aucsesp_11202019.csv")
roc_C_intensity_result <- roc_result





##
names(BCsample_ratio)[1:10]
proteins <- names(BCsample_ratio)[-c(1:7)]

#time point B
data <- BCsample_ratio %>% filter(TimePoint=="B") %>% droplevels()
#data[1:10, 1:10]
roc_result <- data.frame()
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/dox_soma__B_roc plot_optimalcutoff_10012019.pdf", height = 4, width = 4)
i <- 1
for (protein in proteins){
  cat(i, ", ")
   dm <- data.frame(d=data[,"status"], m=data[,protein] )
   proteinname <- somaid$Target[somaid$SomaId==protein]
  
  rocoptimalcutpoint(dm=dm,title=paste(proteinname, ", B, Abnormal:"),n.cutoff=10,labelsize=3,method="Youden")
  
  roc_result<- rbind(roc_result,rocoptimalsesp(dm=dm,protein=proteinname, method="Youden"))
  
  i <- i+1
}
dev.off()
write.csv(roc_result, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/dox_soma_B_roc_aucsesp_10012019.csv")

#time point C
data <- BCsample_ratio %>% filter(TimePoint=="C") %>% droplevels()
#data[1:10, 1:10]
roc_result <- data.frame()
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/dox_soma__C_roc plot_optimalcutoff_10012019.pdf", height = 4, width = 4)
i <- 1
for (protein in proteins){
  cat(i, ", ")
  dm <- data.frame(d=data[,"status"], m=data[,protein] )
  proteinname <- somaid$Target[somaid$SomaId==protein]
  
  rocoptimalcutpoint(dm=dm,title=paste(proteinname, ", C, Abnormal:"),n.cutoff=10,labelsize=3,method="Youden")
  
  roc_result<- rbind(roc_result,rocoptimalsesp(dm=dm,protein=proteinname, method="Youden"))
  
  i <- i+1
}
dev.off()
write.csv(roc_result, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/dox_soma_C_roc_aucsesp_10012019.csv")









set10 <- read.xlsx("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/FDA-18-001.hybNorm.plateScale.medNorm.calibrate.20180412.adat.xlsx", sheet="doxset10")
set10new <- read.xlsx("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/FDA-19-006 SomaLogic Returned Data/05-21-2019/FDA-19-001.hybNorm.plateScale.medNorm.calibrate.20190515.adat.xlsx", sheet="doxnewset10")
set10new239 <- set10new %>% filter(ID %in% set10$ID)
write.xlsx(set10new239, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/FDA-19-006 SomaLogic Returned Data/05-21-2019/FDA-19-001.hybNorm.plateScale.medNorm.calibrate.20190515.s239.xlsx")
set10new239 <- read.xlsx("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/FDA-19-006 SomaLogic Returned Data/05-21-2019/FDA-19-001.hybNorm.plateScale.medNorm.calibrate.20190515.s239.xlsx")

names(set10)
names(set10new239)
set10new239$ID
proteins <- names(set10new239)[-c(1:28)]
protein <- "CHIP"
rsquare <- list()
pdf("set10 vs set10new plot1.pdf")
for (protein in proteins) {
  lmre <- lm(set10[,protein]~ set10new239[, protein])
  sum <- summary(lmre)
  rsquare[protein] <- sum$r.squared
  plot (set10[,protein], set10new239[, protein], xlab="set10", ylab="set10new",
        main=paste(protein, format(sum$r.squared, digits = 4), sep=" R2: "))
  
  
}
dev.off()
write.csv(t(data.frame(rsquare)), "rsquare.csv")

plot(set10$CHIP, set10new239$CHIP)
sample <- intersect(set10$ID, set10new$ID)
length(sample)

# BAsample_ratio <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/Bsample_ratio.csv")
# CAsample_ratio <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/Csample_ratio.csv")

sampleA <- sample %>% filter(TimePoint=="A")
sampleB <- sample %>% filter(TimePoint=="B")
sampleC <- sample %>% filter(TimePoint=="C")

sampleB_ratio <- BCsample_ratio %>% filter(TimePoint=="B")
sampleC_ratio <- BCsample_ratio %>% filter(TimePoint=="C")
dim(sampleB_ratio)
sampleB_ratio$SampleId
sampleB_ratio$SL019100
categorize(sampleB_ratio$SL019100)

table(sampleB_ratio$TimePoint,categorize(sampleB_ratio$SL019100))
table(sampleA$TimePoint,categorize(sampleA$SL019100))

####survive analysis
# source("https://bioconductor.org/biocLite.R")
# biocLite("RTCGA.clinical")

#library(RTCGA.clinical)
Q25 <- function(x){quantile(x, probs=0.25)}
q25_50_75 <- function(x){data <- round(quantile(x, probs=c(0.25,0.5, 0.75), names=FALSE), digits=2)
                         return (paste(data[2],"(", data[1],"-", data[3], ")"))}
x <- c(1:100)/3
q25_50_75(x)
categorize(x)
categorize <- function(x){
  breaks <- unique(quantile(x,probs = c(0, 0.33333333,0.66666666, 1)))
  xx <- cut(x, breaks=breaks, include.lowest = TRUE, labels = c("P33", "P66", "P100"))
  #xx <- factor(xx, levels=c("p100", "P66","P33" ))
  return (xx)
} 

survivedat <- function(tem=sampleA){
  
  tem1<- log(tem[,-c(1:7)],2) %>%  map(~categorize(.x)) %>% data.frame()
  
  tem2 <- cbind(tem[,c("time","status")],  tem1)
  return (tem2)
}

sampleA[1:10,1:10]
dim(sampleA)
name <- "SL019100"
tier_median <- function(sampleA=sampleA){
  re <- list()
  protein <- names(sampleA)[-c(1:7)]
  for (name in protein){
    tem <- data.frame(tiers=categorize(sampleA[,name]), sampleA[name])
    re[[name]]  <- tem %>% group_by(tiers) %>% summarise_all(funs(q25_50_75)) %>% 
      gather(key="SomaId", value = "Median", - tiers)
  }
  result <- do.call(rbind, re) 
  result1 <- result%>% spread(key=tiers, value = Median) %>% select( SomaId,  P33, P66,P100)
  result2 <- inner_join(somaid,result1)
  
  return (result2)
  
}

sampleA_tier_median <- tier_median(sampleA=sampleA)
sampleB_ratio_tier_median <- tier_median(sampleA=sampleB_ratio)
sampleC_ratio_tier_median <- tier_median(sampleA=sampleC_ratio)

write.xlsx(list(sampleA_tier_median=sampleA_tier_median,sampleB_ratio_tier_median=sampleB_ratio_tier_median,
                sampleC_ratio_tier_median=sampleC_ratio_tier_median), "tier_25_50_75_09252018.xlsx")


names(re[[1]])
class(tem_sum)

names(sampleB_ratio)[1:10]


pro$SL019100
fit <- survfit(Surv(time, status)~SL019100, pro)
cox <- coxph(Surv(time, status)~SL019100, sampleA_time)
tem <- summary(cox)
tem$coefficients[2,1]
length(tem$coefficients[,-2])
tem$conf.int[,-2]
tem1 <- cbind(tem$conf.int,tem$coefficients)
summary(fit)
dim(sampleA)
tem$waldtest
tem$concordance
tem$rsq
tem2 <- c(tem$waldtest,tem$concordance, tem$rsq)

tem3 <- list(a=tem1,b=tem2)
tem3$a
somaid$SomaId

####coxph analysis
sampleA_time <- survivedat(tem=sampleA)
sampleA_cox <- univ_cox(data=sampleA_time)

df <- data.frame(sampleA_cox$betahr, check.names = FALSE) %>% rownames_to_column("SomaId") %>% 
  mutate(SomaId=toupper(SomaId)) %>% separate(SomaId, into = c("SomaId","Tiers"), sep ="P")
A_betahr <-  inner_join(somaid, df)
A_test <- inner_join(somaid, data.frame(sampleA_cox$test, check.names = FALSE) %>% rownames_to_column("SomaId"))
  
sampleB_ratio_time <- survivedat(tem=sampleB_ratio)
sampleB_ratio_cox <- univ_cox(data=sampleB_ratio_time)
df <- data.frame(sampleB_ratio_cox$betahr, check.names = FALSE) %>% rownames_to_column("SomaId") %>% 
  mutate(SomaId=toupper(SomaId)) %>% separate(SomaId, into = c("SomaId","Tiers"), sep ="P")
B_ratio_betahr <-  inner_join(somaid, df)
B_ratio_test <- inner_join(somaid, data.frame(sampleB_ratio_cox$test, check.names = FALSE) %>% rownames_to_column("SomaId"))

sampleC_ratio_time <- survivedat(tem=sampleC_ratio)
sampleC_ratio_cox <- univ_cox(data=sampleC_ratio_time)
df <- data.frame(sampleC_ratio_cox$betahr, check.names = FALSE) %>% rownames_to_column("SomaId") %>% 
  mutate(SomaId=toupper(SomaId)) %>% separate(SomaId, into = c("SomaId","Tiers"), sep ="P")
C_ratio_betahr <-  inner_join(somaid, df)
C_ratio_test <- inner_join(somaid, data.frame(sampleC_ratio_cox$test, check.names = FALSE) %>% rownames_to_column("SomaId"))


write.xlsx(list(A_betahr=A_betahr,A_test=A_test,
                B_ratio_betahr=B_ratio_betahr,B_ratio_test=B_ratio_test,
                C_ratio_betahr=C_ratio_betahr, C_ratio_test=C_ratio_test), 
           "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/dox_coxph_result_p100as_control_09212018.xlsx", overwrite=FALSE)

dim(CAsample_ratio)

names(pro)[1:20]
pro[,c(1,2,48)]
i <- 3
i <- "SL012769"
sampleA_time<- survivedat(tem=sampleA)

pvalue_sampleA_dftest<- do.call(rbind,pvalue_sampleAtest)

sampleA_time<- survivedat(tem=sampleA)
sampleA_time$SL019100

pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/sampleA_cumevent_plot3.pdf", width = 6, height = 6)
data=sampleA_time
pvalue_sampleA <- list()
for (i in 3:1307){
  
  #fit <- survfit(as.formula(paste('Surv(time, status)~', name)), data = data)
  fit <- survfit(Surv(time, status)~ data[,i], data = data)
  
  target <- somaid$Target[somaid$SomaId==names(data)[i]]
  
  pvalue_sampleA[[paste(target)]] <- surv_pvalue(fit, data=data)$pval
  
  title=paste("strata:",target)
  # Visualize with survminer
  p1 <- ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    data =data,  # data used to fit survival curves. 
    risk.table = TRUE,       # show risk table.
    pval = TRUE,             # show p-value of log-rank test.
    fun="event",
    #conf.int = TRUE,         # show confidence intervals for 
    # point estimaes of survival curves.
    ylim = c(0,0.8),        # present narrower X axis, but not affect
    # survival estimates.
    legend.title=title,
    legend.labs=c("P33","P66","P100"),
    break.time.by = 1,     # break X axis in time intervals by 500.
    ggtheme = theme_minimal(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
  )
  
  print (p1)
  cat(i," ")
  #pvalue[[i]]
  
}
dev.off()
pvalue_sampleA_df<- do.call(rbind,pvalue_sampleA)
write.csv(pvalue_sampleA_df,"C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/pvalue_sampleA.csv")

sampleB_ratio_time<- survivedat(tem=sampleB_ratio)
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/sampleB_ratio_cumevent_plot4_09212018.pdf", width = 6, height = 6)
data=sampleB_ratio_time
pvalue_sampleB_ratio <- list()
for (i in 3:1307){
  
  #fit <- survfit(as.formula(paste('Surv(time, status)~', name)), data = data)
  fit <- survfit(Surv(time, status)~ data[,i], data = data)
  
  target <- somaid$Target[somaid$SomaId==names(data)[i]]
  
  pvalue_sampleB_ratio[[paste(target)]] <- surv_pvalue(fit, data=data)$pval
  
  title=paste("strata:",target)
  # Visualize with survminer
  p1 <- ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    data =data,  # data used to fit survival curves. 
    risk.table = TRUE,       # show risk table.
    pval = TRUE,             # show p-value of log-rank test.
    fun="event",
    #conf.int = TRUE,         # show confidence intervals for 
    # point estimaes of survival curves.
    ylim = c(0,0.8),        # present narrower X axis, but not affect
    # survival estimates.
    legend.title=title,
    legend.labs=c("P33","P66","P100"),
    break.time.by = 1,     # break X axis in time intervals by 500.
    ggtheme = theme_minimal(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
  )
  
  print (p1)
  cat(i," ")
  #pvalue[[i]]
  
}

dev.off()
pvalue_sampleB_ratio_df<- do.call(rbind,pvalue_sampleB_ratio)
write.csv(pvalue_sampleB_ratio_df,"C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/pvalue_sampleB_ratio_09212018.csv")


sampleC_ratio_time<- survivedat(tem=sampleC_ratio)
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/sampleC_ratio_cumevent_plot3.pdf", width = 6, height = 6)
data=sampleC_ratio_time
pvalue_sampleC_ratio <- list()
for (i in 3:1307){
  
  #fit <- survfit(as.formula(paste('Surv(time, status)~', name)), data = data)
  fit <- survfit(Surv(time, status)~ data[,i], data = data)
  
  target <- somaid$Target[somaid$SomaId==names(data)[i]]
  
  pvalue_sampleC_ratio[[paste(target)]] <- surv_pvalue(fit, data=data)$pval
  
  title=paste("strata:",target)
  # Visualize with survminer
  p1 <- ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    data =data,  # data used to fit survival curves. 
    risk.table = TRUE,       # show risk table.
    pval = TRUE,             # show p-value of log-rank test.
    fun="event",
    #conf.int = TRUE,         # show confidence intervals for 
    # point estimaes of survival curves.
    ylim = c(0,0.8),        # present narrower X axis, but not affect
    # survival estimates.
    legend.title=title,
    legend.labs=c("P33","P66","P100"),
    break.time.by = 1,     # break X axis in time intervals by 500.
    ggtheme = theme_minimal(), # customize plot and risk table with a theme.
    risk.table.y.text.col = T, # colour risk table text annotations.
    risk.table.y.text = FALSE # show bars instead of names in text annotations
    # in legend of risk table
  )
  
  print (p1)
  cat(i," ")
  #pvalue[[i]]
  
}

dev.off()
pvalue_sampleC_ratio_df<- do.call(rbind,pvalue_sampleC_ratio)
write.csv(pvalue_sampleC_ratio_df,"C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/pvalue_sampleC_ratio.csv")

pvalue_all <- data.frame(SampleA=pvalue_sampleA_df, smapleB_ratio=pvalue_sampleB_ratio_df,sampleC_ratio=pvalue_sampleC_ratio_df)
write.csv(pvalue_all,"C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/logrank_pvalues_timepointABC.csv" )


pvalue_all <- pvalue_sampleC_ratio_df$pval


i <- 1
name <- "SL019100"
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/test5d5.pdf", width = 5.5, height = 5.5)

survplot <- function(data=pro){
  pvalue <- list()
  

for (i in 3:6){

#fit <- survfit(as.formula(paste('Surv(time, status)~', name)), data = data)
fit <- survfit(Surv(time, status)~ data[,i], data = data)

target <- somaid$Target[somaid$SomaId==names(data)[i]]

title=paste("strata:",target)
# Visualize with survminer
p1 <- ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  data =data,  # data used to fit survival curves. 
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  fun="event",
  #conf.int = TRUE,         # show confidence intervals for 
  # point estimaes of survival curves.
  ylim = c(0,0.8),        # present narrower X axis, but not affect
  # survival estimates.
  legend.title=title,
  legend.labs=c("P33","P66","P100"),
   break.time.by = 1,     # break X axis in time intervals by 500.
  ggtheme = theme_minimal(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

print (p1)
#pvalue[[paste(target)]] <- surv_pvalue(fit, data=data)
cat(i," ")
#pvalue[[i]]

}
  #return (pvalue)
}



dev.off()









###################################
####09142018 data  analysis
###################################

sample <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/sample_replace_duplicate_with_mean_04142018_subject33labchanged.csv")

CV <- function(x){sd(x)/mean(x)*100}
Q25 <- function(x){quantile(x, probs=0.25)}
Q75 <- function(x){quantile(x, probs=0.75)}

sample_summary <- allsample[,-c(1:4)] %>% gather(key="ID", value="value", -Time.Status) %>% group_by(ID, Time.Status) %>% 
  summarise_all(funs(mean, sd, CV, median,Q25,  Q75 ))
sample_summary1 <- inner_join(somaid,sample_summary, by= c("SomaId" ="ID" ))

write.csv(sample_summary1, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/dox_summary_09142018.csv")


sample_ratio_summary <- BCsample_ratio[,-c(1:4)] %>% gather(key="ID", value="value", -Time.Status) %>% group_by(ID, Time.Status) %>% 
  summarise_all(funs(mean, sd, CV, median,Q25,  Q75 ))
sample_ratio_summary1 <- inner_join(somaid,sample_ratio_summary, by= c("SomaId" ="ID" ))
write.csv(sample_ratio_summary1, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/dox_ratio_summary_09142018.csv")

###############
##PCA analysis
###############
y1 <- sample_mean1$Status
y2 <- sample_mean1$SubjectID
y3 <- sample_mean1$TimePoint
y4 <- sample_mean1$SampleId
sample_pca<-pca(log(sample_mean1[,-c(1:4)],10),ncomp=10,scale=TRUE)
pcaplot(analysis=sample_pca,y=y1, title="samplelog10_04142018.pca",ind.names=y2, cex=6)
pcaplot(analysis=sample_pca,y=y1, title="samplelog10_04142018_timepoint.pca",ind.names=y3, cex=6)
pcaplot(analysis=sample_pca,y=y1, title="samplelog10_04142018_sampleid.pca",ind.names=y4, cex=4)

sample_pcaml<-pca(log(sample_mean1[,-c(1:4)],10), ncomp=10,scale=TRUE, multilevel = data.frame(sample=y2))
pcaplot(analysis=sample_pcaml,y=y1, title="samplelog10_04142018.pcaml",ind.names=y2, cex=6)
pcaplot(analysis=sample_pcaml,y=y1, title="samplelog10_04142018_timepoint.pcaml",ind.names=y3, cex=6)
pcaplot(analysis=sample_pcaml,y=y1, title="samplelog10_04142018_sampleid.pcaml",ind.names=y4, cex=4)

y1 <- BCsample_ratio$Status
y2 <- BCsample_ratio$SubjectID
y3 <- BCsample_ratio$TimePoint
y4 <- BCsample_ratio$SampleId

analysis<-pca(log(BCsample_ratio[,-c(1:4)],10),ncomp=10,scale=TRUE)
pcaplot(analysis=analysis,y=y1, title="sample_ratio_log10_04142018.pca",ind.names=y2, cex=6)
pcaplot(analysis=analysis,y=y1, title="sample_ratio_log10_04142018_sampleid.pca",ind.names=y4, cex=6)

analysis<-plsda(log(BCsample_ratio[,-c(1:4)],10), y1, ncomp=10,scale=TRUE)
pcaplot(analysis=analysis,y=y1, title="sample_ratio_log10_04142018.plsda",ind.names=y2, cex=6)
pcaplot(analysis=analysis,y=y1, title="sample_ratio_log10_04142018_sampleid.plsda",ind.names=y4, cex=6)

sig <- sigFDRnosig_all(data=BC_ratio_re,title="BC_ratio_04142018_fold12p05", fc=1.2,p=0.05, fdr=2, savere=TRUE)

sigwilcox <- sigFDRnosig_all(data=BC_ratio_w_re,title="BC_ratio_04142018_wilcox_fold12p05", fc=1.2,p=0.05, fdr=2, savere=TRUE)


analysis<-plsda(log(BCsample_ratio[,sig$all$analytes],10), y1, ncomp=10,scale=TRUE)
pcaplot(analysis=analysis,y=y1, title="sample_ratio_log10_04142018sig170.plsda",ind.names=y2, cex=6)
pcaplot(analysis=analysis,y=y1, title="sample_ratio_log10_04142018sig170_sampleid.plsda",ind.names=y4, cex=6)


dev.off()


#t-test

sample <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/sample_replace_duplicate_with_mean_04142018_subject33labchanged.csv")
allsample_ratio <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/allsample_ratio.csv")
allsample <- sample
dim(allsample)
BCsample_ratio <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/BCsample_ratio.csv")
BAsample_ratio <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/Bsample_ratio.csv")
CAsample_ratio <- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/Csample_ratio.csv")

allsample[allsample$SubjectID==33,1:6]
#### time point without sample
Anosubject <- c(15, 44, 45)
Bnosubject <- c(5, 14, 15, 31, 44,45, 48, 53, 54, 60) #5, 14,31,48,53,54,60
Cnosubject <- c(3, 8, 12, 15, 44, 45, 48, 52, 53, 54, 60, 62) #3,8,12, 44,45, 48,52, 53, 54, 60, 62
BCnosubject <- c(15, 44, 45)
BAsample <- allsample[!(allsample$SubjectID %in%Bnosubject), ]
CAsample <- allsample[!(allsample$SubjectID %in%Cnosubject), ]
BCsample <-  allsample[!(allsample$SubjectID %in%BCnosubject), ]

set10_noA <- set10[!(set10$SampleDescription %in%Anosubject), ]
set10new239_noA <- set10new239[!(set10new239$SampleDescription %in%Anosubject), ]

ratio <- function(x){x/x[1]}

Bsample_ratio <- Bsample %>% group_by(SubjectID) %>% mutate_if(is.numeric, ratio) %>% data.frame()
Csample_ratio <- Csample %>% group_by(SubjectID) %>% mutate_if(is.numeric, ratio) %>% data.frame()
BCsample_ratio <- BCsample %>% group_by(SubjectID) %>% mutate_if(is.numeric, ratio) %>% data.frame()
#allsample_ratio <- allsample %>% group_by(SubjectID) %>% mutate_if(is.numeric, ratio) %>% data.frame()
set10_ratio <- set10_noA %>% group_by(SampleDescription) %>% mutate_if(is.numeric, ratio) %>% data.frame()
set10new239_ratio <- set10new239_noA %>% group_by(SampleDescription) %>% mutate_if(is.numeric, ratio) %>% data.frame()
write.xlsx(list(set10_ratio=set10_ratio, set10new239_ratio=set10new239_ratio), "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/08212019/set10_ratio.xlsx")

set10_re <- ttest(sample_mean1=set10,locprotein = c(31:1347))
set10_ratio_re <- ratio_ttest(sample_mean1=set10_ratio, locprotein = c(31:1347))

set10new239_re <- ttest(sample_mean1=set10new239,locprotein = c(29:1345))
set10new239_ratio_re <- ratio_ttest(sample_mean1=set10new239_ratio, locprotein = c(29:1345))

write.xlsx(list(set10_re=set10_re,set10_ratio_re=set10_ratio_re,set10new239_re=set10new239_re,set10new239_ratio_re=set10new239_ratio_re), "set10_result.xlsx")

set10_ratio_re$analytes[1:40]
dim(set10_ratio_re)
names(set10_ratio_re)
title1 <- "set10_ratio_re"
pdf("set10_ratio_re_volcanoplot.pdf", width = 8, height = 8)
#volcanoplotfdr(set10_ratio_re[26:1342, c("analytes", "A.s.t.p", "A.s.fc", "A.s.t.fdr")], title=paste(title1,"A.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10_ratio_re[26:1342, c("analytes", "B.s.t.p", "B.s.fc", "B.s.t.fdr")], title=paste(title1,"B.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10_ratio_re[26:1342, c("analytes", "C.s.t.p", "C.s.fc", "C.s.t.fdr")], title=paste(title1,"C.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)

volcanoplotfdr(set10_ratio_re[26:1342, c("analytes", "N.b.t.p", "N.b.fc", "N.b.t.fdr")], title=paste(title1,"N.b"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10_ratio_re[26:1342, c("analytes", "N.c.t.p", "N.c.fc", "N.c.t.fdr")], title=paste(title1,"N.c"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10_ratio_re[26:1342, c("analytes", "S.b.t.p", "S.b.fc", "S.b.t.fdr")], title=paste(title1,"S.b"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10_ratio_re[26:1342, c("analytes", "S.c.t.p", "S.c.fc", "S.c.t.fdr")], title=paste(title1,"S.c"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)

dev.off()

title1 <- "set10new239_ratio_re"
pdf("set10new239_ratio_re_volcanoplot.pdf", width = 8, height = 8)
#volcanoplotfdr(set10new239_ratio_re[24:1340, c("analytes", "A.s.t.p", "A.s.fc", "A.s.t.fdr")], title=paste(title1,"A.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10new239_ratio_re[24:1340, c("analytes", "B.s.t.p", "B.s.fc", "B.s.t.fdr")], title=paste(title1,"B.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10new239_ratio_re[24:1340, c("analytes", "C.s.t.p", "C.s.fc", "C.s.t.fdr")], title=paste(title1,"C.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)

volcanoplotfdr(set10new239_ratio_re[24:1340, c("analytes", "N.b.t.p", "N.b.fc", "N.b.t.fdr")], title=paste(title1,"N.b"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10new239_ratio_re[24:1340, c("analytes", "N.c.t.p", "N.c.fc", "N.c.t.fdr")], title=paste(title1,"N.c"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10new239_ratio_re[24:1340, c("analytes", "S.b.t.p", "S.b.fc", "S.b.t.fdr")], title=paste(title1,"S.b"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10new239_ratio_re[24:1340, c("analytes", "S.c.t.p", "S.c.fc", "S.c.t.fdr")], title=paste(title1,"S.c"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)

dev.off()

#intensity

title1 <- "set10_re"
pdf("set10_re_volcanoplot.pdf", width = 8, height = 8)
volcanoplotfdr(set10_re[26:1342, c("analytes", "A.s.t.p", "A.s.fc", "A.s.t.fdr")], title=paste(title1,"A.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10_re[26:1342, c("analytes", "B.s.t.p", "B.s.fc", "B.s.t.fdr")], title=paste(title1,"B.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10_re[26:1342, c("analytes", "C.s.t.p", "C.s.fc", "C.s.t.fdr")], title=paste(title1,"C.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)

volcanoplotfdr(set10_re[26:1342, c("analytes", "N.b.t.p", "N.b.fc", "N.b.t.fdr")], title=paste(title1,"N.b"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10_re[26:1342, c("analytes", "N.c.t.p", "N.c.fc", "N.c.t.fdr")], title=paste(title1,"N.c"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10_re[26:1342, c("analytes", "S.b.t.p", "S.b.fc", "S.b.t.fdr")], title=paste(title1,"S.b"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10_re[26:1342, c("analytes", "S.c.t.p", "S.c.fc", "S.c.t.fdr")], title=paste(title1,"S.c"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)

dev.off()

title1 <- "set10new239_re"
pdf("set10new239_re_volcanoplot.pdf", width = 8, height = 8)
volcanoplotfdr(set10new239_re[24:1340, c("analytes", "A.s.t.p", "A.s.fc", "A.s.t.fdr")], title=paste(title1,"A.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10new239_re[24:1340, c("analytes", "B.s.t.p", "B.s.fc", "B.s.t.fdr")], title=paste(title1,"B.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10new239_re[24:1340, c("analytes", "C.s.t.p", "C.s.fc", "C.s.t.fdr")], title=paste(title1,"C.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)

volcanoplotfdr(set10new239_re[24:1340, c("analytes", "N.b.t.p", "N.b.fc", "N.b.t.fdr")], title=paste(title1,"N.b"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10new239_re[24:1340, c("analytes", "N.c.t.p", "N.c.fc", "N.c.t.fdr")], title=paste(title1,"N.c"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10new239_re[24:1340, c("analytes", "S.b.t.p", "S.b.fc", "S.b.t.fdr")], title=paste(title1,"S.b"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(set10new239_re[24:1340, c("analytes", "S.c.t.p", "S.c.fc", "S.c.t.fdr")], title=paste(title1,"S.c"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)

dev.off()




dim(set10)
names(set10)[1:40]
names(allsample)[1:10]

dim(set10new239)
names(set10new239)[1:40]


allsample[,c("SubjectID","Status")]

str(Csample)
write.csv(Bsample_ratio, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/Bsample_ratio.csv")
write.csv(Csample_ratio, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/Csample_ratio.csv")
write.csv(BCsample_ratio, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/BCsample_ratio.csv")
write.csv(allsample_ratio, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/allsample_ratio.csv")

Bsample_ratio[1:175, 1:6]
Bsample_ratio$TimePoint
Bsample_ratio$SubjectID
names(Csample)[1:10]
dim(allsample)
dim(Bsample)
dim(Csample)
dim(BCsample)


note2 <- c(6, 26, 32, 60, 62, 64, 67, 72, 75)
note3 <- c(22, 29, 55)
note4 <- c(5, 46)
dim(allsample)
names(allsample)
note2_sample <- allsample %>% filter(!(SubjectID %in% c(note3,note4)))

note3_sample <- allsample %>% filter(!(SubjectID %in% c(note2,note4)))

note4_sample <- allsample %>% filter(!(SubjectID %in% c(note2,note3)))
str(somaid)

####boxplot: all samples
i <- 1
name <-  "SL017128"
names(allsample)[1:10]
proteins <- names(allsample)[-c(1:5)]
length(proteins)
dim(allsample)
allsample$Status <- factor(allsample$Status, levels=c("Normal", "Abnormal"))

pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/dox_log2boxplot_09142018a.pdf", width = 6,height = 5)
for (name in proteins){
 p1 <- ggplot(allsample, aes(TimePoint, log(get(name), 2),fill=Status))+
   theme_bw()+
   scale_fill_manual(values=c("chartreuse", "indianred1"))+
  geom_boxplot(position=position_dodge(0.9), outlier.alpha = 0)+
  geom_jitter(alpha=.8,size=1, position=position_jitterdodge(
    jitter.width=0.7,jitter.height=0,dodge.width =0.9))+
  labs(y=paste("log2",somaid$Target[somaid$SomaId==paste(name)]))+
   theme(axis.title=element_text(size=12),axis.text=element_text(size=8, color="black"))

print (p1)
cat(i," ")
i <- i+1
}
dev.off() 

note2_sample_ratio <- BCsample_ratio %>% filter(!(SubjectID %in% c(note3,note4)))
note3_sample_ratio <- BCsample_ratio %>% filter(!(SubjectID %in% c(note2,note4)))
note4_sample_ratio <- BCsample_ratio %>% filter(!(SubjectID %in% c(note2,note3)))

dim(allsample)
dim(note2_sample)
dim(note3_sample)
dim(note4_sample)

dim(BCsample_ratio)
dim(note2_sample_ratio)
dim(note3_sample_ratio)
dim(note4_sample_ratio)

allsample$SubjectID %in% c(note2,note3)
str(allsample)[1:9]
note4_sample[,1:5]

allsample_re <- ttest(sample_mean1=allsample)
write.csv(allsample_re, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/doxfoldp_ttest_allsample_06212018a.csv")

note2_re <- ttest(sample_mean1=note2_sample)
write.csv(note2_re, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/doxfoldp_ttest_note2_06212018.csv")

note3_re <- ttest(sample_mean1=note3_sample)
write.csv(note3_re, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/doxfoldp_ttest_note3_06212018.csv")

note4_re <- ttest(sample_mean1=note4_sample)
write.csv(note4_re, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/doxfoldp_ttest_note4_06212018.csv")

####

BCsample_ratio_re <- ratio_ttest(sample_mean1=BCsample_ratio)
write.csv(BCsample_ratio_re, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/doxfoldp_ttest_BCsample_ratio_06212018.csv")

note2_ratio_re <- ratio_ttest(sample_mean1=note2_sample_ratio)
write.csv(note2_ratio_re, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/doxfoldp_ttest_note2_ratio_06212018.csv")

note3_ratio_re <- ratio_ttest(sample_mean1=note3_sample_ratio)
write.csv(note3_ratio_re, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/doxfoldp_ttest_note3_ratio_06212018.csv")

note4_ratio_re <- ratio_ttest(sample_mean1=note4_sample_ratio)
write.csv(note4_ratio_re, "C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/06212018/doxfoldp_ttest_note4_ratio_06212018.csv")

###paired t and wilcox test
BA_paired_re <- BA_ttest(sample_mean1=BAsample)
BA_ratio_paired_re <- BA_ttest(sample_mean1=BAsample_ratio)

CA_paired_re <- CA_ttest(sample_mean1=CAsample)
CA_ratio_paired_re <- CA_ttest(sample_mean1=CAsample_ratio)

write_excel_paired(BA_re=BA_paired_re, CA_re=CA_paired_re, savename = "paired_p05fc1d2_09142018.xlsx", pvalue = 0.05, fc = 1.2)
write_excel_paired(BA_re=BA_ratio_paired_re, CA_re=CA_ratio_paired_re, savename = "ratio_paired_p05fc1d2_09142018.xlsx", pvalue = 0.05, fc = 1.2)

library(openxlsx)






pdf(paste(title1,"_volcanoplot.pdf", sep=""), width = 8, height = 8)
#volcanoplotfdr(sample_re[, c("Target", "A.s.t.p", "A.s.fc", "A.s.t.fdr")], title=paste(title1,"A.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
#volcanoplotfdr(sample_re[, c("Target", "B.s.t.p", "B.s.fc", "B.s.t.fdr")], title=paste(title1,"B.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(sample_re[, c("Target", "C.s.t.p", "C.s.fc", "C.s.t.fdr")], title=paste(title1,"C.s"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)

volcanoplotfdr(sample_re[, c("Target", "N.b.t.p", "N.b.fc", "N.b.t.fdr")], title=paste(title1,"N.b"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(sample_re[, c("Target", "N.c.t.p", "N.c.fc", "N.c.t.fdr")], title=paste(title1,"N.c"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
#volcanoplotfdr(sample_re[, c("Target", "S.b.t.p", "S.b.fc", "S.b.t.fdr")], title=paste(title1,"S.b"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)
volcanoplotfdr(sample_re[, c("Target", "S.c.t.p", "S.c.fc", "S.c.t.fdr")], title=paste(title1,"S.c"),fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=2)

dev.off()




#############################
####t test ratio
#############################
ratio_ttest <- function(sample_mean1=BCsample_ratio, locprotein){
  
  resultfoldp<-data.frame("analytes"=character(),
                          "B.s.t.p"=numeric(),"B.s.w.p"=numeric(),"B.s.fc"=numeric(),
                          "C.s.t.p"=numeric(),"C.s.w.p"=numeric(),"C.s.fc"=numeric(),
                          "N.b.t.p"=numeric(),"N.b.w.p"=numeric(),"N.b.fc"=numeric(),"N.c.t.p"=numeric(),"N.c.w.p"=numeric(),"N.c.fc"=numeric(),
                          "S.b.t.p"=numeric(),"S.b.w.p"=numeric(),"S.b.fc"=numeric(),"S.c.t.p"=numeric(),"S.c.w.p"=numeric(),"S.c.fc"=numeric(),stringsAsFactors=FALSE)
  for (i in locprotein){
    #i=6
    B.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="B",i]
    
    B.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="B",i]
    
    B.s.fc<-mean(B.s)/mean(B.n)
    
    if (length(B.s)>=2){
      B.s.t<-t.test(log(B.s,10), log(B.n,10),paired=FALSE,var.qual=FALSE)$p.value
      B.s.w<-wilcox.test(log(B.s,10), log(B.n,10),paired=FALSE,var.qual=FALSE)$p.value
    }    else {B.s.t=NA
               B.s.w=NA}
    
    C.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="C",i]
    
    C.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="C",i]
    
    
    C.s.fc<-mean(C.s)/mean(C.n)
    
    C.s.t<-t.test(log(C.s,10), log(C.n,10),paired=FALSE,var.qual=FALSE)$p.value
    C.s.w<-wilcox.test(log(C.s,10), log(C.n,10),paired=FALSE,var.qual=FALSE)$p.value
    
    N.b.fc<-mean(B.n)
    N.c.fc<-mean(C.n)
    N.b.t<-t.test(log(B.n,10), mu=0, paired=FALSE,var.qual=FALSE)$p.value
    N.b.w<-wilcox.test(log(B.n,10), mu=0, paired=FALSE,var.qual=FALSE)$p.value
    N.c.t<-t.test(log(C.n,10), mu=0, paired=FALSE,var.qual=FALSE)$p.value
    N.c.w<-wilcox.test(log(C.n,10), mu=0,paired=FALSE,var.qual=FALSE)$p.value
    
    S.b.fc<-mean(B.s)
    S.c.fc<-mean(C.s)
    
    if(length(B.s)>=2){
      S.b.t<-t.test(log(B.s,10), mu=0,paired=FALSE,var.qual=FALSE)$p.value
      S.b.w<-wilcox.test(log(B.s,10), mu=0,paired=FALSE,var.qual=FALSE)$p.value
    }      else {S.b.t=NA
                 S.b.w=NA}
    S.c.t<-t.test(log(C.s,10), mu=0,paired=FALSE,var.qual=FALSE)$p.value
    S.c.w<-wilcox.test(log(C.s,10), mu=0,paired=FALSE,var.qual=FALSE)$p.value
    
    resultfoldp[i-5,]<-c(names(sample_mean1)[i],
                         B.s.t,B.s.w, B.s.fc,
                         C.s.t, C.s.w, C.s.fc,
                         N.b.t, N.b.w, N.b.fc, N.c.t,N.c.w, N.c.fc,
                         S.b.t, S.b.w, S.b.fc, S.c.t,S.c.w, S.c.fc)
    
    
  }
  
  resultfoldp[,-1] <- sapply(resultfoldp[,-1], as.numeric)
  
  resultfoldp[,"B.s.t.fdr"]<-p.adjust(resultfoldp$B.s.t.p,"BH")
  
  resultfoldp[,"C.s.t.fdr"]<-p.adjust(resultfoldp$C.s.t.p,"BH")
  
  resultfoldp[,"N.b.t.fdr"]<-p.adjust(resultfoldp$N.b.t.p,"BH")
  resultfoldp[,"N.c.t.fdr"]<-p.adjust(resultfoldp$N.c.t.p,"BH")
  
  resultfoldp[,"S.b.t.fdr"]<-p.adjust(resultfoldp$S.b.t.p,"BH")
  resultfoldp[,"S.c.t.fdr"]<-p.adjust(resultfoldp$S.c.t.p,"BH")
  
  #resultfoldp1 <- inner_join(somaid,resultfoldp, by= c("SomaId" ="analytes" ))
  return (resultfoldp)
}



###########################
###intensity no pair
###########################
ttest <- function(sample_mean1=allsample, locprotein){
  
  resultfoldp<-data.frame("analytes"=character(),"A.s.t.p"=numeric(),"A.s.w.p"=numeric(),"A.s.fc"=numeric(),
                          "B.s.t.p"=numeric(),"B.s.w.p"=numeric(),"B.s.fc"=numeric(),
                          "C.s.t.p"=numeric(),"C.s.w.p"=numeric(),"C.s.fc"=numeric(),
                          "N.b.t.p"=numeric(),"N.b.w.p"=numeric(),"N.b.fc"=numeric(),"N.c.t.p"=numeric(),"N.c.w.p"=numeric(),"N.c.fc"=numeric(),
                          "S.b.t.p"=numeric(),"S.b.w.p"=numeric(),"S.b.fc"=numeric(),"S.c.t.p"=numeric(),"S.c.w.p"=numeric(),"S.c.fc"=numeric(),stringsAsFactors=FALSE)
  for (i in locprotein){
    #i=6
    A.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="A",i]
    
    A.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="A",i]
    
    
    A.s.fc<-mean(A.s)/mean(A.n)
    
    A.s.t<-t.test(log(A.s,10), log(A.n,10),paired=FALSE,var.qual=FALSE)$p.value
    A.s.w<-wilcox.test(log(A.s,10), log(A.n,10),paired=FALSE,var.qual=FALSE)$p.value
    
    B.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="B",i]
    
    B.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="B",i]
    
    B.s.fc<-mean(B.s)/mean(B.n)
    
    if (length(B.s)>=2){
    B.s.t<-t.test(log(B.s,10), log(B.n,10),paired=FALSE,var.qual=FALSE)$p.value
    B.s.w<-wilcox.test(log(B.s,10), log(B.n,10),paired=FALSE,var.qual=FALSE)$p.value
    }    else {B.s.t=NA
               B.s.w=NA}
    
    C.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="C",i]
    
    C.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="C",i]
    
    
    C.s.fc<-mean(C.s)/mean(C.n)
    
    C.s.t<-t.test(log(C.s,10), log(C.n,10),paired=FALSE,var.qual=FALSE)$p.value
    C.s.w<-wilcox.test(log(C.s,10), log(C.n,10),paired=FALSE,var.qual=FALSE)$p.value


      
      N.b.fc<-mean(B.n)/mean(A.n)
      N.c.fc<-mean(C.n)/mean(A.n)
      N.b.t<-t.test(log(B.n,10), log(A.n,10),paired=FALSE,var.qual=FALSE)$p.value
      N.b.w<-wilcox.test(log(B.n,10), log(A.n,10),paired=FALSE,var.qual=FALSE)$p.value
      N.c.t<-t.test(log(C.n,10), log(A.n,10),paired=FALSE,var.qual=FALSE)$p.value
      N.c.w<-wilcox.test(log(C.n,10), log(A.n,10),paired=FALSE,var.qual=FALSE)$p.value
      
      S.b.fc<-mean(B.s)/mean(A.s)
      S.c.fc<-mean(C.s)/mean(A.s)
      
      if(length(B.s)>=2){
      S.b.t<-t.test(log(B.s,10), log(A.s,10),paired=FALSE,var.qual=FALSE)$p.value
      S.b.w<-wilcox.test(log(B.s,10), log(A.s,10),paired=FALSE,var.qual=FALSE)$p.value
      }      else {S.b.t=NA
                   S.b.w=NA}
      S.c.t<-t.test(log(C.s,10), log(A.s,10),paired=FALSE,var.qual=FALSE)$p.value
      S.c.w<-wilcox.test(log(C.s,10), log(A.s,10),paired=FALSE,var.qual=FALSE)$p.value

    
    
    resultfoldp[i-5,]<-c(names(sample_mean1)[i],A.s.t,A.s.w, A.s.fc,
                         B.s.t,B.s.w, B.s.fc,
                         C.s.t, C.s.w, C.s.fc,
                         N.b.t, N.b.w, N.b.fc, N.c.t,N.c.w, N.c.fc,
                         S.b.t, S.b.w, S.b.fc, S.c.t,S.c.w, S.c.fc)

    
  }
  
  resultfoldp[,-1] <- sapply(resultfoldp[,-1], as.numeric)
 
  
  resultfoldp[,"A.s.t.fdr"]<-p.adjust(resultfoldp$A.s.t.p,"BH")
  
  resultfoldp[,"B.s.t.fdr"]<-p.adjust(resultfoldp$B.s.t.p,"BH")
  
  resultfoldp[,"C.s.t.fdr"]<-p.adjust(resultfoldp$C.s.t.p,"BH")
  
  resultfoldp[,"N.b.t.fdr"]<-p.adjust(resultfoldp$N.b.t.p,"BH")
  resultfoldp[,"N.c.t.fdr"]<-p.adjust(resultfoldp$N.c.t.p,"BH")
  
  
  resultfoldp[,"S.b.t.fdr"]<-p.adjust(resultfoldp$S.b.t.p,"BH")
  resultfoldp[,"S.c.t.fdr"]<-p.adjust(resultfoldp$S.c.t.p,"BH")
  
  #resultfoldp1 <- inner_join(somaid,resultfoldp, by= c("SomaId" ="analytes" ))
  return (resultfoldp)
  
}

#######pair test#################3
BA_ttest <- function(sample_mean1=Bsample){
  
  resultfoldp<-data.frame("analytes"=character(),
                          "N.b.t.p"=numeric(),"N.b.w.p"=numeric(),"N.b.fc"=numeric(),
                          "S.b.t.p"=numeric(),"S.b.w.p"=numeric(),"S.b.fc"=numeric(),stringsAsFactors=FALSE)
  for (i in 6:1310){
    #i=6
    A.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="A",i]
    A.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="A",i]
    B.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="B",i]
    B.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="B",i]
    
    
    N.b.fc<-mean(B.n)/mean(A.n)
    N.b.t<-t.test(log(B.n,10), log(A.n,10),paired=TRUE,var.qual=FALSE)$p.value
    N.b.w<-wilcox.test(log(B.n,10), log(A.n,10),paired=TRUE,var.qual=FALSE)$p.value
 
    S.b.fc<-mean(B.s)/mean(A.s)
  
    if(length(B.s)>=2){
      S.b.t<-t.test(log(B.s,10), log(A.s,10),paired=TRUE,var.qual=FALSE)$p.value
      S.b.w<-wilcox.test(log(B.s,10), log(A.s,10),paired=TRUE,var.qual=FALSE)$p.value
    }      else {S.b.t=NA
                 S.b.w=NA}
  
    resultfoldp[i-5,]<-c(names(sample_mean1)[i],
                         N.b.t, N.b.w, N.b.fc, 
                         S.b.t, S.b.w, S.b.fc)
    
  }
  
  resultfoldp[,-1] <- sapply(resultfoldp[,-1], as.numeric)
  resultfoldp[,"N.b.t.fdr"]<-p.adjust(resultfoldp$N.b.t.p,"BH")
  resultfoldp[,"S.b.t.fdr"]<-p.adjust(resultfoldp$S.b.t.p,"BH")
  resultfoldp1 <- inner_join(somaid,resultfoldp, by= c("SomaId" ="analytes" ))
  
}

##########

CA_ttest <- function(sample_mean1=CAsample){
  
  resultfoldp<-data.frame("analytes"=character(),
                          "N.c.t.p"=numeric(),"N.c.w.p"=numeric(),"N.c.fc"=numeric(),
                          "S.c.t.p"=numeric(),"S.c.w.p"=numeric(),"S.c.fc"=numeric(),stringsAsFactors=FALSE)
  for (i in 6:1310){
    #i=6
    A.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="A",i]
    A.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="A",i]
    C.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="C",i]
    C.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="C",i]
    
    N.c.fc<-mean(C.n)/mean(A.n)
    N.c.t<-t.test(log(C.n,10), log(A.n,10),paired=TRUE,var.qual=FALSE)$p.value
    N.c.w<-wilcox.test(log(C.n,10), log(A.n,10),paired=TRUE,var.qual=FALSE)$p.value
    
    S.c.fc<-mean(C.s)/mean(A.s)
    S.c.t<-t.test(log(C.s,10), log(A.s,10),paired=TRUE,var.qual=FALSE)$p.value
    S.c.w<-wilcox.test(log(C.s,10), log(A.s,10),paired=TRUE,var.qual=FALSE)$p.value
    
    resultfoldp[i-5,]<-c(names(sample_mean1)[i],
                          N.c.t,N.c.w, N.c.fc,
                          S.c.t,S.c.w, S.c.fc)
    }
  
  resultfoldp[,-1] <- sapply(resultfoldp[,-1], as.numeric)
  resultfoldp[,"N.c.t.fdr"]<-p.adjust(resultfoldp$N.c.t.p,"BH")
  resultfoldp[,"S.c.t.fdr"]<-p.adjust(resultfoldp$S.c.t.p,"BH")
  
  resultfoldp1 <- inner_join(somaid,resultfoldp, by= c("SomaId" ="analytes" ))
  
}

###ratio paired test

BA_ratio_ttest <- function(sample_mean1=Bsample){
  
  resultfoldp<-data.frame("analytes"=character(),
                          "N.b.t.p"=numeric(),"N.b.w.p"=numeric(),"N.b.fc"=numeric(),
                          "S.b.t.p"=numeric(),"S.b.w.p"=numeric(),"S.b.fc"=numeric(),stringsAsFactors=FALSE)
  for (i in 6:1310){
    #i=6
    A.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="A",i]
    A.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="A",i]
    B.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="B",i]
    B.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="B",i]
    
    
    N.b.fc<-mean(B.n)/mean(A.n)
    N.b.t<-t.test(log(B.n,10), mu=0, paired=TRUE,var.qual=FALSE)$p.value
    N.b.w<-wilcox.test(log(B.n,10), mu=0, paired=TRUE, var.qual=FALSE)$p.value
    
    S.b.fc<-mean(B.s)/mean(A.s)
    
    if(length(B.s)>=2){
      S.b.t<-t.test(log(B.s,10), mu=0, paired=TRUE,var.qual=FALSE)$p.value
      S.b.w<-wilcox.test(log(B.s,10), mu=0, paired=TRUE,var.qual=FALSE)$p.value
    }      else {S.b.t=NA
    S.b.w=NA}
    
    resultfoldp[i-5,]<-c(names(sample_mean1)[i],
                         N.b.t, N.b.w, N.b.fc, 
                         S.b.t, S.b.w, S.b.fc)
    
  }
  
  resultfoldp[,-1] <- sapply(resultfoldp[,-1], as.numeric)
  resultfoldp[,"N.b.t.fdr"]<-p.adjust(resultfoldp$N.b.t.p,"BH")
  resultfoldp[,"S.b.t.fdr"]<-p.adjust(resultfoldp$S.b.t.p,"BH")
  resultfoldp1 <- inner_join(somaid,resultfoldp, by= c("SomaId" ="analytes" ))
  
}

##########

CA_ratio_ttest <- function(sample_mean1=CAsample){
  
  resultfoldp<-data.frame("analytes"=character(),
                          "N.c.t.p"=numeric(),"N.c.w.p"=numeric(),"N.c.fc"=numeric(),
                          "S.c.t.p"=numeric(),"S.c.w.p"=numeric(),"S.c.fc"=numeric(),stringsAsFactors=FALSE)
  for (i in 6:1310){
    #i=6
    A.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="A",i]
    A.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="A",i]
    C.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="C",i]
    C.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="C",i]
    
    N.c.fc<-mean(C.n)/mean(A.n)
    N.c.t<-t.test(log(C.n,10), mu=0, paired=TRUE,var.qual=FALSE)$p.value
    N.c.w<-wilcox.test(log(C.n,10),mu=0, paired=TRUE,var.qual=FALSE)$p.value
    
    S.c.fc<-mean(C.s)/mean(A.s)
    S.c.t<-t.test(log(C.s,10), mu=0, paired=TRUE,var.qual=FALSE)$p.value
    S.c.w<-wilcox.test(log(C.s,10), mu=0, paired=TRUE,var.qual=FALSE)$p.value
    
    resultfoldp[i-5,]<-c(names(sample_mean1)[i],
                         N.c.t,N.c.w, N.c.fc,
                         S.c.t,S.c.w, S.c.fc)
  }
  
  resultfoldp[,-1] <- sapply(resultfoldp[,-1], as.numeric)
  resultfoldp[,"N.c.t.fdr"]<-p.adjust(resultfoldp$N.c.t.p,"BH")
  resultfoldp[,"S.c.t.fdr"]<-p.adjust(resultfoldp$S.c.t.p,"BH")
  
  resultfoldp1 <- inner_join(somaid,resultfoldp, by= c("SomaId" ="analytes" ))
  
}



















###########################
###intensity paired
###########################
ttest <- function(sample_mean1=Bsample, BC="B"){

resultfoldp<-data.frame("analytes"=character(),"A.s.t.p"=numeric(),"A.s.fc"=numeric(),
                        "B.s.t.p"=numeric(),"B.s.fc"=numeric(),
                        "C.s.t.p"=numeric(),"C.s.fc"=numeric(),
                        "N.b.t.p"=numeric(),"N.b.fc"=numeric(),"N.c.t.p"=numeric(),"N.c.fc"=numeric(),
                        "S.b.t.p"=numeric(),"S.b.fc"=numeric(),"S.c.t.p"=numeric(),"S.c.fc"=numeric(),stringsAsFactors=FALSE)
for (i in 5:1039){

  A.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="A",i]

  A.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="A",i]
  
 
  A.s.fc<-mean(A.s)/mean(A.n)
 
  A.s.t<-t.test(log(A.s,10), log(A.n,10),paired=FALSE,var.qual=FALSE)
  
  B.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="B",i]

  B.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="B",i]
  

  B.s.fc<-mean(B.s)/mean(B.n)

  B.s.t<-t.test(log(B.s,10), log(B.n,10),paired=FALSE,var.qual=FALSE)
  
  C.n<-sample_mean1[sample_mean1$Status=="Normal"&sample_mean1$TimePoint=="C",i]
 
  C.s<-sample_mean1[sample_mean1$Status=="Abnormal"&sample_mean1$TimePoint=="C",i]
  

  C.s.fc<-mean(C.s)/mean(C.n)
 
  C.s.t<-t.test(log(C.s,10), log(C.n,10),paired=FALSE,var.qual=FALSE)
  
  if (BC=="B") {
 
   N.b.fc<-mean(B.n)/mean(A.n)
   N.c.fc<-mean(C.n)/mean(A.n)
   N.b.t<-t.test(log(B.n,10), log(A.n,10),paired=TRUE,var.qual=FALSE)
   N.c.t<-t.test(log(C.n,10), log(A.n,10),paired=FALSE,var.qual=FALSE)
  

   S.b.fc<-mean(B.s)/mean(A.s)
   S.c.fc<-mean(C.s)/mean(A.s)
   S.b.t<-t.test(log(B.s,10), log(A.s,10),paired=TRUE,var.qual=FALSE)
   S.c.t<-t.test(log(C.s,10), log(A.s,10),paired=FALSE,var.qual=FALSE)
   
  }
  
  
  if (BC=="C") {
    
    N.b.fc<-mean(B.n)/mean(A.n)
    N.c.fc<-mean(C.n)/mean(A.n)
    N.b.t<-t.test(log(B.n,10), log(A.n,10),paired=FALSE,var.qual=FALSE)
    N.c.t<-t.test(log(C.n,10), log(A.n,10),paired=TRUE,var.qual=FALSE)
    
    
    S.b.fc<-mean(B.s)/mean(A.s)
    S.c.fc<-mean(C.s)/mean(A.s)
    S.b.t<-t.test(log(B.s,10), log(A.s,10),paired=FALSE,var.qual=FALSE)
    S.c.t<-t.test(log(C.s,10), log(A.s,10),paired=TRUE,var.qual=FALSE)
    
  }
  
  
  resultfoldp[i-4,]<-c(names(sample_mean1)[i],A.s.t$p.value,A.s.fc,
                     B.s.t$p.value,B.s.fc,
                     C.s.t$p.value,C.s.fc,
                     N.b.t$p.value,N.b.fc,N.c.t$p.value,N.c.fc,
                     S.b.t$p.value,S.b.fc,S.c.t$p.value,S.c.fc)
  
  
}

resultfoldp[,-1] <- sapply(resultfoldp[,-1], as.numeric)
# str(resultfoldp)
# class(resultfoldp)
# 
# names(resultfoldp)
# head(resultfoldp,20)
# str(resultfoldp)

resultfoldp[,"A.s.t.fdr"]<-p.adjust(resultfoldp$A.s.t.p,"BH")

resultfoldp[,"B.s.t.fdr"]<-p.adjust(resultfoldp$B.s.t.p,"BH")

resultfoldp[,"C.s.t.fdr"]<-p.adjust(resultfoldp$C.s.t.p,"BH")

resultfoldp[,"N.b.t.fdr"]<-p.adjust(resultfoldp$N.b.t.p,"BH")
resultfoldp[,"N.c.t.fdr"]<-p.adjust(resultfoldp$N.c.t.p,"BH")


resultfoldp[,"S.b.t.fdr"]<-p.adjust(resultfoldp$S.b.t.p,"BH")
resultfoldp[,"S.c.t.fdr"]<-p.adjust(resultfoldp$S.c.t.p,"BH")

resultfoldp1 <- left_join(resultfoldp, somaid,by= c("analytes" = "SomaId"))

}



#vocanoplot
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/soma_dox_volcanoplot_fc1d2_04142018_w7h7.pdf", width=7, height=7)

volcanoplot(data.frame(re$Target,re$A.s.t.p,re$A.s.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="A: Abnormal vs Normal")

volcanoplot(data.frame(re$Target,re$B.s.t.p,re$B.s.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="B: Abnormal vs Normal")

volcanoplot(data.frame(re$Target,re$C.s.t.p,re$C.s.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="C: Abnormal vs Normal")

volcanoplot(data.frame(re$Target,re$N.b.t.p,re$N.b.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Normal: B vs A")

volcanoplot(data.frame(re$Target,re$N.c.t.p,re$N.c.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Normal C vs A")

volcanoplot(data.frame(re$Target,re$S.b.t.p,re$S.b.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Abnormal: B vs A")

volcanoplot(data.frame(re$Target,re$S.c.t.p,re$S.c.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Abnormal C vs A")

dev.off()


pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/soma_dox_volcanoplot_paired_fc1d2_04142018_w7h7a.pdf", width=7, height=7)

volcanoplot(data.frame(B_re$Target, B_re$N.b.t.p, B_re$N.b.fc), fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Normal: B vs A paired")

volcanoplot(data.frame(C_re$Target, C_re$N.c.t.p, C_re$N.c.fc), fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Normal C vs A paired")

volcanoplot(data.frame(B_re$Target, B_re$S.b.t.p, B_re$S.b.fc), fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Abnormal: B vs A paired")

volcanoplot(data.frame(C_re$Target, C_re$S.c.t.p, C_re$S.c.fc), fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Abnormal C vs A paired")

dev.off()


pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/soma_dox_volcanoplot_paired_fc1d2_04142018_w7h7a_wilcox.pdf", width=7, height=7)

volcanoplot(data.frame(B_w_re$Target, B_w_re$N.b.t.p, B_w_re$N.b.fc), fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Normal: B vs A paired")

volcanoplot(data.frame(C_w_re$Target, C_w_re$N.c.t.p, C_w_re$N.c.fc), fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Normal C vs A paired")

volcanoplot(data.frame(B_w_re$Target, B_w_re$S.b.t.p, B_w_re$S.b.fc), fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Abnormal: B vs A paired")

volcanoplot(data.frame(C_w_re$Target, C_w_re$S.c.t.p, C_w_re$S.c.fc), fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Abnormal C vs A paired")

dev.off()


pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/soma_dox_volcanoplot_BC_ratio_fc1d2_04142018_w7h7.pdf", width=7, height=7)


volcanoplot(data.frame(BC_ratio_re$Target,BC_ratio_re$B.s.t.p,BC_ratio_re$B.s.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="B: Abnormal vs Normal")

volcanoplot(data.frame(BC_ratio_re$Target,BC_ratio_re$C.s.t.p,BC_ratio_re$C.s.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="C: Abnormal vs Normal")

volcanoplot(data.frame(BC_ratio_re$Target,BC_ratio_re$N.b.t.p,BC_ratio_re$N.b.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Normal: B vs A")

volcanoplot(data.frame(BC_ratio_re$Target,BC_ratio_re$N.c.t.p,BC_ratio_re$N.c.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Normal C vs A")

volcanoplot(data.frame(BC_ratio_re$Target,BC_ratio_re$S.b.t.p,BC_ratio_re$S.b.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Abnormal: B vs A")

volcanoplot(data.frame(BC_ratio_re$Target,BC_ratio_re$S.c.t.p,BC_ratio_re$S.c.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Abnormal C vs A")

dev.off()


pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/04142018/soma_dox_volcanoplot_BC_ratio_fc1d2_04142018_w7h7_wilcox.pdf", width=7, height=7)


volcanoplot(data.frame(BC_ratio_w_re$Target,BC_ratio_w_re$B.s.t.p,BC_ratio_w_re$B.s.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="B: Abnormal vs Normal")

volcanoplot(data.frame(BC_ratio_w_re$Target,BC_ratio_w_re$C.s.t.p,BC_ratio_w_re$C.s.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="C: Abnormal vs Normal")

volcanoplot(data.frame(BC_ratio_w_re$Target,BC_ratio_w_re$N.b.t.p,BC_ratio_w_re$N.b.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Normal: B vs A")

volcanoplot(data.frame(BC_ratio_w_re$Target,BC_ratio_w_re$N.c.t.p,BC_ratio_w_re$N.c.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Normal C vs A")

volcanoplot(data.frame(BC_ratio_w_re$Target,BC_ratio_w_re$S.b.t.p,BC_ratio_w_re$S.b.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Abnormal: B vs A")

volcanoplot(data.frame(BC_ratio_w_re$Target,BC_ratio_w_re$S.c.t.p,BC_ratio_w_re$S.c.fc),fccutoff=1.2, pcutoff=0.05, label="labelsig", title="Abnormal C vs A")

dev.off()



######################################################
##########extract significantly changed data
######################################################

sigFDRnosig_all<-function(data=BC_ratio_re,title="tem", fc=1.2,p=0.05, fdr=2, savere=TRUE){
  A<-data%>%filter((B.s.fc>=fc|B.s.fc<=1/fc)& B.s.t.p<p & B.s.t.fdr<fdr)%>%droplevels()
  B<-data%>%filter((C.s.fc>=fc|C.s.fc<=1/fc)& C.s.t.p<p & C.s.t.fdr<fdr)%>%droplevels()
  
  C<-data%>%filter((N.b.fc>=fc|N.b.fc<=1/fc)& N.b.t.p<p & N.b.t.fdr<fdr)%>%droplevels()
  D<-data%>%filter((N.c.fc>=fc|N.c.fc<=1/fc)& N.c.t.p<p & N.c.t.fdr<fdr)%>%droplevels()
  
  E<-data%>%filter((S.b.fc>=fc|S.b.fc<=1/fc)& S.b.t.p<p & S.b.t.fdr<fdr)%>%droplevels()
  FF<-data%>%filter((S.c.fc>=fc|S.c.fc<=1/fc)& S.c.t.p<p & S.c.t.fdr<fdr)%>%droplevels()
  
  data[,"A"] <- ifelse(data$analyte %in% A$analyte, 1,0)
  data[,"B"] <- ifelse(data$analyte %in% B$analyte, 1,0)
  data[,"C"] <- ifelse(data$analyte %in% C$analyte, 1,0)
  data[,"D"] <- ifelse(data$analyte %in% D$analyte, 1,0)
  data[,"E"] <- ifelse(data$analyte %in% E$analyte, 1,0)
  data[,"FF"] <- ifelse(data$analyte %in% FF$analyte, 1,0)
  
  data[,"ABCDEF"] <- paste(data$A, data$B, data$C, data$D, data$E, data$FF, sep="")
  
  nosig <- list(B.s=data$analyte[data$A==0],C.s=data$analyte[data$B==0], N.b=data$analyte[data$C==0],
                N.c=data$analyte[data$D==0],S.b=data$analyte[data$E==0], S.c=data$analyte[data$FF==0])

  sig <- list(B.s=data$analyte[data$A==1],C.s=data$analyte[data$B==1], N.b=data$analyte[data$C==1],
                N.c=data$analyte[data$D==1],S.b=data$analyte[data$E==1], S.c=data$analyte[data$FF==1])
  
  abcdef <- unique(c(paste(A$analyte),paste(B$analyte),paste(C$analyte),paste(D$analyte),paste(E$analyte), paste(FF$analyte)))
  all<- data[data$analyte%in%abcdef, ]
  all[,"A"] <- ifelse(all$analyte %in% A$analyte, 1,0)
  all[,"B"] <- ifelse(all$analyte %in% B$analyte, 1,0)
  all[,"C"] <- ifelse(all$analyte %in% C$analyte, 1,0)
  all[,"D"] <- ifelse(all$analyte %in% D$analyte, 1,0)
  all[,"E"] <- ifelse(all$analyte %in% E$analyte, 1,0)
  all[,"FF"] <- ifelse(all$analyte %in% FF$analyte, 1,0)
  all[,"ABCDEF"] <- paste(all$A, all$B, all$C, all$D, all$E, all$FF, sep="")
  
  #paste(c(1,1,1,1),c(0,0,0,0),c(1,1,1,1),sep="")
  
  if (savere){
    write.xlsx2(all, paste(title, ".xlsx"), sheetName="all",append=TRUE )
    
    write.xlsx2(all[all$analyte%in%A$analyte, ], paste(title, ".xlsx"), sheetName="B.s",append=TRUE )
    write.xlsx2(all[all$analyte%in%B$analyte, ], paste(title, ".xlsx"), sheetName="C.s",append=TRUE )
    write.xlsx2(all[all$analyte%in%C$analyte, ], paste(title, ".xlsx"), sheetName="N.b",append=TRUE )
    write.xlsx2(all[all$analyte%in%D$analyte, ], paste(title, ".xlsx"), sheetName="N.c",append=TRUE )
    write.xlsx2(all[all$analyte%in%E$analyte, ], paste(title, ".xlsx"), sheetName="S.b",append=TRUE )
    write.xlsx2(all[all$analyte%in%FF$analyte, ], paste(title, ".xlsx"), sheetName="S.c",append=TRUE )
    
    write.xlsx2(data, paste(title, ".xlsx"), sheetName="result",append=TRUE )
    #write.xlsx2(nosig, paste(title, ".xlsx"), sheetName="nosig",append=TRUE )
    #write.xlsx2(sig, paste(title, ".xlsx"), sheetName="sig",append=TRUE )
  }
  return (list(all=all,A=A,B=B,C=C,D=D,E=E,FF=FF, result=data,sig=sig, nosig=nosig))
}












library(xlsx)
fc <- 1.2
dox_A.s <- re%>% filter(A.s.t.p<0.05&(A.s.fc<=1/fc|A.s.fc>=fc)) %>% droplevels()
dox_B.s <- re%>% filter(B.s.t.p<0.05&(B.s.fc<=1/fc|B.s.fc>=fc)) %>% droplevels()
dox_C.s <- re%>% filter(C.s.t.p<0.05&(C.s.fc<=1/fc|C.s.fc>=fc)) %>% droplevels()

dox_N.b <- re%>% filter(N.b.t.p<0.05&(N.b.fc<=1/fc|N.b.fc>=fc)) %>% droplevels()
dox_N.c <- re%>% filter(N.c.t.p<0.05&(N.c.fc<=1/fc|N.c.fc>=fc)) %>% droplevels()

dox_S.b <- re%>% filter(S.b.t.p<0.05&(S.b.fc<=1/fc|S.b.fc>=fc)) %>% droplevels()
dox_S.c <- re%>% filter(S.c.t.p<0.05&(S.c.fc<=1/fc|S.c.fc>=fc)) %>% droplevels()

write.xlsx2(dox_A.s, "dox_fold1d2p05_04142018.xlsx", sheetName="dox_A.s",append=TRUE )
write.xlsx2(dox_B.s, "dox_fold1d2p05_04142018.xlsx", sheetName="dox_B.s",append=TRUE )
write.xlsx2(dox_C.s, "dox_fold1d2p05_04142018.xlsx", sheetName="dox_C.s",append=TRUE )

write.xlsx2(dox_N.b, "dox_fold1d2p05_04142018.xlsx", sheetName="dox_N.b",append=TRUE )
write.xlsx2(dox_N.c, "dox_fold1d2p05_04142018.xlsx", sheetName="dox_N.c",append=TRUE )

write.xlsx2(dox_S.b, "dox_fold1d2p05_04142018.xlsx", sheetName="dox_S.b",append=TRUE )
write.xlsx2(dox_S.c, "dox_fold1d2p05_04142018.xlsx", sheetName="dox_S.c",append=TRUE )

####Boxplot for HCG

names(sample_mean1)[1:4]
HCG<- somaid$SomaId[somaid$Target=="HCG"]
sample_mean1$Status <- factor(sample_mean1$Status, levels=c("Normal", "Abnormal"))
ggplot(sample_mean1, aes(TimePoint, SL001766, color=Status))+geom_boxplot(color=c("black", "red"))+
  labs(y="HCG")
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/11-29-2016/HCG boxplot.pdf", width=4, height=3)
p1<-ggplot(sample_mean1, aes(x=TimePoint, y=SL001766,color=Status))+
  theme_classic()+ 
  #facet_wrap( ~ protein, scales="free", nrow=2, ncol=5)+ #
  #stat_boxplot(geom ='errorbar',width=0.4,coef=100)+
  geom_boxplot(outlier.alpha=1, position=position_dodge(0.8), show.legend = T)+
  #facet_grid(protein ~Time, scales="free")+
  scale_color_manual(name="Status",values = c("black", "red"))+
  labs(y="HCG")+ 
  theme(axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.text.x=element_text(vjust=0.6,angle=0, color="black"),
        axis.text.y=element_text(color="black"))+
  #theme(axis.title=element_text(size=16),axis.text.y=element_text(size=14),axis.text.x=element_blank(),axis.text.x=element_text(angle=45))+
  theme(plot.title=element_text(size=10))+
  #theme(plot.margin=unit(c(10,15,10,15),"line"))+
  theme(legend.position="right")+
  theme(strip.background = element_blank(),strip.text = element_text(size=8))+
  theme(axis.ticks = element_line(color="black"))
p1
dev.off()

PGRPS<- somaid$SomaId[somaid$Target=="PGRP-S"]
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/11-29-2016/PGRP-Slog10 boxplot.pdf", width=4, height=3)
p1<-ggplot(sample_mean1, aes(x=TimePoint, y=log(SL004515,10),color=Status))+
  theme_classic()+ 
  #facet_wrap( ~ protein, scales="free", nrow=2, ncol=5)+ #
  #stat_boxplot(geom ='errorbar',width=0.4,coef=100)+
  geom_boxplot(outlier.alpha=1, position=position_dodge(0.8), show.legend = T)+
  #facet_grid(protein ~Time, scales="free")+
  scale_color_manual(name="Status",values = c("black", "red"))+
  labs(y="log10 PGRP-S")+ 
  theme(axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.text.x=element_text(vjust=0.6,angle=0, color="black"),
        axis.text.y=element_text(color="black"))+
  #theme(axis.title=element_text(size=16),axis.text.y=element_text(size=14),axis.text.x=element_blank(),axis.text.x=element_text(angle=45))+
  theme(plot.title=element_text(size=10))+
  #theme(plot.margin=unit(c(10,15,10,15),"line"))+
  theme(legend.position="right")+
  theme(strip.background = element_blank(),strip.text = element_text(size=8))+
  theme(axis.ticks = element_line(color="black"))
p1
dev.off()


Renin<- somaid$SomaId[somaid$Target=="Renin"]

pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/11-29-2016/Renin boxplot.pdf", width=4, height=3)
p1<-ggplot(sample_mean1, aes(x=TimePoint, y=SL000565,color=Status))+
  theme_classic()+ 
  #facet_wrap( ~ protein, scales="free", nrow=2, ncol=5)+ #
  #stat_boxplot(geom ='errorbar',width=0.4,coef=100)+
  geom_boxplot(outlier.alpha=1, position=position_dodge(0.8), show.legend = T)+
  #facet_grid(protein ~Time, scales="free")+
  scale_color_manual(name="Status",values = c("black", "red"))+
  labs(y="Renin")+ 
  theme(axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.text.x=element_text(vjust=0.6,angle=0, color="black"),
        axis.text.y=element_text(color="black"))+
  #theme(axis.title=element_text(size=16),axis.text.y=element_text(size=14),axis.text.x=element_blank(),axis.text.x=element_text(angle=45))+
  theme(plot.title=element_text(size=10))+
  #theme(plot.margin=unit(c(10,15,10,15),"line"))+
  theme(legend.position="right")+
  theme(strip.background = element_blank(),strip.text = element_text(size=8))+
  theme(axis.ticks = element_line(color="black"))
p1
dev.off()




Proteinase3<- somaid$SomaId[somaid$Target=="Proteinase-3"]

pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/11-29-2016/proteinase3 boxplot.pdf", width=4, height=3)
p1<-ggplot(sample_mean1, aes(x=TimePoint, y=SL004008,color=Status))+
  theme_classic()+ 
  #facet_wrap( ~ protein, scales="free", nrow=2, ncol=5)+ #
  #stat_boxplot(geom ='errorbar',width=0.4,coef=100)+
  geom_boxplot(outlier.alpha=1, position=position_dodge(0.8), show.legend = T)+
  #facet_grid(protein ~Time, scales="free")+
  scale_color_manual(name="Status",values = c("black", "red"))+
  labs(y="Proteinase-3")+ 
  theme(axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.text.x=element_text(vjust=0.6,angle=0, color="black"),
        axis.text.y=element_text(color="black"))+
  #theme(axis.title=element_text(size=16),axis.text.y=element_text(size=14),axis.text.x=element_blank(),axis.text.x=element_text(angle=45))+
  theme(plot.title=element_text(size=10))+
  #theme(plot.margin=unit(c(10,15,10,15),"line"))+
  theme(legend.position="right")+
  theme(strip.background = element_blank(),strip.text = element_text(size=8))+
  theme(axis.ticks = element_line(color="black"))
p1
dev.off()

names(sample_mean1)[1:4]
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/11-29-2016/proteinase3log10_status_time_boxplot.pdf", width=4, height=3)
p1<-ggplot(sample_mean1, aes(x=Status, y=log(SL004008, 10),color=TimePoint))+
  theme_classic()+ 
  #facet_wrap( ~ protein, scales="free", nrow=2, ncol=5)+ #
  #stat_boxplot(geom ='errorbar',width=0.4,coef=100)+
  geom_boxplot(outlier.alpha=1, position=position_dodge(0.8), show.legend = T)+
  #facet_grid(protein ~Time, scales="free")+
  scale_color_manual(name="TimePoint",values = c("black", "orange","red"))+
  labs(y="log10 Proteinase-3")+ 
  theme(axis.title=element_text(size=10),axis.text=element_text(size=8),
        axis.text.x=element_text(vjust=0.6,angle=0, color="black"),
        axis.text.y=element_text(color="black"))+
  #theme(axis.title=element_text(size=16),axis.text.y=element_text(size=14),axis.text.x=element_blank(),axis.text.x=element_text(angle=45))+
  theme(plot.title=element_text(size=10))+
  #theme(plot.margin=unit(c(10,15,10,15),"line"))+
  theme(legend.position="right")+
  theme(strip.background = element_blank(),strip.text = element_text(size=8))+
  theme(axis.ticks = element_line(color="black"))
p1
dev.off()










##check what outliers are:subject 11a ,b,c, 33c
p1<-ggplot(sample_mean1, aes(x=TimePoint, y=SL001766,color=Status,shape=SubjectID))+
  scale_shape_manual(values=c(97:119, 0:26))+
  #geom_dotplot(binaxis = "y", dotsize=0.4, position=position_dodge(0.8))+
  geom_jitter(alpha=0.5,size=3)+
  labs(color="status", shape="subject")
p1
sample_mean1$SubjectID <- factor(sample_mean1$SubjectID)
sample_mean1[sample_mean1$SubjectID==10,1:5]


boxplot(c(1:101))
##################################################################################
######09-27-2016 data analysis
##################################################################################
names(duplicates)[1:30]
names(protein)
names(sample1)
sample1<-sample[,26:1347]
sample$SampleDescription<-as.factor(sample$SampleDescription)

protein<-write.csv(protein,"C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/09-27-2016/FDA-16-003.hybNorm.plateScale.medNorm.calibrate.20160927protein1.csv") # all proteins
#sample quality check

sample_pca<-pca(log(sample1),ncomp=3,scale=TRUE)
colors<-c("orange","black","red","purple","gray", "black", "gold","violet", "tan")
plotIndiv(sample_pca,ind.names=sample$SampleId, group=sample$Status, ellipse = F,
          ellipse.level = 0.95,col.per.group = colors[1:3],cex=5,legend=TRUE, title="pca")
plotVar(res_pca, comp = c(1,2), var.names  =list(names(xlog)),title="pca")


#names(sample)
#str(protein)
#sample<-sample %>% filter(RowCheck=="PASS") %>% droplevels(SampleId) %>% summarise_each(funs(mean,sd))
sample1<-sample[,26:1347]
length(names(sample1))
#names(sample1)

#t-test




resultfoldp<-data.frame("analytes"=character(),"A.m.t.p"=numeric(),"A.m.fc"=numeric(),"A.s.t.p"=numeric(),"A.s.fc"=numeric(),
                                               "B.m.t.p"=numeric(),"B.m.fc"=numeric(),"B.s.t.p"=numeric(),"B.s.fc"=numeric(),
                                               "C.m.t.p"=numeric(),"C.m.fc"=numeric(),"C.s.t.p"=numeric(),"C.s.fc"=numeric(),
                                               "N.b.t.p"=numeric(),"N.b.fc"=numeric(),"N.c.t.p"=numeric(),"N.c.fc"=numeric(),
                                               "M.b.t.p"=numeric(),"M.b.fc"=numeric(),"M.c.t.p"=numeric(),"M.c.fc"=numeric(),
                                               "S.b.t.p"=numeric(),"S.b.fc"=numeric(),"S.c.t.p"=numeric(),"S.c.fc"=numeric(),stringsAsFactors=FALSE)
for (i in 1:1322){

  
  A.n<-sample1[sample$TreatStatus=="Normal_A",i]
  A.m<-sample1[sample$TreatStatus=="Moderate_A",i]
  A.s<-sample1[sample$TreatStatus=="Severe_A",i]
  A.m.fc<-mean(A.m)/mean(A.n)
  A.s.fc<-mean(A.s)/mean(A.n)
  A.m.t<-t.test(log(A.m,10), log(A.n,10),paired=FALSE,var.qual=FALSE)
  A.s.t<-t.test(log(A.s,10), log(A.n,10),paired=FALSE,var.qual=FALSE)
  
  B.n<-sample1[sample$TreatStatus=="Normal_B",i]
  B.m<-sample1[sample$TreatStatus=="Moderate_B",i]
  B.s<-sample1[sample$TreatStatus=="Severe_B",i]
  B.m.fc<-mean(B.m)/mean(B.n)
  B.s.fc<-mean(B.s)/mean(B.n)
  B.m.t<-t.test(log(B.m,10), log(B.n,10),paired=FALSE,var.qual=FALSE)
  B.s.t<-t.test(log(B.s,10), log(B.n,10),paired=FALSE,var.qual=FALSE)
  
  C.n<-sample1[sample$TreatStatus=="Normal_C",i]
  C.m<-sample1[sample$TreatStatus=="Moderate_C",i]
  C.s<-sample1[sample$TreatStatus=="Severe_C",i]
  C.m.fc<-mean(C.m)/mean(C.n)
  C.s.fc<-mean(C.s)/mean(C.n)
  C.m.t<-t.test(log(C.m,10), log(C.n,10),paired=FALSE,var.qual=FALSE)
  C.s.t<-t.test(log(C.s,10), log(C.n,10),paired=FALSE,var.qual=FALSE)
  
  
  #normalA<-sample$SampleId[sample$TreatStatus=="Normal_A"]
  #normalB<-sample$SampleId[sample$TreatStatus=="Normal_B"]
  #normalC<-sample$SampleId[sample$TreatStatus=="Normal_C"]
  
  #normalA[c(-5,-18,-21)]
  #normalB[c(-5,-15)]
 
  #normalA[c(-3,-17,-22)]
  #normalC[c(-5,-15)]
  
  N.b.fc<-mean(B.n[c(-5,-15)])/mean(A.n[c(-5,-18,-21)])
  N.c.fc<-mean(C.n[c(-5,-15)])/mean(A.n[c(-3,-17,-22)])
  N.b.t<-t.test(log(B.n[c(-5,-15)],10), log(A.n[c(-5,-18,-21)],10),paired=TRUE,var.qual=FALSE)
  N.c.t<-t.test(log(C.n[c(-5,-15)],10), log(A.n[c(-3,-17,-22)],10),paired=TRUE,var.qual=FALSE)
  
  M.b.fc<-mean(B.m)/mean(A.m)
  M.c.fc<-mean(C.m)/mean(A.m)
  M.b.t<-t.test(log(B.m,10), log(A.m,10),paired=TRUE,var.qual=FALSE)
  M.c.t<-t.test(log(C.m,10), log(A.m,10),paired=TRUE,var.qual=FALSE)
  
  S.b.fc<-mean(B.s)/mean(A.s)
  S.c.fc<-mean(C.s)/mean(A.s)
  S.b.t<-t.test(log(B.s,10), log(A.s,10),paired=TRUE,var.qual=FALSE)
  S.c.t<-t.test(log(C.s,10), log(A.s,10),paired=TRUE,var.qual=FALSE)
  
  
  resultfoldp[i,]<-c(names(sample1)[i],A.m.t$p.value,A.m.fc,A.s.t$p.value,A.s.fc,
                                       B.m.t$p.value,B.m.fc,B.s.t$p.value,B.s.fc,
                                       C.m.t$p.value,C.m.fc,C.s.t$p.value,C.s.fc,
                                       N.b.t$p.value,N.b.fc,N.c.t$p.value,N.c.fc,
                                       M.b.t$p.value,M.b.fc,M.c.t$p.value,M.c.fc,
                                       S.b.t$p.value,S.b.fc,S.c.t$p.value,S.c.fc)
                     
  
}








names(resultfoldp)

resultfoldp[,"A.m.t.fdr"]<-p.adjust(resultfoldp$A.m.t.p,"BH")
resultfoldp[,"A.s.t.fdr"]<-p.adjust(resultfoldp$A.s.t.p,"BH")
resultfoldp[,"B.m.t.fdr"]<-p.adjust(resultfoldp$B.m.t.p,"BH")
resultfoldp[,"B.s.t.fdr"]<-p.adjust(resultfoldp$B.s.t.p,"BH")
resultfoldp[,"C.m.t.fdr"]<-p.adjust(resultfoldp$C.m.t.p,"BH")
resultfoldp[,"C.s.t.fdr"]<-p.adjust(resultfoldp$C.s.t.p,"BH")

resultfoldp[,"N.b.t.fdr"]<-p.adjust(resultfoldp$N.b.t.p,"BH")
resultfoldp[,"N.c.t.fdr"]<-p.adjust(resultfoldp$N.c.t.p,"BH")
resultfoldp[,"M.b.t.fdr"]<-p.adjust(resultfoldp$M.b.t.p,"BH")
resultfoldp[,"M.c.t.fdr"]<-p.adjust(resultfoldp$M.c.t.p,"BH")
resultfoldp[,"S.b.t.fdr"]<-p.adjust(resultfoldp$S.b.t.p,"BH")
resultfoldp[,"S.c.t.fdr"]<-p.adjust(resultfoldp$S.c.t.p,"BH")



write.csv(resultfoldp,"C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/09-27-2016/doxfoldp_paired.csv")

re<-read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/09-27-2016/doxfoldp_paired.csv")
re[,"ColCheck"]<-protein$ColCheck[which(re$analytes %in% protein$Target)]

write.csv(re,"C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/09-27-2016/doxfoldp_paired.csv")
str(re)
names(re)

###box plot
fc<-1.2
p<-0.05
re<-re[re$ColCheck=="PASS",]
nrow(re)
names(re)
doxAmp05fc1d2<-re%>% filter((re[,4]<=p&(re[,5]>=fc|re[,5]<=1/fc)))%>%droplevels()  
doxAsp05fc1d2<-re%>% filter(re[,6]<=p&(re[,7]>=fc|re[,7]<=1/fc)) %>%droplevels()   
doxAp05fc1d2<-re%>% filter((re[,4]<=p&(re[,5]>=fc|re[,5]<=1/fc))|(re[,6]<=p&(re[,7]>=fc|re[,7]<=1/fc))) %>%droplevels()  
doxAmsp05fc1d2<-re%>% filter((re[,4]<=p&(re[,5]>=fc|re[,5]<=1/fc))&(re[,6]<=p&(re[,7]>=fc|re[,7]<=1/fc))) %>%droplevels()
                          
doxBmp05fc1d2<-re%>% filter((re[,8]<=p&(re[,9]>=fc|re[,9]<=1/fc))) %>%droplevels()
doxBsp05fc1d2<-re%>% filter(re[,10]<=p&(re[,11]>=fc|re[,11]<=1/fc))%>%droplevels()
doxBp05fc1d2<-re%>% filter((re[,8]<=p&(re[,9]>=fc|re[,9]<=1/fc))|(re[,10]<=p&(re[,11]>=fc|re[,11]<=1/fc))) %>%droplevels()
doxBmsp05fc1d2<-re%>% filter((re[,8]<=p&(re[,9]>=fc|re[,9]<=1/fc))&(re[,10]<=p&(re[,11]>=fc|re[,11]<=1/fc))) %>%droplevels()



doxCmp05fc1d2<-re%>% filter((re[,12]<=p&(re[,13]>=fc|re[,13]<=1/fc)))%>%droplevels()
doxCsp05fc1d2<-re%>% filter(re[,14]<=p&(re[,15]>=fc|re[,15]<=1/fc)) %>%droplevels()
doxCp05fc1d2<-re%>% filter((re[,12]<=p&(re[,13]>=fc|re[,13]<=1/fc))|(re[,14]<=p&(re[,15]>=fc|re[,15]<=1/fc))) %>%droplevels()
doxCmsp05fc1d2<-re%>% filter((re[,12]<=p&(re[,13]>=fc|re[,13]<=1/fc))&(re[,14]<=p&(re[,15]>=fc|re[,15]<=1/fc))) %>%droplevels()



doxnBp05fc1d2<-re%>% filter((re[,16]<=p&(re[,17]>=fc|re[,17]<=1/fc))) %>%droplevels()
doxnCp05fc1d2<-re%>% filter(re[,18]<=p&(re[,19]>=fc|re[,19]<=1/fc))%>%droplevels()
doxnp05fc1d2<-re%>% filter((re[,16]<=p&(re[,17]>=fc|re[,17]<=1/fc))|(re[,18]<=p&(re[,19]>=fc|re[,19]<=1/fc))) %>%droplevels()
doxnBCp05fc1d2<-re%>% filter((re[,16]<=p&(re[,17]>=fc|re[,17]<=1/fc))&(re[,18]<=p&(re[,19]>=fc|re[,19]<=1/fc))) %>%droplevels()

doxmBp05fc1d2<-re%>% filter((re[,20]<=p&(re[,21]>=fc|re[,21]<=1/fc))) %>%droplevels()
doxmCp05fc1d2<-re%>% filter((re[,22]<=p&(re[,23]>=fc|re[,23]<=1/fc))) %>%droplevels()
doxmp05fc1d2<-re%>% filter((re[,20]<=p&(re[,21]>=fc|re[,21]<=1/fc))|(re[,22]<=p&(re[,23]>=fc|re[,23]<=1/fc))) %>%droplevels()
doxmBCp05fc1d2<-re%>% filter((re[,20]<=p&(re[,21]>=fc|re[,21]<=1/fc))&(re[,22]<=p&(re[,23]>=fc|re[,23]<=1/fc))) %>%droplevels()


doxsBp05fc1d2<-re%>% filter((re[,24]<=p&(re[,25]>=fc|re[,25]<=1/fc))) %>%droplevels()
doxsCp05fc1d2<-re%>% filter((re[,26]<=p&(re[,27]>=fc|re[,27]<=1/fc))) %>%droplevels()
doxsp05fc1d2<-re%>% filter((re[,24]<=p&(re[,25]>=fc|re[,25]<=1/fc))|(re[,26]<=p&(re[,27]>=fc|re[,27]<=1/fc))) %>%droplevels()
doxsBCp05fc1d2<-re%>% filter((re[,24]<=p&(re[,25]>=fc|re[,25]<=1/fc))&(re[,26]<=p&(re[,27]>=fc|re[,27]<=1/fc))) %>%droplevels()




nrow(doxsBCp05fc1d2)


doxp05fc1d2<-re%>% filter((re[,4]<=p&(re[,5]>=fc|re[,5]<=1/fc))|(re[,6]<=p&(re[,7]>=fc|re[,7]<=1/fc))|
                             (re[,8]<=p&(re[,9]>=fc|re[,9]<=1/fc))|(re[,10]<=p&(re[,11]>=fc|re[,11]<=1/fc))|
                             (re[,12]<=p&(re[,13]>=fc|re[,13]<=1/fc))|(re[,14]<=p&(re[,15]>=fc|re[,15]<=1/fc))|
                             (re[,16]<=p&(re[,17]>=fc|re[,17]<=1/fc))|(re[,18]<=p&(re[,19]>=fc|re[,19]<=1/fc))| 
                             (re[,20]<=p&(re[,21]>=fc|re[,21]<=1/fc))|(re[,22]<=p&(re[,23]>=fc|re[,23]<=1/fc))| 
                             (re[,24]<=p&(re[,25]>=fc|re[,25]<=1/fc))|(re[,26]<=p&(re[,27]>=fc|re[,27]<=1/fc))) %>%droplevels()






write.csv(doxp.05fc2,"C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/09-27-2016/doxfoldp_pairedp05fc2.csv")

y<-sample[,paste(doxp05fc1d2$analytes)]
names(y)
dim(y)
dim(sample)


library(gridExtra)


pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/09-27-2016/doxp05fc1.2.pdf", paper="a4r", width = 0, height = 0)


for (i in 1:618){

    Amtp<-format(re$A.m.t.p[re$analytes==names(y)[i]],digits=1,nsmall=3)
    Amtfdr<-format(re$A.m.t.fdr[re$analytes==names(y)[i]],digits=1,nsmall=3)
    
    Astp<-format(re$A.s.t.p[re$analytes==names(y)[i]],digits=1,nsmall=3)
    Astfdr<-format(re$A.s.t.fdr[re$analytes==names(y)[i]],digits=1,nsmall=3)
    
    Bmtp<-format(re$B.m.t.p[re$analytes==names(y)[i]],digits=1,nsmall=3)
    Bmtfdr<-format(re$B.m.t.fdr[re$analytes==names(y)[i]],digits=1,nsmall=3)
 
    Bstp<-format(re$B.s.t.p[re$analytes==names(y)[i]],digits=1,nsmall=3)
    Bstfdr<-format(re$B.s.t.fdr[re$analytes==names(y)[i]],digits=1,nsmall=3)
    
    Cmtp<-format(re$C.m.t.p[re$analytes==names(y)[i]],digits=1,nsmall=3)
    Cmtfdr<-format(re$C.m.t.fdr[re$analytes==names(y)[i]],digits=1,nsmall=3)
    
    Cstp<-format(re$C.s.t.p[re$analytes==names(y)[i]],digits=1,nsmall=3)
    Cstfdr<-format(re$C.s.t.fdr[re$analytes==names(y)[i]],digits=1,nsmall=3)
    
 
    Nbtp<-format(re$N.b.t.p[re$analytes==names(y)[i]],digits=1,nsmall=3)
    Nbtfdr<-format(re$N.b.t.fdr[re$analytes==names(y)[i]],digits=1,nsmall=3)
    
    Nctp<-format(re$N.c.t.p[re$analytes==names(y)[i]],digits=1,nsmall=3)
    Nctfdr<-format(re$N.c.t.fdr[re$analytes==names(y)[i]],digits=1,nsmall=3)
    
    Mbtp<-format(re$M.b.t.p[re$analytes==names(y)[i]],digits=1,nsmall=3)
    Mbtfdr<-format(re$M.b.t.fdr[re$analytes==names(y)[i]],digits=1,nsmall=3)
    
    Mctp<-format(re$M.c.t.p[re$analytes==names(y)[i]],digits=1,nsmall=3)
    Mctfdr<-format(re$M.c.t.fdr[re$analytes==names(y)[i]],digits=1,nsmall=3)
    
    Sbtp<-format(re$S.b.t.p[re$analytes==names(y)[i]],digits=1,nsmall=3)
    Sbtfdr<-format(re$S.b.t.fdr[re$analytes==names(y)[i]],digits=1,nsmall=3)
    
    Sctp<-format(re$S.c.t.p[re$analytes==names(y)[i]],digits=1,nsmall=3)
    Sctfdr<-format(re$S.c.t.fdr[re$analytes==names(y)[i]],digits=1,nsmall=3)

  
  pvalue<-paste(c("group:", "A_m.n","A_s.n","B_m.n","B_s.n","C_m.n","C_s.n","n_B.A","n_C.A","m_B.A","m_C.A","s_B.A","s_C.A", 
                  "\np:      ",    Amtp,  Astp,  Bmtp,  Bstp,   Cmtp,   Cstp,   Nbtp,   Nctp,  Mbtp,     Mctp,  Sbtp,  Sctp,
                  "\nfdr:",Amtfdr, Astfdr, Bmtfdr, Bstfdr, Cmtfdr, Cstfdr, Nbtfdr, Nctfdr, Mbtfdr, Mctfdr, Sbtfdr, Sctfdr)
                     , collapse="  ", sep="      ")
 


  #plot for concentration
  #shape c(0:25,33:38)
  
  plot1<-ggplot(sample, aes(x=TreatStatus, y=y[, i], group=TreatStatus, shape=SampleDescription))+
    theme_bw()+ 
    geom_boxplot(outlier.color=NA, show.legend = F)+
    geom_jitter(size=2, fill="Red")+
    ggtitle(paste("  ",pvalue))+
    scale_shape_manual(name="SubjectID",values=c(0:25,65:72), labels=levels(sample$SampleDescription))+
    labs(x="", y=names(y)[i])+ 
    theme(axis.title=element_text(size=16),axis.text=element_text(size=14),axis.text.x=element_text(angle=90))+
    theme(plot.title=element_text(size=14))+
    theme(plot.margin=unit(c(1,1,1,1),"line"))

  grid.arrange(plot1, nrow=1)
 
}


dev.off() 

#
#boxplot overlap
#


protein5<- read.csv("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/09-27-2016/dox_fiveprotein_boxplot.csv", sep=",", header=TRUE)
names(protein5)
head(protein5,5)
y<-protein5[,6:10]
str(protein5)
ggplot(protein5, aes(x=Time, y=y[, 1], group=TreatStatus,shape=TreatStatus1))+
  theme_bw()+ 
  geom_boxplot(aes(fill=TreatStatus1))+
  geom_points(aes(y=y[, 1], group=TreatStatus),size=4, fill="Red")+
  #ggtitle(paste("  ",pvalue))+
  scale_shape_manual(name="SubjectID",values=c(19,21,12), labels=levels(protein5$TreatStatus1))+
  labs(x="", y=names(y)[i])+ 
  theme(axis.title=element_text(size=16),axis.text=element_text(size=14),axis.text.x=element_text(angle=90))+
  theme(plot.title=element_text(size=14))+
  theme(plot.margin=unit(c(1,1,1,1),"line"))
protein5$TreatStatus<-factor(protein5$TreatStatus, levels=c("Normal_A", "Moderate_A","Severe_A","Normal_B", "Moderate_B", "Severe_B","Normal_C","Moderate_C","Severe_C"),
                             labels=c("A_n", "A_m","A_s","B_n", "B_m","B_s","C_n", "C_m","C_s"))
protein5$TreatStatus1<-factor(protein5$TreatStatus1, levels=c("Normal", "Moderate","Severe"))
library(gridExtra)
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/09-27-2016/dox_5proteinslog10.pdf", paper="a4r", width = 0, height = 0)

for (i in 1:5){
i=1
  plot1<-ggplot(protein5, aes(x=TreatStatus1, y=log(y[, i],10),fill=TreatStatus1))+
  theme_bw()+ 
  geom_boxplot(position=position_dodge(0.8), show.legend = T)+
  #geom_dotplot(binaxis = "y", dotsize=0.4, position=position_dodge(0.8))+
  geom_jitter(alpha=0.5,position=position_jitterdodge(0.8))+
  #facet_wrap( ~ Time, scales="free")+
  #geom_point(position=position_dodge(0.8))+
  scale_fill_manual(name="Status",values = c("green", "pink",  "coral"))+
  #ggtitle(paste("  ",pvalue))+
  #scale_shape_manual(name="SubjectID",values=c(19,21,12), labels=levels(protein5$TreatStatus1))+
  labs(x=" ", y=paste("log10 ",names(y)[i]))+ 
  facet_grid(.~Time)+
  theme(axis.title=element_text(size=16),axis.text=element_text(size=14),axis.text.x=element_text(angle=45))+
  theme(plot.title=element_text(size=14))+
  theme(plot.margin=unit(c(8,15,8,15),"line"))

grid.arrange(plot1, nrow=1)
}
dev.off() 
pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/09-27-2016/dox_5proteinslog10_grid.pdf", paper="a4r", width = 0, height = 0)
for (i in 1:5){
  
  plot1<-ggplot(protein5, aes(x=TreatStatus1, y=log(y[, i],10),fill=TreatStatus1))+
    theme_bw()+ 
    geom_boxplot(position=position_dodge(NULL), show.legend = T)+
    #geom_dotplot(binaxis = "y", dotsize=0.4, position=position_dodge(0.8))+
    geom_jitter(alpha=0.5,position=position_jitterdodge(jitter.width=3,jitter.height=0,dodge.width = NULL))+
    #facet_wrap( ~ Time, scales="free")+
    #geom_point(position=position_dodge(0.8))+
    scale_fill_manual(name="Status",values = c("green", "pink",  "coral"))+
    #ggtitle(paste("  ",pvalue))+
    #scale_shape_manual(name="SubjectID",values=c(19,21,12), labels=levels(protein5$TreatStatus1))+
    labs(x=" ", y=paste("log10 ",names(y)[i]))+ 
    facet_grid(.~Time)+
    #theme(axis.title=element_text(size=16),axis.text=element_text(size=14),axis.text.x=element_text(angle=45))+
    theme(axis.title=element_text(size=16),axis.text.y=element_text(size=14),axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),axis.text.x=element_text(angle=45))+
    theme(plot.title=element_text(size=14))+
    theme(plot.margin=unit(c(8,15,8,15),"line"))+
    theme(legend.position="top",legend.title=element_blank())
  
  grid.arrange(plot1, nrow=1)
}
dev.off() 

pdf("C:/zhijuncao/cardiotoxicity/SomaLogic Returned Data/09-27-2016/dox_5proteinslog10_dodge.pdf", paper="a4r", width = 0, height = 0)
for (i in 1:5){

  plot1<-ggplot(protein5, aes(x=Time, y=log(y[, i],10),fill=TreatStatus1))+
    theme_bw()+ 
    geom_boxplot(position=position_dodge(0.8), show.legend = T)+
    #geom_dotplot(binaxis = "y", dotsize=0.4, position=position_dodge(0.8))+
    geom_jitter(alpha=0.5, position=position_jitterdodge(jitter.width=1,jitter.height=0,dodge.width = 0.8))+
    #facet_wrap( ~ Time, scales="free")+
    #geom_point(position=position_dodge(0.8))+
    scale_fill_manual(name="Status",values = c("green", "pink",  "coral"))+
    #ggtitle(paste("  ",pvalue))+
    #scale_shape_manual(name="SubjectID",values=c(19,21,12), labels=levels(protein5$TreatStatus1))+
    labs(x=" ", y=paste("log10 ",names(y)[i]))+ 
    #facet_grid(.~Time)+
    theme(axis.title=element_text(size=16),axis.text=element_text(size=14),axis.text.x=element_text(angle=0))+
    #theme(axis.title=element_text(size=16),axis.text.y=element_text(size=14),axis.text.x=element_blank(),axis.text.x=element_text(angle=45))+
    theme(plot.title=element_text(size=14))+
    theme(plot.margin=unit(c(10,15,10,15),"line"))+
    theme(legend.position="right")
  
  grid.arrange(plot1, nrow=1)
}
dev.off() 

protein5$TreatStatus
protein5[protein5$Time=="A","Timecolor"]<-"black"
protein5[protein5$Time=="B","Timecolor"]<-"green"
protein5[protein5$Time=="C","Timecolor"]<-"red"
head(protein5,5)
ggplot(protein5, aes(x=TreatStatus, y=log(y[, i],10),shape=protein5$Timecolor,fill=TreatStatus1))+
  theme_bw()+ 
  geom_boxplot(position=position_dodge(0.8), show.legend = T)+
  #geom_dotplot(binaxis = "y", dotsize=0.4, position=position_dodge(0.8))+
  geom_jitter(size=2)+
  #facet_wrap( ~ Time, scales="free")+
  #geom_point(position=position_dodge(0.8))+
  scale_fill_manual(name="Status",values = c("green", "pink",  "coral"))+
  #ggtitle(paste("  ",pvalue))+
  scale_shape_manual(name="SubjectID",values=c(19,21,12), labels=levels(protein5$Time))+
  scale_x_discrete(labels=c("","A"," "," ", "B"," ","","C"," "))+
  labs(x="", y=paste("log10 ",names(y)[i]))+ 
  theme(axis.title=element_text(size=16),axis.text=element_text(size=14),axis.text.x=element_text(angle=0))+
  theme(plot.title=element_text(size=14))+
  theme(plot.margin=unit(c(8,5,8,5),"line"))





#classfication analysis
sampleA<-sample %>% filter(SampleGroup=="A")%>%droplevels()
sampleB<-sample %>% filter(SampleGroup=="B")%>%droplevels()
sampleC<-sample %>% filter(SampleGroup=="C")%>%droplevels()

samplen<-sample %>% filter(Status=="Normal")%>%droplevels()
samplem<-sample %>% filter(Status=="Moderate")%>%droplevels()
samples<-sample %>% filter(Status=="Severe")%>%droplevels()

xlogA<-log(sampleA[,paste(doxAp05fc1d2$analytes[doxAp05fc1d2$ColCheck=="PASS"])],10)
xlogB<-log(sampleB[,paste(doxBp05fc1d2$analytes[doxBp05fc1d2$ColCheck=="PASS"])],10)
xlogC<-log(sampleC[,paste(doxCp05fc1d2$analytes[doxCp05fc1d2$ColCheck=="PASS"])],10)

xlogn<-log(samplen[,paste(doxnp05fc1d2$analytes[doxnp05fc1d2$ColCheck=="PASS"])],10)
xlogm<-log(samplem[,paste(doxmp05fc1d2$analytes[doxmp05fc1d2$ColCheck=="PASS"])],10)
xlogs<-log(samples[,paste(doxsp05fc1d2$analytes[doxsp05fc1d2$ColCheck=="PASS"])],10)

colors<-c("orange","black","red","green","cyan", "gold","violet", "purple","gray")

#####heat status
#sampleA

res_pcaA<-pca(xlogA, ncomp=3,scale=TRUE)

plotIndiv(res_pcaA,ind.names=sampleA$SampleId, group=sampleA$TreatStatus, ellipse = TRUE,
          ellipse.level = 0.95,col.per.group = colors[1:3],cex=5,legend=TRUE, title="pca")
plotVar(res_pca, comp = c(1,2), var.names  =list(names(xlog)),title="pca")

#sampleB
res_pcaB<-pca(xlogB, ncomp=3,scale=TRUE)
plotIndiv(res_pcaB,ind.names=sampleB$SampleId, group=sampleB$TreatStatus, ellipse = TRUE,
          ellipse.level = 0.95,col.per.group = colors[1:3],cex=5,legend=TRUE, title="pca")
plotVar(res_pca, comp = c(1,2), var.names  =list(names(xlog)),title="pca")

#sampleC
res_pcaC<-pca(xlogC, ncomp=3,scale=TRUE)
plotIndiv(res_pcaC,ind.names=sampleC$SampleId, group=sampleC$TreatStatus, ellipse = TRUE,
          ellipse.level = 0.95,col.per.group = colors[1:3],cex=5,legend=TRUE, title="pca")

#plsda multilevel
re.splsdaA<-splsda(xlogA,Y=sampleA$TreatStatus,ncomp=3, keepX = rep(202,3),logratio = "none", multilevel=NULL)
plotIndiv(re.splsdaA, comp = 1:2, ind.names = sampleA$SampleId, ellipse = TRUE,
          ellipse.level = 0.95,pch=16, col.per.group=colors[1:3], 
          cex=6,group=sampleA$TreatStatus,legend=T,title="pls-da")
plotVar(re.splsdaA, comp = c(1,2), var.names  =list(names(xlogA)),title="A: pls-da", overlap=T)


re.splsdaAerror<-perf(re.splsdaA,validation = "loo")
re.splsdaAerror$error.rate
re.splsdaAerror$error.rate.class

re.splsdaB<-splsda(xlogB,Y=sampleB$TreatStatus,ncomp=3, keepX = rep(96,3),logratio = "none", multilevel=NULL)
plotIndiv(re.splsdaB, comp = 1:2, ind.names = sampleB$SampleId, ellipse = TRUE,
          ellipse.level = 0.95,pch=16, col.per.group=colors[1:3], 
          cex=6,group=sampleB$TreatStatus,legend=T,title="pls-da")
plotVar(re.splsdaB, comp = c(1,2), var.names  =list(names(xlogB)),title="B: pls-da", overlap=T)


re.splsdaBerror<-perf(re.splsdaB,validation = "loo")
re.splsdaBerror$error.rate
re.splsdaBerror$error.rate.class

re.splsdaC<-splsda(xlogC,Y=sampleC$TreatStatus,ncomp=3, keepX = rep(92,3),logratio = "none", multilevel=NULL)
plotIndiv(re.splsdaC, comp = 1:2, ind.names = sampleC$SampleId, ellipse = TRUE,
          ellipse.level = 0.95,pch=16, col.per.group=colors[1:3], 
          cex=6,group=sampleC$TreatStatus,legend=T,title="pls-da")
plotVar(re.splsdaC, comp = c(1,2), var.names  =list(names(xlogC)),title="C: pls-da", overlap=T)

re.splsdaCerror<-perf(re.splsdaC,validation = "loo")
re.splsdaCerror$error.rate
re.splsdaCerror$error.rate.class



####time course
#samplen
designn <- data.frame(sample2 = samplen$SampleDescription)
res_pcan<-pca(xlogn, ncomp=3,scale=TRUE, multilevel=designn)
plotIndiv(res_pcan,ind.names=samplen$SampleId, group=samplen$TreatStatus, ellipse = TRUE,
          ellipse.level = 0.95,col.per.group = colors[1:3],cex=5,legend=TRUE, title="multilevel pca")

plotVar(res_pca, comp = c(1,2), var.names  =list(names(xlog)),title="multilevel pca")

#samplem
designm <- data.frame(sample2 = samplem$SampleDescription)
res_pcam<-pca(xlogm, ncomp=3,scale=TRUE, multilevel=designm)
plotIndiv(res_pcam,ind.names=samplem$SampleId, group=samplem$TreatStatus, ellipse = TRUE,
          ellipse.level = 0.95,col.per.group = colors[1:3],cex=5,legend=TRUE, title="multilevel pca")
plotVar(res_pca, comp = c(1,2), var.names  =list(names(xlog)),title="multilevel pca")

#samples
designs <- data.frame(sample2 = samples$SampleDescription)
res_pcas<-pca(xlogs, ncomp=3,scale=TRUE, multilevel=designs)
plotIndiv(res_pcas,ind.names=samples$SampleId, group=samples$TreatStatus, ellipse = TRUE,
          ellipse.level = 0.95,col.per.group = colors[1:3],cex=5,legend=TRUE, title="multilevel pca")

#plsda
designn <- data.frame(sample2 = samplen$SampleDescription)
re.splsdan<-splsda(xlogn,Y=samplen$SampleGroup,ncomp=3, keepX = rep(268,3),logratio = "none", multilevel=designn)
plotIndiv(re.splsdan, comp = 1:2, ind.names = samplen$SampleId, ellipse = TRUE,
          ellipse.level = 0.95,pch=16, col.per.group=colors[1:3], 
          cex=6,group=samplen$TreatStatus,legend=T,title="multilevel pls-da")
plotVar(re.splsdan, comp = c(1,2), var.names  =list(names(xlogn)),title="Normal: multilevel pls-da", overlap=T)

re.splsdanerror<-perf(re.splsdan,validation = "loo")
re.splsdanerror$error.rate
re.splsdanerror$error.rate.class

designm <- data.frame(sample2 = samplem$SampleDescription)
re.splsdam<-splsda(xlogm,Y=samplem$SampleGroup,ncomp=3, keepX = rep(95,3),logratio = "none", multilevel=designm)
plotIndiv(re.splsdam, comp = 1:2, ind.names = samplem$SampleId, ellipse = TRUE,
          ellipse.level = 0.95,pch=16, col.per.group=colors[1:3], 
          cex=6,group=samplem$TreatStatus,legend=T,title="multilevel pls-da")
plotVar(re.splsdam, comp = c(1,2), var.names  =list(names(xlogm)),title="Moderate: multilevel pls-da", overlap=T)

re.splsdamerror<-perf(re.splsdam,validation = "loo")
re.splsdamerror$error.rate
re.splsdamerror$error.rate.class

designs <- data.frame(sample2 = samples$SampleDescription)
re.splsdas<-splsda(xlogs,Y=samples$SampleGroup,ncomp=3, keepX = rep(111,3),logratio = "none", multilevel=designs)
plotIndiv(re.splsdas, comp = 1:2, ind.names = samples$SampleId, ellipse = TRUE,
          ellipse.level = 0.95,pch=16, col.per.group=colors[1:3], 
          cex=6,group=samples$TreatStatus,legend=T,title="multilevel pls-da")
plotVar(re.splsdas, comp = c(1,2), var.names  =list(names(xlogs)),title="Severe: multilevel pls-da", overlap=T)

re.splsdaserror<-perf(re.splsdas,validation = "loo")
re.splsdaserror$error.rate
re.splsdaserror$error.rate.class









plotVar(re.splsda, comp = c(1,2), var.names  =list(names(xlog)),title="pls-da", overlap=F)
plotLoadings(re.splsda, contrib="max", method="median", comp=2, title="pls-da")


#day8
group<-sampleday8$status
xlog<-log(sampleday8[,paste(day8.p05$analytes[day8.p05$ColCheck=="PASS"])],10)
res_pca<-pca(xlog, ncomp=3,scale=TRUE)
plotIndiv(res_pca,ind.names=group, group=sampleday8$status, ellipse = TRUE,
          ellipse.level = 0.95,pch=16,col.per.group = c("red", "green"),cex=6,legend=TRUE, title="pca")
plotVar(res_pca, comp = c(1,2), var.names  =list(names(xlog)),title="pca")
plotLoadings(res_pca, contrib="max", method="median", comp=1)


res_ipca<-ipca(xlog, ncomp=3,scale=TRUE)
plotIndiv(res_ipca,ind.names=group, group=group, ellipse = TRUE,
          ellipse.level = 0.95,pch=16,col.per.group = c("red", "green"),cex=6,legend=TRUE, title="ipca")

re.splsda<-splsda(xlog,Y=group,ncomp=3, keepX = c(153,153,153),logratio = "none", multilevel=NULL)
plotIndiv(re.splsda, comp = 1:2, ind.names = group, ellipse = TRUE,
          ellipse.level = 0.95,pch=16,col.per.group=c("red", "green"), 
          cex=6,group=group,legend=T,title="pls-da")
plotVar(re.splsda, comp = c(1,2), var.names  =list(names(xlog)),title="pls-da")
plotLoadings(re.splsda, contrib="max", method="median", comp=1, ndisplay=50, title="pls-da")
plotLoadings(re.splsda, contrib="max", method="median", comp=2, ndisplay=50, title="pls-da")



#day1 and 8
xlog<-log(sample1[,paste(day1.8.p05$analytes[day1.8.p05$ColCheck=="PASS"])],10)
names(xlog)
group<-sample$status
res_pca<-pca(xlog, ncomp=3,scale=TRUE)
plotIndiv(res_pca,ind.names=sample$RowCheck, group=sample$status, ellipse = TRUE,
          ellipse.level = 0.95,pch=16, col.per.group = c("red", "green", "blue", "purple"),cex=6,legend=TRUE, title="pca")

plotVar(res_pca, comp = c(1,2), var.names  =list(names(xlog)),overlap=F)
plotLoadings(res_pca, contrib="max", method="median", comp=1)


res_ipca<-ipca(xlog, ncomp=3,scale=TRUE)
plotIndiv(res_ipca,ind.names=group, group=group, ellipse = TRUE,
          ellipse.level = 0.95,pch=16,col.per.group = c("red", "green"),cex=6,legend=TRUE, title="ipca")

re.splsda<-splsda(xlog,Y=group,ncomp=3, keepX = c(199,199,199),logratio = "none", multilevel=NULL)
plotIndiv(re.splsda, comp = 1:2, ind.names = group, ellipse = TRUE,
          ellipse.level = 0.95,pch=16,col.per.group=c("red", "green","blue","purple"), 
          cex=6,group=group,legend=T,title="pls-da")


#heatmap
heatmap.2(as.matrix(sampleday1a[,names(y)]),srtCol=0,  adjCol = c(0.5,1), 
          key=F,margins = c(4, 10), dendrogram="column",
          reorderfun=function(d, w) reorder(d, w, agglo.FUN =mean),
          labRow=sampleday1$status,
          cexRow = 0.1+ 1/log10(nrow(sample[,23:73])),
          cexCol = 0.2 + 1/log10(ncol(sample[,23:73])))

## volcanoplot
names(re)

fc<-1.2
sigcut <- data.frame(xx=seq(-2.5, 2.5, length.out=20), yy=(-log10(0.05)))
sigcutfc<-data.frame(xxp=log2(1.2),xxn=log2(1/1.2),yyp=seq(0, 3, length.out=40))

## A_m.n

re[re$A.m.t.p  >= 0.05, "colgroup"] <- "black"
re[re$A.m.t.p  < 0.05, "colgroup"] <- "orange"
re[re$A.m.t.p  < 0.05 & re$A.m.fc >=fc, "colgroup"] <- "red" 
re[re$A.m.t.p  < 0.05 & re$A.m.fc <= 1/fc , "colgroup"] <- "green" 


  ggplot(data=re, aes(log2(A.m.fc), -log10(A.m.t.p)))+
  geom_point(size=2,color=re$colgroup)+
  geom_text_repel(aes(label=analytes), data=subset(re,colgroup!="black"), size=3, col='black')+
  geom_line(data=sigcut, aes(xx, yy),lty=3,lwd=1)+
  geom_line(data=sigcutfc, aes(xxp, yyp),lty=3,lwd=1)+
  geom_line(data=sigcutfc, aes(xxn, yyp),lty=3,lwd=1)+
  ggtitle("A: Moderate vs Normal (p 0.05, fc 1.2)") +
  labs(x="log2(Moderate/Normal)",y="-log10(p value)")+
  theme(plot.title = element_text(color="black", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="black", face="bold", size=18))

## A_s.n
re[re$A.s.t.p  >= 0.05, "colgroup"] <- "black"
re[re$A.s.t.p  < 0.05, "colgroup"] <- "orange"
re[re$A.s.t.p  < 0.05 & re$A.s.fc >=fc, "colgroup"] <- "red" 
re[re$A.s.t.p  < 0.05 & re$A.s.fc <= 1/fc , "colgroup"] <- "green" 

  ggplot(data=re, aes(log2(A.s.fc), -log10(A.s.t.p)))+
  geom_point(size=2,color=re$colgroup)+
  geom_text_repel(aes(label=analytes), data=subset(re,colgroup!="black"), size=3, col='black')+
  geom_line(data=sigcut, aes(xx, yy),lty=3,lwd=1)+
  geom_line(data=sigcutfc, aes(xxp, yyp),lty=3,lwd=1)+
  geom_line(data=sigcutfc, aes(xxn, yyp),lty=3,lwd=1)+
  ggtitle("A: Severe vs Normal (p 0.05, fc 1.2)") +
  labs(x="log2(Severe/Normal)",y="-log10(p value)")+
  theme(plot.title = element_text(color="black", face="bold", size=18, hjust=0)) +
  theme(axis.title = element_text(color="black", face="bold", size=18))

  
  ## B_m.n
  
  re[re$B.m.t.p  >= 0.05, "colgroup"] <- "black"
  re[re$B.m.t.p  < 0.05, "colgroup"] <- "orange"
  re[re$B.m.t.p  < 0.05 & re$B.m.fc >=fc, "colgroup"] <- "red" 
  re[re$B.m.t.p  < 0.05 & re$B.m.fc <= 1/fc , "colgroup"] <- "green" 
  
  
  ggplot(data=re, aes(log2(B.m.fc), -log10(B.m.t.p)))+
    geom_point(size=2,color=re$colgroup)+
    geom_text_repel(aes(label=analytes), data=subset(re,colgroup!="black"), size=3, col='black')+
    geom_line(data=sigcut, aes(xx, yy),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxp, yyp),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxn, yyp),lty=3,lwd=1)+
    ggtitle("B: Moderate vs Normal (p 0.05, fc 1.2)") +
    labs(x="log2(Moderate/Normal)",y="-log10(p value)")+
    theme(plot.title = element_text(color="black", face="bold", size=18, hjust=0)) +
    theme(axis.title = element_text(color="black", face="bold", size=18))
  
  ## B_s.n
  re[re$B.s.t.p  >= 0.05, "colgroup"] <- "black"
  re[re$B.s.t.p  < 0.05, "colgroup"] <- "orange"
  re[re$B.s.t.p  < 0.05 & re$B.s.fc >=fc, "colgroup"] <- "red" 
  re[re$B.s.t.p  < 0.05 & re$B.s.fc <= 1/fc , "colgroup"] <- "green" 
  
  ggplot(data=re, aes(log2(B.s.fc), -log10(B.s.t.p)))+
    geom_point(size=2,color=re$colgroup)+
    geom_text_repel(aes(label=analytes), data=subset(re,colgroup!="black"), size=3, col='black')+
    geom_line(data=sigcut, aes(xx, yy),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxp, yyp),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxn, yyp),lty=3,lwd=1)+
    ggtitle("B: Severe vs Normal (p 0.05, fc 1.2)") +
    labs(x="log2(Severe/Normal)",y="-log10(p value)")+
    theme(plot.title = element_text(color="black", face="bold", size=18, hjust=0)) +
    theme(axis.title = element_text(color="black", face="bold", size=18))

  ## C_m.n
  
  re[re$C.m.t.p  >= 0.05, "colgroup"] <- "black"
  re[re$C.m.t.p  < 0.05, "colgroup"] <- "orange"
  re[re$C.m.t.p  < 0.05 & re$C.m.fc >=fc, "colgroup"] <- "red" 
  re[re$C.m.t.p  < 0.05 & re$C.m.fc <= 1/fc , "colgroup"] <- "green" 
  
  
  ggplot(data=re, aes(log2(C.m.fc), -log10(C.m.t.p)))+
    geom_point(size=2,color=re$colgroup)+
    geom_text_repel(aes(label=analytes), data=subset(re,colgroup!="black"), size=3, col='black')+
    geom_line(data=sigcut, aes(xx, yy),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxp, yyp),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxn, yyp),lty=3,lwd=1)+
    ggtitle("C: Moderate vs Normal (p 0.05, fc 1.2)") +
    labs(x="log2(Moderate/Normal)",y="-log10(p value)")+
    theme(plot.title = element_text(color="black", face="bold", size=18, hjust=0)) +
    theme(axis.title = element_text(color="black", face="bold", size=18))
  
  ## C_s.n
  re[re$C.s.t.p  >= 0.05, "colgroup"] <- "black"
  re[re$C.s.t.p  < 0.05, "colgroup"] <- "orange"
  re[re$C.s.t.p  < 0.05 & re$C.s.fc >=fc, "colgroup"] <- "red" 
  re[re$C.s.t.p  < 0.05 & re$C.s.fc <= 1/fc , "colgroup"] <- "green" 
  
  ggplot(data=re, aes(log2(C.s.fc), -log10(C.s.t.p)))+
    geom_point(size=2,color=re$colgroup)+
    geom_text_repel(aes(label=analytes), data=subset(re,colgroup!="black"), size=3, col='black')+
    geom_line(data=sigcut, aes(xx, yy),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxp, yyp),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxn, yyp),lty=3,lwd=1)+
    ggtitle("C: Severe vs Normal (p 0.05, fc 1.2)") +
    labs(x="log2(Severe/Normal)",y="-log10(p value)")+
    theme(plot.title = element_text(color="black", face="bold", size=18, hjust=0)) +
    theme(axis.title = element_text(color="black", face="bold", size=18))

  
  
  ## n_B.A
  
  re[re$N.b.t.p  >= 0.05, "colgroup"] <- "black"
  re[re$N.b.t.p  < 0.05, "colgroup"] <- "orange"
  re[re$N.b.t.p   < 0.05 & re$N.b.fc>=fc, "colgroup"] <- "red" 
  re[re$N.b.t.p   < 0.05 & re$N.b.fc <= 1/fc , "colgroup"] <- "green" 
  
  
  ggplot(data=re, aes(log2(N.b.fc), -log10(N.b.t.p)))+
    geom_point(size=2,color=re$colgroup)+
    geom_text_repel(aes(label=analytes), data=subset(re,colgroup!="black"), size=3, col='black')+
    geom_line(data=sigcut, aes(xx, yy),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxp, yyp),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxn, yyp),lty=3,lwd=1)+
    ggtitle("Normal: B vs A (p 0.05, fc 1.2)") +
    labs(x="log2(B/A)",y="-log10(p value)")+
    theme(plot.title = element_text(color="black", face="bold", size=18, hjust=0)) +
    theme(axis.title = element_text(color="black", face="bold", size=18))
  
  ## n_C.A
  re[re$N.c.t.p  >= 0.05, "colgroup"] <- "black"
  re[re$N.c.t.p  < 0.05, "colgroup"] <- "orange"
  re[re$N.c.t.p   < 0.05 & re$N.c.fc>=fc, "colgroup"] <- "red" 
  re[re$N.c.t.p   < 0.05 & re$N.c.fc <= 1/fc , "colgroup"] <- "green" 
  
  
  ggplot(data=re, aes(log2(N.c.fc), -log10(N.c.t.p)))+
    geom_point(size=2,color=re$colgroup)+
    geom_text_repel(aes(label=analytes), data=subset(re,colgroup!="black"), size=3, col='black')+
    geom_line(data=sigcut, aes(xx, yy),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxp, yyp),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxn, yyp),lty=3,lwd=1)+
    ggtitle("Normal: C vs A (p 0.05, fc 1.2)") +
    labs(x="log2(C/A)",y="-log10(p value)")+
    theme(plot.title = element_text(color="black", face="bold", size=18, hjust=0)) +
    theme(axis.title = element_text(color="black", face="bold", size=18))
  
  
  ## m_B.A
  
  re[re$M.b.t.p  >= 0.05, "colgroup"] <- "black"
  re[re$M.b.t.p  < 0.05, "colgroup"] <- "orange"
  re[re$M.b.t.p   < 0.05 & re$M.b.fc>=fc, "colgroup"] <- "red" 
  re[re$M.b.t.p   < 0.05 & re$M.b.fc <= 1/fc , "colgroup"] <- "green" 
  
  
  ggplot(data=re, aes(log2(M.b.fc), -log10(M.b.t.p)))+
    geom_point(size=2,color=re$colgroup)+
    geom_text_repel(aes(label=analytes), data=subset(re,colgroup!="black"), size=3, col='black')+
    geom_line(data=sigcut, aes(xx, yy),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxp, yyp),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxn, yyp),lty=3,lwd=1)+
    ggtitle("Moderate: B vs A (p 0.05, fc 1.2)") +
    labs(x="log2(B/A)",y="-log10(p value)")+
    theme(plot.title = element_text(color="black", face="bold", size=18, hjust=0)) +
    theme(axis.title = element_text(color="black", face="bold", size=18))
  
  ## m_C.A
  re[re$M.c.t.p  >= 0.05, "colgroup"] <- "black"
  re[re$M.c.t.p  < 0.05, "colgroup"] <- "orange"
  re[re$M.c.t.p   < 0.05 & re$M.c.fc>=fc, "colgroup"] <- "red" 
  re[re$M.c.t.p   < 0.05 & re$M.c.fc <= 1/fc , "colgroup"] <- "green" 
  
  
  ggplot(data=re, aes(log2(M.c.fc), -log10(M.c.t.p)))+
    geom_point(size=2,color=re$colgroup)+
    geom_text_repel(aes(label=analytes), data=subset(re,colgroup!="black"), size=3, col='black')+
    geom_line(data=sigcut, aes(xx, yy),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxp, yyp),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxn, yyp),lty=3,lwd=1)+
    ggtitle("Moderate: C vs A (p 0.05, fc 1.2)") +
    labs(x="log2(C/A)",y="-log10(p value)")+
    theme(plot.title = element_text(color="black", face="bold", size=18, hjust=0)) +
    theme(axis.title = element_text(color="black", face="bold", size=18))
  
  
  ## s_B.A
  
  re[re$S.b.t.p  >= 0.05, "colgroup"] <- "black"
  re[re$S.b.t.p  < 0.05, "colgroup"] <- "orange"
  re[re$S.b.t.p   < 0.05 & re$S.b.fc>=fc, "colgroup"] <- "red" 
  re[re$S.b.t.p   < 0.05 & re$S.b.fc <= 1/fc , "colgroup"] <- "green" 
  
  
  ggplot(data=re, aes(log2(S.b.fc), -log10(S.b.t.p)))+
    geom_point(size=2,color=re$colgroup)+
    geom_text_repel(aes(label=analytes), data=subset(re,colgroup!="black"), size=3, col='black')+
    geom_line(data=sigcut, aes(xx, yy),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxp, yyp),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxn, yyp),lty=3,lwd=1)+
    ggtitle("Severe: B vs A (p 0.05, fc 1.2)") +
    labs(x="log2(B/A)",y="-log10(p value)")+
    theme(plot.title = element_text(color="black", face="bold", size=18, hjust=0)) +
    theme(axis.title = element_text(color="black", face="bold", size=18))
  
  ## s_C.A
  re[re$S.c.t.p  >= 0.05, "colgroup"] <- "black"
  re[re$S.c.t.p  < 0.05, "colgroup"] <- "orange"
  re[re$S.c.t.p   < 0.05 & re$S.c.fc>=fc, "colgroup"] <- "red" 
  re[re$S.c.t.p   < 0.05 & re$S.c.fc <= 1/fc , "colgroup"] <- "green" 
  
  
  ggplot(data=re, aes(log2(S.c.fc), -log10(S.c.t.p)))+
    geom_point(size=2,color=re$colgroup)+
    geom_text_repel(aes(label=analytes), data=subset(re,colgroup!="black"), size=3, col='black')+
    geom_line(data=sigcut, aes(xx, yy),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxp, yyp),lty=3,lwd=1)+
    geom_line(data=sigcutfc, aes(xxn, yyp),lty=3,lwd=1)+
    ggtitle("Severe: C vs A (p 0.05, fc 1.2)") +
    labs(x="log2(C/A)",y="-log10(p value)")+
    theme(plot.title = element_text(color="black", face="bold", size=18, hjust=0)) +
    theme(axis.title = element_text(color="black", face="bold", size=18))
