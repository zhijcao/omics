library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)

rm(list=ls())

qcid<- read.csv("C:/zhijuncao/metabolomics/lcmsre/serum_pos_skyline_result.csv")
qcid<- read.csv("C:/zhijuncao/metabolomics/lcmsre/serum_neg_skyline_result.csv")

qcid<- read.csv("C:/zhijuncao/metabolomics/lcmsre/urine_pos_skyline_result.csv")
qcid<- read.csv("C:/zhijuncao/metabolomics/lcmsre/urine_neg_skyline_result.csv")

qcid<- read.csv("C:/zhijuncao/metabolomics/lcmsre/nist_pos_skyline_result.csv")
qcid<- read.csv("C:/zhijuncao/metabolomics/lcmsre/nist_neg_skyline_result.csv")


qcid<- read.csv("C:/zhijuncao/metabolomics/lcmsre/serum_neg_feature_table_long.csv")
qcid<- read.csv("C:/zhijuncao/metabolomics/lcmsre/serum_pos_feature_table_long.csv")

qcid<- read.csv("C:/zhijuncao/metabolomics/lcmsre/nist_pos_feature_table_long.csv")
qcid<- read.csv("C:/zhijuncao/metabolomics/lcmsre/nist_pos_FT_skyline_result.csv")

title <- "nist_pos"
title <- "nist_neg"
title <- "serum_neg_FT"
title <- "serum_pos_FT"
title <- "nist_pos_FT"

names(qcid)
str(qcid)
dim(qcid)
CV <- function(x) sd(x, na.rm=TRUE)/mean(x,na.rm=TRUE)*100


area_sample_sum <- qcid %>% group_by(Replicate) %>% summarise_at("Total.Area", funs(sum, mean, median), na.rm = TRUE)

qcid_sum <- left_join(qcid, area_sample_sum, by="Replicate")
qcid_normalization <- qcid_sum %>% mutate(toSum=Total.Area/sum,toMean=Total.Area/mean,toMedian=Total.Area/median)


#####cross all samples cv
qcid_normalization_cv <- qcid_normalization %>% group_by(metaID) %>% 
  summarise_at(c("Total.Area","toSum", "toMean","toMedian"),CV) %>% 
  gather(key="Normalization", value="CV", -metaID)
qcid_normalization_cv1 <- qcid_normalization %>% group_by(metaID) %>% 
  summarise_at(c("Total.Area","toSum", "toMean","toMedian"),CV)
qcid_normalization_cv$Normalization <- factor(qcid_normalization_cv$Normalization, levels=c("Total.Area", "toSum", "toMean", "toMedian"))

write.csv(qcid_normalization_cv1, paste( "C:/zhijuncao/metabolomics/lcmsre/", title, "_skyline_cv.csv", sep=""))

#####by project cv
qcid_normalization_project_cv <- qcid_normalization %>% group_by(metaID, Project) %>% 
  summarise_at(c("Total.Area","toSum", "toMean","toMedian"),CV) %>% 
  gather(key="Normalization", value="CV",-metaID, -Project)
qcid_normalization_project_cv1 <- qcid_normalization %>% group_by(metaID, Project) %>% 
  summarise_at(c("Total.Area","toSum", "toMean","toMedian"),CV)
qcid_normalization_project_cv$Normalization <- factor(qcid_normalization_project_cv$Normalization, levels=c("Total.Area", "toSum", "toMean", "toMedian"))


write.csv(qcid_normalization_project_cv1, paste( "C:/zhijuncao/metabolomics/lcmsre/", title, "_skyline_project_cv.csv", sep=""))

##########################result plot###############################

pdf(paste("C:/zhijuncao/metabolomics/lcmsre/",title,  "_skyline_cv_boxplot.pdf", sep=""), width = 8, height = 6)
ggplot(qcid_normalization_cv, aes(Normalization, CV))+
  theme_classic()+
  geom_boxplot()+
  labs(x="")+
  theme(text=element_text(size=12, colour = "black"))+
  theme(axis.text =element_text(color = "black"))+
  theme(axis.text.x =element_text(angle=0,hjust=0.5))

ggplot(qcid_normalization_project_cv, aes(Normalization, CV, fill=Project))+
  theme_classic()+
  geom_boxplot()+
  labs(x="")+
  theme(text=element_text(size=12, colour = "black"))+
  theme(axis.text =element_text(color = "black"))+
  theme(axis.text.x =element_text(angle=0,hjust=0.5))

dev.off()



########xcms feature plot
pdf(paste("C:/zhijuncao/metabolomics/lcmsre/",title,  "_skyline_area_RT_line.pdf", sep=""), width = 8, height = 6)
ggplot(qcid, aes(Replicate, rtmed, fill=metaID))+
  theme_classic()+
  geom_crossbar(aes(ymin = rtmin, ymax =rtmax), width = 0.4, show.legend=F, linetype="blank")+
  labs(x="", y="Retention Time")+
  theme(text=element_text(size=10, colour = "black"))+
  theme(axis.text =element_text(color = "black"))+
  theme(axis.text.x =element_text(angle=45,hjust=1))

ggplot(qcid_normalization, aes(Replicate,Total.Area, color=Project, group=metaID))+
  theme_classic()+
  geom_line()+
  labs(x="")+
  theme(text=element_text(size=10, colour = "black"))+
  theme(axis.text =element_text(color = "black"))+
  theme(axis.text.x =element_text(angle=45,hjust=1))

ggplot(qcid_normalization, aes(Replicate,toSum, color=Project, group=metaID))+
  theme_classic()+
  geom_line()+
  labs(x="")+
  theme(text=element_text(size=10, colour = "black"))+
  theme(axis.text =element_text(color = "black"))+
  theme(axis.text.x =element_text(angle=45,hjust=1))

ggplot(qcid_normalization, aes(Replicate,toMean, color=Project, group=metaID))+
  theme_classic()+
  geom_line()+
  labs(x="")+
  theme(text=element_text(size=10, colour = "black"))+
  theme(axis.text =element_text(color = "black"))+
  theme(axis.text.x =element_text(angle=45,hjust=1))

ggplot(qcid_normalization, aes(Replicate,toMedian, color=Project, group=metaID))+
  theme_classic()+
  geom_line()+
  labs(x="")+
  theme(text=element_text(size=10, colour = "black"))+
  theme(axis.text =element_text(color = "black"))+
  theme(axis.text.x =element_text(angle=45,hjust=1))

dev.off()





########
pdf(paste("C:/zhijuncao/metabolomics/lcmsre/",title,  "_skyline_area_RT_line.pdf", sep=""), width = 8, height = 6)
ggplot(qcid, aes(Replicate, Best.Retention.Time, fill=metaID))+
  theme_classic()+
  geom_crossbar(aes(ymin = Min.Start.Time, ymax =Max.End.Time), width = 0.4, show.legend=F, linetype="blank")+
  labs(x="", y="Retention Time")+
  theme(text=element_text(size=10, colour = "black"))+
  theme(axis.text =element_text(color = "black"))+
  theme(axis.text.x =element_text(angle=45,hjust=1))

ggplot(qcid_normalization, aes(Replicate,Total.Area, color=Project, group=metaID))+
  theme_classic()+
  geom_line()+
  labs(x="")+
  theme(text=element_text(size=10, colour = "black"))+
  theme(axis.text =element_text(color = "black"))+
  theme(axis.text.x =element_text(angle=45,hjust=1))

ggplot(qcid_normalization, aes(Replicate,toSum, color=Project, group=metaID))+
  theme_classic()+
  geom_line()+
  labs(x="")+
  theme(text=element_text(size=10, colour = "black"))+
  theme(axis.text =element_text(color = "black"))+
  theme(axis.text.x =element_text(angle=45,hjust=1))

ggplot(qcid_normalization, aes(Replicate,toMean, color=Project, group=metaID))+
  theme_classic()+
  geom_line()+
  labs(x="")+
  theme(text=element_text(size=10, colour = "black"))+
  theme(axis.text =element_text(color = "black"))+
  theme(axis.text.x =element_text(angle=45,hjust=1))

ggplot(qcid_normalization, aes(Replicate,toMedian, color=Project, group=metaID))+
  theme_classic()+
  geom_line()+
  labs(x="")+
  theme(text=element_text(size=10, colour = "black"))+
  theme(axis.text =element_text(color = "black"))+
  theme(axis.text.x =element_text(angle=45,hjust=1))

dev.off()






