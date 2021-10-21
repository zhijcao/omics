library(ggplot2)
library(plotROC)
library(OptimalCutpoints)
library(pROC)
library(ggrepel)
#dm<-sampleday1[,c("status2", "FGF23")]
#x=0.5
#y=0.25
#title="test"

#dm<-data.frame(d=1:100,m=1:100)
#dm<-within(dm, Rank <- as.factor(cut(m, quantile(m,probs=seq(0,1,0.33)), include.lowest=TRUE, labels=F)))

###############################
#######ROC 
###############################


ROC<-function(dm=dm, title="title",x=0.5,y=0.25){
  names(dm)<-c("d","m")
  p1<-ggplot(dm)+geom_roc(aes(d = d, m =m))
  auc<-round(calc_auc(p1)$AUC,2) 
  
  p<-ggplot(dm)+geom_roc(aes(d = d, m = m), n.cuts=10)+
    style_roc(theme = theme_grey,ylab="Sensitivity", xlab = "1 - Specificity")+
    annotate("text",x=x,y=y+0.2, label=paste("AUC :",auc))+
    ggtitle(title)
  p
}

###############################
#######Ranked ROC end
###############################

###############################
#######ROC marker decrease
###############################


ROC_decrease<-function(dm=dm, title="title",x=0.5,y=0.25){
  names(dm)<-c("d","m")
  p1<-ggplot(dm)+geom_roc(aes(d = d, m =-m))
  auc<-round(calc_auc(p1)$AUC,2) 
  
  p<-ggplot(dm)+geom_roc(aes(d = d, m = -m),n.cuts=10)+
    style_roc(theme = theme_grey,ylab="Sensitivity", xlab = "1 - Specificity")+
    annotate("text",x=x,y=y+0.2, label=paste("AUC :",auc))+
    ggtitle(title)
  p
}

###############################
#######Ranked ROC end
###############################

###############################
#######ROC optimalcutpoint auto directon 
###############################

dm<-A_moderate[, c("labPos","CXCL6")]
dm<-A_moderate[, c("labPos","CCL21")]

dm<-A_moderate[, c("labPos","CCL21")]
dm<-B_moderate[, c("labPos","CCL21")]
dm<-C_moderate[, c("labPos","CCL21")]


title="ROC"
n.cutoff=10
labelsize=3
method="Youden"

rocoptimalcutpoint<-function(dm=dm,title="ROC",n.cutoff=10,labelsize=3,method="Youden"){
  names(dm)<-c("d","m")   #d:lable 0 or 1; m: measurement
  #tem1<-roc(dm$d, dm$m, direction="auto")
  
  if (median(dm[dm$d==1, ]$m, na.rm=TRUE)>=median(dm[dm$d==0, ]$m, na.rm=TRUE))
    direction<-"<"   else direction<-">"
    
  tem<- optimal.cutpoints(X = "m", status = "d", tag.healthy = 0,
                                               methods = method, data = dm, pop.prev = NULL, categorical.cov = NULL, direction=direction,
                                               control = control.cutpoints(), ci.fit = FALSE, conf.level = 0.95, trace = FALSE)
  
  sespcut<-tem$Youden$Global$measures.acc
  optimal<-tem$Youden$Global$optimal.cutoff
  auc<-round(sespcut$AUC[[1]],2)
  optcutoff<-round(optimal$cutoff[1],1)
  optse<- round(optimal$Se[[1]],2)
  optsp<- round(optimal$Sp[[1]],2)
  optj<-round(tem$Youden$Global$optimal.criterion,2)
  # str(sespcut)
  # str(tem)
  if (direction==">")
  df_roc<-data.frame(fpr=c(0,1-as.numeric(sespcut$Sp)), sensitivity=c(0,as.numeric(sespcut$Se)), 
                     cutoff=c("inf", sespcut$cutoffs)) else df_roc<-data.frame(fpr=c(1-as.numeric(sespcut$Sp), 0), sensitivity=c(as.numeric(sespcut$Se), 0), 
                          cutoff=c(sespcut$cutoffs, "inf")) 
  # jindex=tem$Youden$Global$criterion)
  # df_j<-data.frame(x=c(1-optsp,1-optsp),y=c(1-optsp,optse))
  df_j<-data.frame(x=1-optsp,y=optse)
   direction1<-ifelse(direction==">","-","+")
   notes<-paste(" auc:",auc, " J:", optj, "\n", " cutoff:", optcutoff,  "\n"," sen:", optse, " spe:",optsp)
   
  p<-ggplot(data=df_roc, aes(fpr, sensitivity))+
    theme_bw()+
    geom_path(color="black",size=0.5)+
    #geom_abline(linetype="1F")+
    #geom_line(aes(x,y),linetype="dotted",size=0.5, data=df_j)+
    #geom_point(aes(x, y),data=df_j, size=1)+
    labs(x="1-Specificity",y="Sensitivity")+
    ggtitle(paste(title, direction1, sep=" "))+
    annotate("text",x=0.6, y=0.15, label=notes, size = 3, color="black")+
    theme(axis.title=element_text(size=12,color="black"),axis.text=element_text(size=10, color="black"))+
    theme(plot.title=element_text(size=12, color="black"))+
    #theme(legend.position="right")+
    theme(plot.margin=unit(c(1,1,1,1),"line"))+
    coord_cartesian(xlim=c(0, 1), ylim=c(0, 1))+
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
    theme(axis.ticks = element_line(color="black"))
  print(p)
  
  if (n.cutoff>0){
    n=nrow(df_roc)
    if (n.cutoff>n) n.cutoff=n
    cutoffat<-seq(1,n,ceiling(n/n.cutoff))
    maxj<-as.numeric(rownames(df_roc[df_roc$jindex==max(df_roc$jindex), ]))
    cutoffat<-unique(c(maxj,cutoffat))
    p+geom_text_repel(aes(fpr, sensitivity,label=cutoff), data = df_roc[cutoffat,], size=labelsize)

   }
}

###############################
#######ROC auto direction end
###############################


rocoptimalsesp<-function(dm=dm,analyte="name", method="Youden"){
  names(dm)<-c("d","m") #d:lable 0 or 1; m: measurement
  #tem1<-roc(dm$d, dm$m, direction="auto")
  
  if (median(dm[dm$d==1, ]$m, na.rm=TRUE)>=median(dm[dm$d==0, ]$m, na.rm=TRUE))
    direction<-"<"   else direction<-">"
    
    tem<- optimal.cutpoints(X = "m", status = "d", tag.healthy = 0,
                            methods = method, data = dm, pop.prev = NULL, categorical.cov = NULL, direction=direction,
                            control = control.cutpoints(), ci.fit = FALSE, conf.level = 0.95, trace = FALSE)
    
    sespcut<-tem$Youden$Global$measures.acc
    optimal<-tem$Youden$Global$optimal.cutoff
    auc<-round(sespcut$AUC[[1]],2)
    optcutoff<-round(optimal$cutoff[1],1)
    optse<- round(optimal$Se[[1]],2)
    optsp<- round(optimal$Sp[[1]],2)
    optj<-round(tem$Youden$Global$optimal.criterion,2)
   
    sesp<-data.frame(analyte=analyte, auc=auc,optse=optse, optsp=optsp, optj=optj)
    
    return (sesp)
}

###############################
#######ROC auto direction end
###############################








###############################
#######ROC auto directon 
###############################


rocauto<-function(dm=dm,title="ROC",n.cutoff=10,labelsize=3){
  names(dm)<-c("d","m")
  tem<-roc(dm$d, dm$m, direction="auto")
  df_roc<-data.frame(fpr=1-tem$specificities,sensitivity=tem$sensitivities,cutoff=tem$thresholds)
  auc<-round(tem$auc[1],2)
  direction<-ifelse(tem$direction==">","-","+")
  
  p<-ggplot(data=df_roc, aes(fpr, sensitivity))+
    theme_grey()+ 
    geom_point(size=1)+
    geom_path()+
    labs(x="1-Specificity",y="Sensitivity")+
    ggtitle(paste(title, direction, sep=": "))+
    annotate("text",x=0.75, y=0.1, label=paste("AUC :",auc), size = 3)+
    theme(plot.title = element_text(color="black", face="bold", size=10, hjust=0))+
    theme(axis.title = element_text(color="black", face="bold", size=10))+
    theme(plot.margin=unit(c(1,1,1,1),"line"))
  
  if (n.cutoff>0){
    n=nrow(df_roc)
    cutoffat<-seq(1,n,ceiling(n/n.cutoff))
    p+geom_text_repel(aes(label=round(cutoff,1)), data = df_roc[cutoffat,], size=labelsize)
  }
}

###############################
#######ROC auto directon end
###############################




###############################
#######Ranked ROC 
###############################


RankROC<-function(dm=dm, title="title",x=0.5,y=0.25){
names(dm)<-c("d","m")
dm[,"total"]<-"0-100th"
dm<-within(dm, Rank <- as.factor(cut(m, quantile(m), include.lowest=TRUE, labels=F)))
dm$Rank<-factor(dm$Rank, labels=c("0-25th","25-50th","50-75th","75-100th"))

p1<-ggplot(dm)+geom_roc(aes(d = d, m =m))
p2<-ggplot(dm)+geom_roc(aes(d = d, m = m,color=Rank),n.cuts = 0) 

rank_auc<-round(calc_auc(p2)$AUC,2)
total_auc<-round(calc_auc(p1)$AUC,2) 

p<-ggplot(dm)+geom_roc(aes(d = d, m = m,color=Rank),n.cuts = 0)+
  geom_roc(aes(d = d, m = m,shape=total),cutoffs.at=quantile(dm$m))+
  style_roc(theme = theme_grey,ylab="Sensitivity", xlab = "1 - Specificity")+
  annotate("text",x=x,y=y, label=paste("0-25th AUC :",rank_auc[1]))+
  annotate("text",x=x,y=y+0.05, label=paste("25-50th AUC :",rank_auc[2]))+
  annotate("text",x=x,y=y+0.1, label=paste("50-75th AUC :",rank_auc[3]))+
  annotate("text",x=x,y=y+0.15, label=paste("75-100th AUC :",rank_auc[4]))+
  annotate("text",x=x,y=y+0.2, label=paste("0-100th AUC :",total_auc))+
  ggtitle(title)
p
}

###############################
#######Ranked with tertile ROC end 
###############################
#m<-1:100
#data.frame(m=m,rank=Rank)

RankROC<-function(dm=dm, title="title",x=0.5,y=0.25){
  names(dm)<-c("d","m")
  dm[,"total"]<-"0-100th"
  dm<-within(dm, Rank <- as.factor(cut(m, quantile(m,probs = seq(0, 1, 1/3),names=F), include.lowest=TRUE, labels=F)))
  dm$Rank<-factor(dm$Rank, labels=c("1st tertile","2nd tertile","3rd tertile"))
  
  p1<-ggplot(dm)+geom_roc(aes(d = d, m =m))
  p2<-ggplot(dm)+geom_roc(aes(d = d, m = m,color=Rank),n.cuts = 0) 
  
  rank_auc<-round(calc_auc(p2)$AUC,2)
  total_auc<-round(calc_auc(p1)$AUC,2) 
  
  p<-ggplot(dm)+geom_roc(aes(d = d, m = m,color=Rank),n.cuts = 0)+
    #geom_roc(aes(d = d, m = m,shape=total),cutoffs.at=quantile(dm$m))+
    style_roc(theme = theme_grey,ylab="Sensitivity", xlab = "1 - Specificity")+
    annotate("text",x=x,y=y, label=paste("1st tertile AUC :",rank_auc[1]))+
    annotate("text",x=x,y=y+0.05, label=paste("2nd tertile AUC :",rank_auc[2]))+
    annotate("text",x=x,y=y+0.1, label=paste("3rd tertile AUC :",rank_auc[3]))+
    #annotate("text",x=x,y=y+0.15, label=paste("0-100th AUC :",total_auc))+
    ggtitle(title)
  p
}

###############################
#######Ranked ROC end
###############################





