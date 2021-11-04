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
library(gplots)
library(stringr)

# source("https://bioconductor.org/biocLite.R")
# biocLite("hmdbQuery")

library(hmdbQuery)
browseVignettes("hmdbQuery")
lk1 = HmdbEntry(prefix = "http://www.hmdb.ca/metabolites/", 
                id = "HMDB0000001")

library(multcomp)
library(nlme)
library(lme4)
library(MASS)
library(broom)
rm(list=ls())
fm1 <- lme(distance ~ age, data = Orthodont) 
Orthodont.df <- as.data.frame(Orthodont)
fm2 <- lme(distance ~ factor(age), data = Orthodont.df) 

summary(fm1)

x <- "alpha.tocopherol"
lm_fit<-function(x,time1=twodrug_TIME1, time2=twodrug_TIME2){
  
   fml <- paste("log(", x, ", 2) ~ TREATMENT*DIET")
  fml <-  as.formula(fml)
  re_lm1<-lm(fml,   data=time1)
  re_lm2<-lm(fml,   data=time2)
  t1_df <- tidy(re_lm1) %>% mutate(Time="T1")
  t2_df <- tidy(re_lm2) %>% mutate(Time="T2")
  t1t2 <- bind_rows(t1_df, t2_df) %>% add_column(Meta=x, .before = 1)
  return (t1t2)
}

lm_result <- lapply(names(twodrug_TIME1)[-c(1:15)], function(x)lm_fit(x))

lm_result_bind <- bind_rows(lm_result)
x <- c(0.01,0.04,1,5,4)
sig01 <- function(x) ifelse(x<0.05, 1, 0)

lm_result_bind_p <- lm_result_bind %>% spread_n(key = term, values = c(estimate, std.error, statistic, p.value))
lm_result_bind_p_sig <- lm_result_bind_p %>% dplyr::select(contains("p.value")) %>% mutate_all(sig01) %>% unite(sig, sep="")
lm_result_bind_p <- add_column(lm_result_bind_p, sig=lm_result_bind_p_sig$sig, .before=3) 
                                                                
lm_result_bind_p_only <- lm_result_bind_p %>% dplyr::select(Meta, Time, sig, contains("p.value"))
dim(lm_result_bind_p05)
names(lm_result_bind_p)


write.xlsx(lm_result_bind_p, "C:/zhijuncao/sun/diet-drug/lm_result_wid.xlsx")
write.xlsx(lm_result_bind_p_only, "C:/zhijuncao/sun/diet-drug/lm_result_wid_p.xlsx")




#sample <- read.csv("C:/zhijuncao/sun/diet-drug/scaledata_1.csv")
#sample <- read.csv("C:/zhijuncao/sun/diet-drug/twodrug_scaledata_1_remove2samples.csv")
sample <- read.csv("C:/zhijuncao/sun/diet-drug/twodrug_scaledata_1_remove2samples_568metabolites.csv")
sample <- read.csv("C:/zhijuncao/sun/diet-drug/twodrug_scaledata_1_remove2samples_565metabolites.csv")
metaid <- read.csv("C:/zhijuncao/sun/diet-drug/metaID_1.csv")

twodrug <- sample %>% filter(TREATMENT=="Ctrl"|TREATMENT=="Tgt"|TREATMENT=="Tmx") %>% droplevels()
twodrug_TIME1 <- twodrug %>% filter(TIME2=="T1") %>% droplevels()
twodrug_TIME2 <- twodrug %>% filter(TIME2=="T2") %>% droplevels()


write.csv(twodrug, "C:/zhijuncao/sun/diet-drug/twodrug_scaledata_1.csv")
twodrug <- read.csv( "C:/zhijuncao/sun/diet-drug/twodrug_scaledata_1.csv")
names(twodrug)[1:30]

test = sample[, c(5,6, 9, 16:19 )]
names(test)
names(sample)[1:20]

##pathway enrichment analysis 
#twodrug_re <- read.xlsx("C:/zhijuncao/sun/diet-drug/112018/two_Drug_all_p_fc_fdr_sig_01label.xlsx", sheet="meta568")
twodrug_re <- read.xlsx("C:/zhijuncao/sun/diet-drug/01282019/two_Drug_all_p_fc_fdr_sig_01label_020619.xlsx", sheet = "meta565")
twodrug_re <- read.xlsx("C:/zhijuncao/sun/diet-drug/01282019/two_drug_all_p_fc_fdr_01lable_02132019.xlsx")
dim(twodrug_re)
sample_names<- c("Ctrl_T1_CD_HFD","Ctrl_T2_CD_HFD", "Tgt_T1_CD_HFD","Tgt_T2_CD_HFD",
                 "Tmx_T1_CD_HFD","Tmx_T2_CD_HFD","CD_T1_Ctrl_Tgt","CD_T2_Ctrl_Tgt", 
                 "HFD_T1_Ctrl_Tgt", "HFD_T2_Ctrl_Tgt","CD_T1_Ctrl_Tmx","CD_T2_Ctrl_Tmx",
                 "HFD_T1_Ctrl_Tmx","HFD_T2_Ctrl_Tmx")
names(twodrug_re)
dim(twodrug_re)
pathway_enrich_re <- list()
for (i in 1:14){
  pathway_enrich_re[[i]]<- pathway_enrich(i=i)
}

pathway_enrich_df <- do.call(cbind, pathway_enrich_re) %>% data.frame()
names(pathway_enrich_df)<- sample_names
#write.csv(pathway_enrich_df, "C:/zhijuncao/sun/diet-drug/01282019/two_drug_sub_pathway_enrich_pvalue565.csv")
pathway_enrich_df <- read.csv( "C:/zhijuncao/sun/diet-drug/01282019/two_drug_sub_pathway_enrich_pvalue565.csv")
names(pathway_enrich_df) <- c( "X","Ctrl_T1:HFD vs CD","Ctrl_T2:HFD vs CD","Tgt_T1:HFD vs CD","Tgt_T2:HFD vs CD","Tmx_T1:HFD vs CD",  
                               "Tmx_T2:HFD vs CD", "CD_T1:Tgt vs Ctrl", "CD_T2:Tgt vs Ctrl","HFD_T1:Tgt vs Ctrl", "HFD_T2:Tgt vs Ctrl", "CD_T1:Tmx vs Ctrl", 
                               "CD_T2:Tmx vs Ctrl",  "HFD_T1:Tmx vs Ctrl", "HFD_T2:Tmx vs Ctrl")
sample_names <-  c("Ctrl_T1:HFD vs CD","Ctrl_T2:HFD vs CD","Tgt_T1:HFD vs CD","Tgt_T2:HFD vs CD","Tmx_T1:HFD vs CD",  
                    "Tmx_T2:HFD vs CD", "CD_T1:Tgt vs Ctrl", "CD_T2:Tgt vs Ctrl","HFD_T1:Tgt vs Ctrl", "HFD_T2:Tgt vs Ctrl", "CD_T1:Tmx vs Ctrl", 
                    "CD_T2:Tmx vs Ctrl",  "HFD_T1:Tmx vs Ctrl", "HFD_T2:Tmx vs Ctrl")
sup_pathway_enrich_re <- list()
for (i in 1:14){
  sup_pathway_enrich_re[[i]]<- sup_pathway_enrich(i=i)
}

sup_pathway_enrich_df <- do.call(cbind, sup_pathway_enrich_re) %>% data.frame()
names(sup_pathway_enrich_df)<- sample_names
#write.csv(sup_pathway_enrich_df, "C:/zhijuncao/sun/diet-drug/01282019/two_drug_sup_pathway_enrich_pvalue565.csv")
sup_pathway_enrich_df <- read.csv( "C:/zhijuncao/sun/diet-drug/01282019/two_drug_sup_pathway_enrich_pvalue565.csv")
names(sup_pathway_enrich_df) <- c( "X","Ctrl_T1:HFD vs CD","Ctrl_T2:HFD vs CD","Tgt_T1:HFD vs CD","Tgt_T2:HFD vs CD","Tmx_T1:HFD vs CD",  
                                   "Tmx_T2:HFD vs CD", "CD_T1:Tgt vs Ctrl", "CD_T2:Tgt vs Ctrl","HFD_T1:Tgt vs Ctrl", "HFD_T2:Tgt vs Ctrl", "CD_T1:Tmx vs Ctrl", 
                                   "CD_T2:Tmx vs Ctrl",  "HFD_T1:Tmx vs Ctrl", "HFD_T2:Tmx vs Ctrl")
unique(twodrug_re$SUB_PATHWAY)
unique(twodrug_re$SUPER_PATHWAY)
str(pathway_enrich_df)
sigmore1 <- rowSums(pathway_enrich_df[,-1]<0.05)
sig_subpathway <- paste(pathway_enrich_df$X[sigmore1>=1])

sup_sigmore1 <- rowSums(sup_pathway_enrich_df[,-1]<0.15)

pdf("sub_pathway enrichment analysis heat map03122019_565_p05.pdf", width = 10, height = 12)
heatmap.2(as.matrix(pathway_enrich_df[sigmore1>=1,-1]),col=redblue, Rowv = TRUE, Colv = TRUE, margins = c(30, 15),
          labCol=names(pathway_enrich_df)[-1],labRow =pathway_enrich_df$X[sigmore1>=1] ,trace="none", 
          na.color = "grey",srtCol=45,cexRow = 0.8,keysize=0.8, key.xlab = "P value",  lhei=c(1,8),
          density.info="none",breaks = seq(0,0.4,0.05))
dev.off()

pdf("sup_pathway enrichment analysis heat map03122019_565.pdf", width = 10, height = 12)
heatmap.2(as.matrix(sup_pathway_enrich_df[sup_sigmore1>=1,-1]),col=redblue, Rowv = TRUE, Colv = TRUE, margins = c(55, 15),
          labCol=names(sup_pathway_enrich_df)[-1],labRow = sup_pathway_enrich_df$X[sup_sigmore1>=1],trace="none", 
          na.color = "grey",srtCol=45, cexRow = 0.8,keysize=0.8, key.xlab = "P value",lhei=c(1,8),
          density.info="none",breaks = seq(0, 0.4, 0.05))
dev.off()

pathway_heatmap <- function(df=sup_pathway_enrich_df,name=sample_names){
  heatmap.2(as.matrix(df[,name]),col=redblue, Rowv = TRUE, Colv = TRUE, margins = c(8, 15+14-min(14,length(name))),
          labCol=name,labRow = df$X,trace="none", #lhei = c(0.5,4),
          na.color = "grey",srtCol=45, cexRow = 0.8, cexCol = 1, key.xlab = "P value", keysize=0.8,
          density.info="none",breaks = seq(0, 0.3, 0.05))
}
sup_pathway_heatmap <- function(df=sup_pathway_enrich_df,name=sample_names){
  heatmap.2(as.matrix(df[,name]),col=redblue, Rowv = TRUE, Colv = TRUE, margins = c(25, 15+14-min(14,length(name))),
            labCol=name,labRow = df$X,trace="none", #lhei = c(0.5,4),
            na.color = "grey",srtCol=45, cexRow = 0.8, cexCol = 1, key.xlab = "P value", keysize=0.8,
            density.info="none",breaks = seq(0, 0.3, 0.05))
}

names(sup_pathway_enrich_df)

t1 <- c("Ctrl_T1:HFD vs CD","Tgt_T1:HFD vs CD","Tmx_T1:HFD vs CD",  
 "CD_T1:Tgt vs Ctrl", "HFD_T1:Tgt vs Ctrl",  "CD_T1:Tmx vs Ctrl", "HFD_T1:Tmx vs Ctrl")
  
            
t2 <- c("Ctrl_T2:HFD vs CD","Tgt_T2:HFD vs CD",  
        "Tmx_T2:HFD vs CD",  "CD_T2:Tgt vs Ctrl", "HFD_T2:Tgt vs Ctrl",  
        "CD_T2:Tmx vs Ctrl",   "HFD_T2:Tmx vs Ctrl")
           
t1_tgt <- c("Ctrl_T1:HFD vs CD","Tgt_T1:HFD vs CD", "CD_T1:Tgt vs Ctrl")

t2_tgt <- c("Ctrl_T2:HFD vs CD","Tgt_T2:HFD vs CD", "CD_T2:Tgt vs Ctrl")

t1_tmx <- c("Ctrl_T1:HFD vs CD","Tmx_T1:HFD vs CD",  "CD_T1:Tmx vs Ctrl", "HFD_T1:Tmx vs Ctrl")

t2_tmx <- c("Ctrl_T2:HFD vs CD","Tmx_T2:HFD vs CD",  "CD_T2:Tmx vs Ctrl", "HFD_T2:Tmx vs Ctrl")
       
pdf("sup_pathway enrichment analysis heat map020619_565_groups_10x8a.pdf", width = 10, height = 8)
sup_pathway_heatmap(df=sup_pathway_enrich_df,name=sample_names)
sup_pathway_heatmap(df=sup_pathway_enrich_df,name=t1)
sup_pathway_heatmap(df=sup_pathway_enrich_df,name=t2)
sup_pathway_heatmap(df=sup_pathway_enrich_df,name=t1_tgt)
sup_pathway_heatmap(df=sup_pathway_enrich_df,name=t2_tgt)
sup_pathway_heatmap(df=sup_pathway_enrich_df,name=t1_tmx)
sup_pathway_heatmap(df=sup_pathway_enrich_df,name=t2_tmx)
dev.off()

pdf("sub_pathway enrichment analysis heat map020619_565_groupsa.pdf", width = 10, height = 12)
pathway_heatmap(df=pathway_enrich_df,name=sample_names)
pathway_heatmap(df=pathway_enrich_df,name=t1)
pathway_heatmap(df=pathway_enrich_df,name=t2)
pathway_heatmap(df=pathway_enrich_df,name=t1_tgt)
pathway_heatmap(df=pathway_enrich_df,name=t2_tgt)
pathway_heatmap(df=pathway_enrich_df,name=t1_tmx)
pathway_heatmap(df=pathway_enrich_df,name=t2_tmx)
dev.off()




#####pathway analysis function
dim(twodrug_re)
pathway_enrich <- function(i=1){
#i <- 14
pathway <- table(ref=twodrug_re$SUB_PATHWAY, experiment=substr(twodrug_re$sig,i,i))
total <- table(substr(twodrug_re$sig,i,i))
pathway_change <- cbind(data.frame(nochange=pathway[,1]), data.frame(change=pathway[,2]))

notinpathway<- apply(pathway_change,1, function(x)(total-x)) %>% data.frame() #should add comments

df <- cbind(pathway_change,data.frame(t(notinpathway)))
fisher <- c()
#name <- "Ascorbate and Aldarate Metabolism"
for (name in rownames(df)){
  
 tem_matrix<- matrix(unlist(df[name,]),2,2, byrow = TRUE)

fisher[name] <- fisher.test(tem_matrix, alternative = "less")$p.value

}
return (fisher)
}

name <- "Alanine and Aspartate Metabolism" 
###################
names(twodrug_re)[1:10]
#####super pathway analysis function
sup_pathway_enrich <- function(i=1){
  pathway <- table(ref=twodrug_re$SUPER_PATHWA, experiment=substr(twodrug_re$sig,i,i))
  total <- table(substr(twodrug_re$sig,i,i))
  pathway_change <- cbind(data.frame(nochange=pathway[,1]), data.frame(change=pathway[,2]))
  
  notinpathway<- apply(pathway_change,1, function(x)(total-x)) %>% data.frame()
  
  df <- cbind(pathway_change,data.frame(t(notinpathway)))
  fisher <- c()
  for (name in rownames(df)){
    
    tem_matrix<- matrix(unlist(df[name,]),2,2, byrow = TRUE)
    
    fisher[name] <- fisher.test(tem_matrix, alternative = "less")$p.value
    
  }
  return (fisher)
}

###################



#twodrug <- read.csv("C:/zhijuncao/sun/diet-drug/twodrug_scaledata_1_remove2samples.csv")

#sample <- read.csv("C:/zhijuncao/sun/diet-drug/twodrug_scaledata_1_remove2samples_568metabolites.csv")
sample <- read.csv("C:/zhijuncao/sun/diet-drug/twodrug_scaledata_1_remove2samples_565metabolites.csv")
metaid <- read.csv("C:/zhijuncao/sun/diet-drug/metaID_565.csv")

sample <- twodrug

cc <- colorRamps::primary.colors(12)
sc <- colorRamps::primary.colors(8)
names(cc) <- levels(twodrug$SampleID) #column color
names(sc) <- levels(metaid$SUPER_PATHWAY) #row color
names(twodrug)[-c(1:15)]
metaid$SUPER_PATHWAY
metaid$SUB_PATHWAY
twodrug$SampleID

tem <- as.matrix(log(twodrug[,-c(1:15)],2))
pdf("heatmap565_log2.pdf", width = 12, height =8)
heatmap.2(t(tem),col=bluered, Rowv = TRUE, Colv = TRUE, margins = c(8, 22),
          labCol=twodrug$SampleID,labRow = FALSE,trace="none", RowSideColors=sc[metaid$SUPER_PATHWAY], 
          ColSideColors=cc[twodrug$SampleID],na.color = "grey", srtCol=45,cexRow = 0.8,keysize=0.8,
          density.info="none",breaks = seq(-2,2,0.2))
legend("topright", legend = names(cc)[1:12],fill = cc,title = "sample", cex = 1)
legend("bottomright", legend = names(sc)[1:8],fill = sc,title = "SUPER_PATHWAY", cex = 1)
dev.off()

#twodrug_re <- read.xlsx("C:/zhijuncao/sun/diet-drug/112018/two_Drug_all_p_fc_fdr_sig_01label.xlsx") 
twodrug_re <- read.xlsx("C:/zhijuncao/sun/diet-drug/01282019/two_Drug_all_p_fc_fdr_sig_01label_020619.xlsx", sheet = "meta565")
dim(twodrug_re)
names(twodrug_re)
write.csv(as.data.frame(table(twodrug_re$SUPER_PATHWAY)), "C:/zhijuncao/sun/diet-drug/01282019/SUPER_PATHWAY_metabolite565.csv")
write.csv(as.data.frame(table(twodrug_re$SUB_PATHWAY)), "C:/zhijuncao/sun/diet-drug/01282019/sub_PATHWAY_metabolite565.csv")

names(twodrug_re)[1:20]
name_order <- c("Name",  "PATHWAY_SORTORDER", "BIOCHEMICAL", "SUPER_PATHWAY",
                 "SUB_PATHWAY", "COMP_ID", "PLATFORM",  "CHEMICAL_ID",  "RI",                 
                 "MASS", "CAS", "PUBCHEM", "CHEMSPIDER",  "KEGG",     
                 "HMDB_ID", "sig",  "CD_HFD", "Tgt", "Tmx", 
                 "Ctrl_T1_CD_HFD.fc","Ctrl_T2_CD_HFD.fc", "Tgt_T1_CD_HFD.fc","Tgt_T2_CD_HFD.fc",
                 "Tmx_T1_CD_HFD.fc","Tmx_T2_CD_HFD.fc","CD_T1_Ctrl_Tgt.fc","CD_T2_Ctrl_Tgt.fc", 
                 "HFD_T1_Ctrl_Tgt.fc", "HFD_T2_Ctrl_Tgt.fc","CD_T1_Ctrl_Tmx.fc","CD_T2_Ctrl_Tmx.fc",
                 "HFD_T1_Ctrl_Tmx.fc","HFD_T2_Ctrl_Tmx.fc",
                 
                 "Ctrl_T1_CD_HFD.p", "Ctrl_T2_CD_HFD.p",  "Tgt_T1_CD_HFD.p",  "Tgt_T2_CD_HFD.p","Tmx_T1_CD_HFD.p","Tmx_T2_CD_HFD.p",
                 "CD_T1_Ctrl_Tgt.p", "CD_T2_Ctrl_Tgt.p","HFD_T1_Ctrl_Tgt.p","HFD_T2_Ctrl_Tgt.p","CD_T1_Ctrl_Tmx.p",
                  "CD_T2_Ctrl_Tmx.p", "HFD_T1_Ctrl_Tmx.p","HFD_T2_Ctrl_Tmx.p",
                 
                 "Ctrl_T1_CD_HFD.FDR","Ctrl_T2_CD_HFD.FDR", "Tgt_T1_CD_HFD.FDR", "Tgt_T2_CD_HFD.FDR",
                  "Tmx_T1_CD_HFD.FDR",  "Tmx_T2_CD_HFD.FDR", "CD_T1_Ctrl_Tgt.FDR",       
                   "CD_T2_Ctrl_Tgt.FDR", "HFD_T1_Ctrl_Tgt.FDR", "HFD_T2_Ctrl_Tgt.FDR","CD_T1_Ctrl_Tmx.FDR",             
                      "CD_T2_Ctrl_Tmx.FDR",  "HFD_T1_Ctrl_Tmx.FDR","HFD_T2_Ctrl_Tmx.FDR")          
fc <- c("Ctrl_T1_CD_HFD.fc","Ctrl_T2_CD_HFD.fc", "Tgt_T1_CD_HFD.fc","Tgt_T2_CD_HFD.fc",
        "Tmx_T1_CD_HFD.fc","Tmx_T2_CD_HFD.fc","CD_T1_Ctrl_Tgt.fc","CD_T2_Ctrl_Tgt.fc", 
        "HFD_T1_Ctrl_Tgt.fc", "HFD_T2_Ctrl_Tgt.fc","CD_T1_Ctrl_Tmx.fc","CD_T2_Ctrl_Tmx.fc",
        "HFD_T1_Ctrl_Tmx.fc","HFD_T2_Ctrl_Tmx.fc")                  
                       
twodrug_fc <- twodrug_re[,name_order]

names(twodrug_fc)
tem <- twodrug_re[,fc]
cc <- colorRamps::primary.colors(14)
sc <- colorRamps::primary.colors(8)
names(cc) <- names(tem) #column color
names(sc) <- levels(metaid$SUPER_PATHWAY) #row color


twodrug_fc <- read.xlsx("C:/zhijuncao/sun/diet-drug/01282019/two_drug_all_p_fc_fdr_01lable_02132019.xlsx", sheet="fc") 
names(twodrug_fc)[c(27:40)] <- sample_names

twodrug <- read.xlsx("C:/zhijuncao/sun/diet-drug/01282019/two_drug_all_p_fc_fdr_01lable_02132019.xlsx", sheet="Sheet 1") 
names(twodrug)
names_order <- c("Name",  "PATHWAY_SORTORDER", "BIOCHEMICAL", "SUPER_PATHWAY",
                "SUB_PATHWAY", "COMP_ID", "PLATFORM",  "CHEMICAL_ID",  "RI",                 
                "MASS", "CAS", "PUBCHEM", "CHEMSPIDER",  "KEGG",     
                "HMDB_ID", "sig",  "CD_HFD","CD_HFD_T1", "CD_HFD_T2", "Tgt", "Tmx", "Tgt_T1","Tgt_T2", "Tmx_T1","Tmx_T2",
                "Ctrl_T1_CD_HFD.fc","Ctrl_T2_CD_HFD.fc", "Tgt_T1_CD_HFD.fc","Tgt_T2_CD_HFD.fc",
                "Tmx_T1_CD_HFD.fc","Tmx_T2_CD_HFD.fc","CD_T1_Ctrl_Tgt.fc","CD_T2_Ctrl_Tgt.fc", 
                "HFD_T1_Ctrl_Tgt.fc", "HFD_T2_Ctrl_Tgt.fc","CD_T1_Ctrl_Tmx.fc","CD_T2_Ctrl_Tmx.fc",
                "HFD_T1_Ctrl_Tmx.fc","HFD_T2_Ctrl_Tmx.fc",
                
                "Ctrl_T1_CD_HFD.p", "Ctrl_T2_CD_HFD.p",  "Tgt_T1_CD_HFD.p",  "Tgt_T2_CD_HFD.p","Tmx_T1_CD_HFD.p","Tmx_T2_CD_HFD.p",
                "CD_T1_Ctrl_Tgt.p", "CD_T2_Ctrl_Tgt.p","HFD_T1_Ctrl_Tgt.p","HFD_T2_Ctrl_Tgt.p","CD_T1_Ctrl_Tmx.p",
                "CD_T2_Ctrl_Tmx.p", "HFD_T1_Ctrl_Tmx.p","HFD_T2_Ctrl_Tmx.p",
                
                "Ctrl_T1_CD_HFD.FDR","Ctrl_T2_CD_HFD.FDR", "Tgt_T1_CD_HFD.FDR", "Tgt_T2_CD_HFD.FDR",
                "Tmx_T1_CD_HFD.FDR",  "Tmx_T2_CD_HFD.FDR", "CD_T1_Ctrl_Tgt.FDR",       
                "CD_T2_Ctrl_Tgt.FDR", "HFD_T1_Ctrl_Tgt.FDR", "HFD_T2_Ctrl_Tgt.FDR","CD_T1_Ctrl_Tmx.FDR",             
                "CD_T2_Ctrl_Tmx.FDR",  "HFD_T1_Ctrl_Tmx.FDR","HFD_T2_Ctrl_Tmx.FDR",
                
                "Ctrl_T1_CD_HFD","Ctrl_T2_CD_HFD","Tgt_T1_CD_HFD", "Tgt_T2_CD_HFD",      
                "Tmx_T1_CD_HFD", "Tmx_T2_CD_HFD", "CD_T1_Ctrl_Tgt", "CD_T2_Ctrl_Tgt", "HFD_T1_Ctrl_Tgt", "HFD_T2_Ctrl_Tgt",   
                "CD_T1_Ctrl_Tmx", "CD_T2_Ctrl_Tmx", "HFD_T1_Ctrl_Tmx", "HFD_T2_Ctrl_Tmx") 

twodrug_fcpfdr <- twodrug[, names_order]
names(twodrug_fcpfdr)

labsample <- c("Ctrl_T1:HFD vs CD","Ctrl_T2:HFD vs CD","Tgt_T1:HFD vs CD","Tgt_T2:HFD vs CD","Tmx_T1:HFD vs CD",  
                "Tmx_T2:HFD vs CD", "CD_T1:Tgt vs Ctrl", "CD_T2:Tgt vs Ctrl","HFD_T1:Tgt vs Ctrl", "HFD_T2:Tgt vs Ctrl", "CD_T1:Tmx vs Ctrl", 
                                   "CD_T2:Tmx vs Ctrl",  "HFD_T1:Tmx vs Ctrl", "HFD_T2:Tmx vs Ctrl")

#fc: 26 to 39; sig: 68 to 81; p: 40 to 53; fdr:54 to 67
####diets associate
tem <- twodrug_fcpfdr[twodrug_fcpfdr$CD_HFD=="111111" ,]
dim(tem)
names(tem)

pdf("C:/zhijuncao/sun/diet-drug/03182019/ratio_heatmap_diet_68common_T1T2.pdf", width = 12, height =12)
heatmap.2(log(as.matrix(tem[,26:39][,1:6]), 2),col=bluered, Rowv = TRUE, Colv = TRUE, margins = c(8, 44),
          labCol=labsample[1:6],labRow = tem$BIOCHEMICAL,trace="none",  cellnote=tem[,68:81][,1:6],
          na.color = "grey", srtCol=45,cexCol=1, cexRow = 1,keysize=0.8, lhei = c(1,8),
          density.info="none",breaks = seq(-1.5,1.5,0.2),key.xlab="log2ratio", main = NULL)
dev.off()

#####Tgt associate
tem <- twodrug_fcpfdr[twodrug_fcpfdr$Tgt=="1111",]
dim(tem)
pdf("C:/zhijuncao/sun/diet-drug/03182019/ratio_heatmap_label_sig_Tgt_25common.pdf", width = 12, height =12)
heatmap.2(log(as.matrix(tem[,26:39][,7:10]),2),col=bluered, Rowv = TRUE, Colv = TRUE, margins = c(8+58-25, 48),
          labCol=labsample[7:10],labRow = tem$BIOCHEMICAL,trace="none",  cellnote=tem[,68:81][,7:10],notecex=1,
          na.color = "grey", srtCol=45,cexCol=1, cexRow = 1,keysize=0.8, lhei = c(1,8),
          density.info="none",breaks = seq(-1.5,1.5,0.2),key.xlab="log2(Ratio)", main = NULL)
dev.off()


####Tmx associate
tem <- twodrug_fcpfdr[twodrug_fcpfdr$Tmx=="1111",]
dim(tem)
pdf("C:/zhijuncao/sun/diet-drug/03182019/ratio_heatmap_label_sig_Tmx_6common.pdf", width = 12, height =12)
heatmap.2(log(as.matrix(tem[,26:39][,11:14]),2),col=bluered, Rowv = TRUE, Colv = TRUE, margins = c(8+58-8, 48),
          labCol=labsample[11:14],labRow = tem$BIOCHEMICAL,trace="none",  cellnote=tem[,68:81][,11:14],notecex=1,
          na.color = "grey", srtCol=45,cexCol=1, cexRow = 1,keysize=0.8, lhei = c(1,8),
          density.info="none",breaks = seq(-1.5,1.5,0.2),key.xlab="log2(Ratio)", main = NULL)
dev.off()

###Diet and drug related heatmap in pathways
twodrug_fc_sig <- twodrug_fcpfdr %>% filter(sig!="00000000000000")
dim(twodrug_fc)
names(twodrug_fc_sig)
sub_pathway <- unique(twodrug_fc_sig$SUB_PATHWAY)
tem <- twodrug_fc_sig 
fc_sub_diet_heatmap <- function(){
    for (PATHWAY in sig_subpathway){
  tem <- twodrug_fc_sig %>% filter(SUB_PATHWAY==PATHWAY)

  if (nrow(tem)>1){
    heatmap.2(log(as.matrix(tem[,26:39][,1:6]),2),col=bluered, Rowv = TRUE, Colv = F, margins = c(12+48-nrow(tem), 44),
              labCol=labsample[1:6],labRow = tem$BIOCHEMICAL,trace="none",  cellnote=tem[,68:81][,1:6],
              na.color = "grey", srtCol=45,cexCol=1,cexRow = 1,keysize=0.8,lhei = c(1,8),symkey=F,
              density.info="none",breaks = seq(-1.5,1.5,0.2),key.xlab="log2(Ratio)", main = PATHWAY)
  }
}
}

pdf("C:/zhijuncao/sun/diet-drug/03182019/ratio_heatmap_label_sig_diet_33sub_pathway_colVfalse.pdf", width = 12, height =12)
par(cex.main=1)
fc_sub_diet_heatmap()
dev.off()

fc_sub_drug_heatmap <- function(){
  for (PATHWAY in sig_subpathway){
    tem <- twodrug_fc_sig %>% filter(SUB_PATHWAY==PATHWAY)
    
    if (nrow(tem)>1){
      heatmap.2(log(as.matrix(tem[,26:39][,7:14]),2),col=bluered, Rowv = TRUE, Colv = F, margins = c(12+48-nrow(tem), 40),
                labCol=labsample[7:14],labRow = tem$BIOCHEMICAL,trace="none",  cellnote=tem[,68:81][,7:14],
                na.color = "grey", srtCol=45,cexCol=1,cexRow = 1,keysize=0.8,lhei = c(1,8),symkey=F,
                density.info="none",breaks = seq(-1.5,1.5,0.2),key.xlab="log2(Ratio)", main = PATHWAY)
    }
  }
}
pdf("C:/zhijuncao/sun/diet-drug/03182019/ratio_heatmap_label_sig_drug_33sub_pathway_colvfalse.pdf", width = 12, height =12)
par(cex.main=1)
fc_sub_drug_heatmap()
dev.off()

fc_sub_tgt_heatmap <- function(){
  for (PATHWAY in sig_subpathway){
    tem <- twodrug_fc_sig %>% filter(SUB_PATHWAY==PATHWAY)
    
    if (nrow(tem)>1){
      heatmap.2(log(as.matrix(tem[,26:39][,7:10]),2),col=bluered, Rowv = TRUE, Colv = T, margins = c(12+48-nrow(tem), 48),
                labCol=labsample[7:10],labRow = tem$BIOCHEMICAL,trace="none",  cellnote=tem[,68:81][,7:10],
                na.color = "grey", srtCol=45,cexCol=1,cexRow = 1,keysize=0.8,lhei = c(1,8),symkey=F,
                density.info="none",breaks = seq(-1.5,1.5,0.2),key.xlab="log2(Ratio)", main = PATHWAY)
    }
  }
}
pdf("C:/zhijuncao/sun/diet-drug/03182019/ratio_heatmap_label_sig_tgt_33sub_pathway.pdf", width = 12, height =12)
par(cex.main=1)
fc_sub_tgt_heatmap()
dev.off()

fc_sub_tmx_heatmap <- function(){
  for (PATHWAY in sig_subpathway){
    tem <- twodrug_fc_sig %>% filter(SUB_PATHWAY==PATHWAY)
    
    if (nrow(tem)>1){
      heatmap.2(log(as.matrix(tem[,26:39][,11:14]),2),col=bluered, Rowv = TRUE, Colv = T, margins = c(12+48-nrow(tem), 48),
                labCol=labsample[11:14],labRow = tem$BIOCHEMICAL,trace="none",  cellnote=tem[,68:81][,11:14],
                na.color = "grey", srtCol=45,cexCol=1,cexRow = 1,keysize=0.8,lhei = c(1,8),symkey=F,
                density.info="none",breaks = seq(-1.5,1.5,0.2),key.xlab="log2(Ratio)", main = PATHWAY)
    }
  }
}
pdf("C:/zhijuncao/sun/diet-drug/03182019/ratio_heatmap_label_sig_tmx_33sub_pathway.pdf", width = 12, height =12)
par(cex.main=1)
fc_sub_tmx_heatmap()
dev.off()




names(twodrug_fc_sig)

pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_sub_pathway_ratioa.pdf", width = 12, height =12)
fc_sub_heatmap(name=sample_names)
dev.off()
pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_sub_pathway_t1_ratioa.pdf", width = 12, height =12)
fc_sub_heatmap(name=t1)
dev.off()
pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_sub_pathway_t2_ratioa.pdf", width = 12, height =12)
fc_sub_heatmap(name=t2)
dev.off()

pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_sub_pathway_t1_tgt_ratioa.pdf", width = 12, height =12)
fc_sub_heatmap(name=t1_tgt)
dev.off()
pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_sub_pathway_t2_tgt_ratioa.pdf", width = 12, height =12)
fc_sub_heatmap(name=t2_tgt)
dev.off()

pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_sub_pathway_t1_tmx_ratioa.pdf", width = 12, height =12)
fc_sub_heatmap(name=t1_tmx)
dev.off()
pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_sub_pathway_t2_tmx_ratioa.pdf", width = 12, height =12)
fc_sub_heatmap(name=t2_tmx)
dev.off()


#margins = c(12+38-nrow(tem), 22+14-min(14,length(name))
fc_sup_heatmap <- function(name=sample_names){
  for (PATHWAY in super_pathway){
    tem <- twodrug_fc_sig %>% filter(SUPER_PATHWAY==PATHWAY)
 
    if (nrow(tem)>1){
      heatmap.2(log(as.matrix(tem[,name]),2),col=bluered, Rowv = TRUE, Colv = TRUE, margins = c(12+48-nrow(tem), 22+14-min(14, length(name))),
                labCol=name,labRow = tem$BIOCHEMICAL,trace="none",  cellnote=round(as.matrix(tem[,name]),2),
                na.color = "grey", srtCol=45,cexRow = 0.8,keysize=0.8,lhei=c(1,20),notecex=0.8,
                density.info="none",breaks = seq(-2,2,0.2),key.xlab="log2(Ratio)", main = PATHWAY)
    }
  }
}

pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_sup_pathway_ratioa.pdf", width = 12, height =12)
fc_sup_heatmap(name=sample_names)
dev.off()
pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_sup_pathway_t1_ratioa.pdf", width = 12, height =12)
fc_sup_heatmap(name=t1)
dev.off()
pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_sup_pathway_t2_ratioa.pdf", width = 12, height =12)
fc_sup_heatmap(name=t2)
dev.off()

pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_sup_pathway_t1_tgt_ratioa.pdf", width = 12, height =12)
fc_sup_heatmap(name=t1_tgt)
dev.off()
pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_sup_pathway_t2_tgt_ratioa.pdf", width = 12, height =12)
fc_sup_heatmap(name=t2_tgt)
dev.off()

pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_sup_pathway_t1_tmx_ratioa.pdf", width = 12, height =12)
fc_sup_heatmap(name=t1_tmx)
dev.off()
pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_sup_pathway_t2_tmx_ratioa.pdf", width = 12, height =12)
fc_sup_heatmap(name=t2_tmx)
dev.off()

#lipid and aa
pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_lipidaa_sup_pathway_ratio.pdf", width = 12, height =24)
fc_sup_heatmap(name=sample_names)
dev.off()
pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_lipidaa_sup_pathway_t1_ratio.pdf", width = 12, height =24)
fc_sup_heatmap(name=t1)
dev.off()
pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_lipidaa_sup_pathway_t2_ratio.pdf", width = 12, height =24)
fc_sup_heatmap(name=t2)
dev.off()

pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_lipidaa_sup_pathway_t1_tgt_ratio.pdf", width = 12, height =24)
fc_sup_heatmap(name=t1_tgt)
dev.off()
pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_lipidaa_sup_pathway_t2_tgt_ratio.pdf", width = 12, height =24)
fc_sup_heatmap(name=t2_tgt)
dev.off()

pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_lipidaa_sup_pathway_t1_tmx_ratio.pdf", width = 12, height =24)
fc_sup_heatmap(name=t1_tmx)
dev.off()
pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_lipidaa_sup_pathway_t2_tmx_ratio.pdf", width = 12, height =24)
fc_sup_heatmap(name=t2_tmx)
dev.off()



PATHWAY <- "Tocopherol Metabolism"
pdf("C:/zhijuncao/sun/diet-drug/112018/heatmap_fold_change_log2_sub_pathway_ratio.pdf", width = 12, height =12)

for (PATHWAY in sub_pathway){
  tem <- twodrug_fc_sig %>% filter(SUB_PATHWAY==PATHWAY)
  if (nrow(tem)>1){
      heatmap.2(log(as.matrix(tem[,21:34]),2),col=bluered, Rowv = TRUE, Colv = TRUE, margins = c(12+48-nrow(tem), 22),
            labCol=names(tem[,21:34]),labRow = tem$BIOCHEMICAL,trace="none",  cellnote=round(as.matrix(tem[,21:34]),2),
            na.color = "grey", srtCol=45,cexRow = 1.2,keysize=0.8,
            density.info="none",breaks = seq(-2,2,0.2),key.xlab="log2(Ratio)", main = PATHWAY)
  }

}
dev.off()


pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_super_pathway_ratio.pdf", width = 12, height =12)
pdf("C:/zhijuncao/sun/diet-drug/01282019/heatmap_fold_change_log2_super_lipid_aa_ratio_12x25.pdf", width = 12, height =25)
for (PATHWAY in super_pathway){
  tem <- twodrug_fc_sig %>% filter(SUPER_PATHWAY==PATHWAY)
  if (nrow(tem)>1){
    heatmap.2(log(as.matrix(tem[,21:34]),2),col=bluered, Rowv = TRUE, Colv = TRUE, margins = c(12+ 38-min(38, nrow(tem)), 22),
              labCol=names(tem[,21:34]),labRow = tem$BIOCHEMICAL,trace="none", cellnote=round(as.matrix(tem[,21:34]),2), 
              na.color = "grey", srtCol=45,cexRow = 0.8,keysize=0.8,lhei=c(1,20),notecex=0.8,
              density.info="none",breaks = seq(-2,2,0.2), key.xlab="log2(Ratio)", main = PATHWAY)
  }
  
}
dev.off()






names(twodrug_fc_sig[,21:34])



pdf("heatmap_fold_change_log2.pdf", width = 12, height =12)
heatmap.2(log(as.matrix(tem),2),col=bluered, Rowv = TRUE, Colv = TRUE, margins = c(12, 22),
          labCol=names(tem),labRow = FALSE,trace="none", RowSideColors=sc[metaid$SUPER_PATHWAY], 
          ColSideColors=cc[names(tem)],na.color = "grey", srtCol=45,cexRow = 0.8,keysize=0.8,
          density.info="none",breaks = seq(-2,2,0.2))
legend("topright", legend = names(cc)[1:12],fill = cc,title = "sample", cex = 1)
legend("bottomright", legend = names(sc)[1:8],fill = sc,title = "SUPER_PATHWAY", cex = 1)
dev.off()


name(re1)

dim(sample)
names(sample)[1:16]

#group data extraction

names(sample)[1:16]
levels(sample$SampleID)
Ctrl_CD_T1 <- sample %>% filter(SampleID=='Ctrl_CD_T1')
Ctrl_CD_T2 <- sample %>% filter(SampleID=='Ctrl_CD_T2')
Ctrl_HFD_T1 <- sample %>% filter(SampleID=='Ctrl_HFD_T1')
Ctrl_HFD_T2 <- sample %>% filter(SampleID=='Ctrl_HFD_T2')

Tgt_CD_T1 <- sample %>% filter(SampleID=='Tgt_CD_T1')
Tgt_CD_T2 <- sample %>% filter(SampleID=='Tgt_CD_T2')
Tgt_HFD_T1 <- sample %>% filter(SampleID=='Tgt_HFD_T1')
Tgt_HFD_T2 <- sample %>% filter(SampleID=='Tgt_HFD_T2')

Tmx_CD_T1 <- sample %>% filter(SampleID=='Tmx_CD_T1')
Tmx_CD_T2 <- sample %>% filter(SampleID=='Tmx_CD_T2')
Tmx_HFD_T1 <- sample %>% filter(SampleID=='Tmx_HFD_T1')
Tmx_HFD_T2 <- sample %>% filter(SampleID=='Tmx_HFD_T2')

####Venn plot########
twodrug_re <- read.xlsx("C:/zhijuncao/sun/diet-drug/01282019/two_Drug_all_p_fc_fdr_sig_01label_020619.xlsx", sheet = "meta565")

#######get charater by position function
getChar <- function(x, positions){
  # x <- "abcdef"
  # positions <- c(1,3,5)
  re <- vector(length = length(positions))
  i <- 1
  for (position in positions){
    re[i] <- substr(x, position, position)
    i <- i+1
  }
  re <- paste(re, collapse = "")
}
vgetChar <- Vectorize(getChar, vectorize.args = "x", USE.NAMES = FALSE)
###################

twodrug_re[, "CD_HFD_T1"] <- vgetChar(twodrug_re$CD_HFD, c(1,3,5))
twodrug_re[, "Ctrl_T1_CD_HFD"] <- vgetChar(twodrug_re$CD_HFD, 1)
twodrug_re[, "Tgt_T1_CD_HFD"] <- vgetChar(twodrug_re$CD_HFD, 3)
twodrug_re[, "Tmx_T1_CD_HFD"] <- vgetChar(twodrug_re$CD_HFD, 5)

twodrug_re[, "CD_HFD_T2"] <- vgetChar(twodrug_re$CD_HFD,c(2,4,6))
twodrug_re[, "Ctrl_T2_CD_HFD"] <- vgetChar(twodrug_re$CD_HFD, 2)
twodrug_re[, "Tgt_T2_CD_HFD"] <- vgetChar(twodrug_re$CD_HFD, 4)
twodrug_re[, "Tmx_T2_CD_HFD"] <- vgetChar(twodrug_re$CD_HFD, 6)


twodrug_re[, "Tgt_T1"] <- vgetChar(twodrug_re$Tgt, c(1,3))
twodrug_re[, "Tgt_T2"] <- vgetChar(twodrug_re$Tgt, c(2,4))
twodrug_re[, "Tmx_T1"] <- vgetChar(twodrug_re$Tmx, c(1,3))
twodrug_re[, "Tmx_T2"] <- vgetChar(twodrug_re$Tmx, c(2,4))

twodrug_re[, "CD_T1_Tgt"] <- vgetChar(twodrug_re$Tgt, 1)
twodrug_re[, "HFD_T1_Tgt"] <- vgetChar(twodrug_re$Tgt, 3)

twodrug_re[, "CD_T2_Tgt"] <- vgetChar(twodrug_re$Tgt, 2)
twodrug_re[, "HFD_T2_Tgt"] <- vgetChar(twodrug_re$Tgt, 4)

twodrug_re[, "CD_T1_Tmx"] <- vgetChar(twodrug_re$Tmx, 1)
twodrug_re[, "HFD_T1_Tmx"] <- vgetChar(twodrug_re$Tmx, 3)

twodrug_re[, "CD_T2_Tmx"] <- vgetChar(twodrug_re$Tmx, 2)
twodrug_re[, "HFD_T2_Tmx"] <- vgetChar(twodrug_re$Tmx, 4)

write.xlsx(twodrug_re, "c:/zhijuncao/sun/diet-drug/01282019/two_drug_all_p_fc_fdr_01lable_02132019a.xlsx")
twodrug_re <- read.xlsx("c:/zhijuncao/sun/diet-drug/01282019/two_drug_all_p_fc_fdr_01lable_02132019a.xlsx")
names(twodrug_re)

CD_HFD_T1_sig <- twodrug_re$Name[twodrug_re$CD_HFD_T1!="000"]
Ctrl_T1_CD_HFD_sig <- twodrug_re$Name[twodrug_re$Ctrl_T1_CD_HFD =="1"]
Tgt_T1_CD_HFD_sig <- twodrug_re$Name[twodrug_re$Tgt_T1_CD_HFD =="1"]
Tmx_T1_CD_HFD_sig <- twodrug_re$Name[twodrug_re$Tmx_T1_CD_HFD =="1"]

CD_HFD_T2_sig <- twodrug_re$Name[twodrug_re$CD_HFD_T2!="000"]
Ctrl_T2_CD_HFD_sig <- twodrug_re$Name[twodrug_re$Ctrl_T2_CD_HFD =="1"]
Tgt_T2_CD_HFD_sig <- twodrug_re$Name[twodrug_re$Tgt_T2_CD_HFD =="1"]
Tmx_T2_CD_HFD_sig <- twodrug_re$Name[twodrug_re$Tmx_T2_CD_HFD =="1"]

vennplot(list(Ctrl=Ctrl_T1_CD_HFD_sig, Tgt=Tgt_T1_CD_HFD_sig, Tmx=Tmx_T1_CD_HFD_sig),savename = "CD_HFD_T1", labelsize = 1.5 )
vennplot(list(Ctrl=Ctrl_T2_CD_HFD_sig, Tgt=Tgt_T2_CD_HFD_sig, Tmx=Tmx_T2_CD_HFD_sig),savename = "CD_HFD_T2", labelsize =1.5 )


Tgt_T1_sig <- twodrug_re$Name[twodrug_re$Tgt_T1!="00"]
Tgt_T2_sig <- twodrug_re$Name[twodrug_re$Tgt_T2!="00"]

Tmx_T1_sig <- twodrug_re$Name[twodrug_re$Tmx_T1!="00"]
Tmx_T2_sig <- twodrug_re$Name[twodrug_re$Tmx_T2!="00"]

CD_T1_Tgt_sig <- twodrug_re$Name[twodrug_re$CD_T1_Tgt=="1"]
HFD_T1_Tgt_sig <- twodrug_re$Name[twodrug_re$HFD_T1_Tgt=="1"]
CD_T2_Tgt_sig <- twodrug_re$Name[twodrug_re$CD_T2_Tgt=="1"]
HFD_T2_Tgt_sig <- twodrug_re$Name[twodrug_re$HFD_T2_Tgt=="1"]

CD_T1_Tmx_sig <- twodrug_re$Name[twodrug_re$CD_T1_Tmx=="1"]
HFD_T1_Tmx_sig <- twodrug_re$Name[twodrug_re$HFD_T1_Tmx=="1"]
CD_T2_Tmx_sig <- twodrug_re$Name[twodrug_re$CD_T2_Tmx=="1"]
HFD_T2_Tmx_sig <- twodrug_re$Name[twodrug_re$HFD_T2_Tmx=="1"]

vennplot(list(CD_HFD_T1=CD_HFD_T1_sig, CD_HFD_T2=CD_HFD_T2_sig),savename = "CD_HFD_T1_T2" )

vennplot(list(Tgt_T1=Tgt_T1_sig, Tgt_T2=Tgt_T2_sig),savename = "Tgt_T1_T2" )
vennplot(list(Tmx_T1=Tmx_T1_sig, Tmx_T2=Tmx_T2_sig),savename = "Tmx_T1_T2" )

vennplot(list(CD_HFD_T1=CD_HFD_T1_sig, Tgt_T1=Tgt_T1_sig, Tmx_T1=Tmx_T1_sig),savename = "T1" )
vennplot(list(CD_HFD_T2=CD_HFD_T2_sig, Tgt_T2=Tgt_T2_sig, Tmx_T2=Tmx_T2_sig),savename = "T2" )

vennplot(list(CD_T1=CD_T1_Tgt_sig, HFD_T1=HFD_T1_Tgt_sig, CD_T2=CD_T2_Tgt_sig, HFD_T2=HFD_T2_Tgt_sig),savename = "Tgt", labelsize = 0.8 )
vennplot(list(CD_T1=CD_T1_Tmx_sig, HFD_T1=HFD_T1_Tmx_sig, CD_T2=CD_T2_Tmx_sig, HFD_T2=HFD_T2_Tmx_sig),savename = "Tmx", labelsize = 0.8 )

vennplot(list(CD_T1=CD_T1_Tgt_sig, HFD_T1=HFD_T1_Tgt_sig),savename = "Tgt_T1" )
vennplot(list(CD_T2=CD_T2_Tgt_sig, HFD_T2=HFD_T2_Tgt_sig),savename = "Tgt_T2" )

vennplot(list(CD_T1=CD_T1_Tmx_sig, HFD_T1=HFD_T1_Tmx_sig),savename = "Tmx_T1" )
vennplot(list( CD_T2=CD_T2_Tmx_sig, HFD_T2=HFD_T2_Tmx_sig),savename = "Tmx_T2" )

vennplot(list(CD_HFD=CD_HFD_T1_sig, 
              CD_Tgt=CD_T1_Tgt_sig, HFD_Tgt=HFD_T1_Tgt_sig,
              CD_Tmx=CD_T1_Tmx_sig, HFD_Tmx=HFD_T1_Tmx_sig),savename = "T1_CD_HFD_Tgt_Tmx", labelsize = 1 )

vennplot(list(CD_HFD=CD_HFD_T2_sig, 
              CD_Tgt=CD_T2_Tgt_sig, HFD_Tgt=HFD_T2_Tgt_sig,
              CD_Tmx=CD_T2_Tmx_sig, HFD_Tmx=HFD_T2_Tmx_sig),savename = "T2_CD_HFD_Tgt_Tmx" , labelsize = 1)



##############################################
#t test
##############################################
#name <- "pinitol"
my.t.test.p.value <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

i <- 39
length(names(sample)[1:16])

re <- list()

for (i in 16:583){
#CD vs HFD
Ctrl_T1_CD_HFD.p <- my.t.test.p.value(log(Ctrl_CD_T1[,i],2),log(Ctrl_HFD_T1[,i],2),var.equal = TRUE, alternative = "two.sided")
Ctrl_T1_CD_HFD.fc <- mean(Ctrl_HFD_T1[,i])/mean(Ctrl_CD_T1[,i])
Ctrl_T2_CD_HFD.p <- my.t.test.p.value(log(Ctrl_CD_T2[,i],2),log(Ctrl_HFD_T2[,i],2), var.equal = TRUE,alternative = "two.sided")
Ctrl_T2_CD_HFD.fc <- mean(Ctrl_HFD_T2[,i])/mean(Ctrl_CD_T2[,i])

Tgt_T1_CD_HFD.p <- my.t.test.p.value(log(Tgt_CD_T1[,i],2),log(Tgt_HFD_T1[,i],2), var.equal = TRUE,alternative = "two.sided")
Tgt_T1_CD_HFD.fc <- mean(Tgt_HFD_T1[,i])/mean(Tgt_CD_T1[,i])
Tgt_T2_CD_HFD.p <- my.t.test.p.value(log(Tgt_CD_T2[,i],2),log(Tgt_HFD_T2[,i],2), var.equal = TRUE,alternative = "two.sided")
Tgt_T2_CD_HFD.fc <- mean(Tgt_HFD_T2[,i])/mean(Tgt_CD_T2[,i])

Tmx_T1_CD_HFD.p <- my.t.test.p.value(log(Tmx_CD_T1[,i],2),log(Tmx_HFD_T1[,i],2), var.equal = TRUE,alternative = "two.sided")
Tmx_T1_CD_HFD.fc <- mean(Tmx_HFD_T1[,i])/mean(Tmx_CD_T1[,i])
Tmx_T2_CD_HFD.p <- my.t.test.p.value(log(Tmx_CD_T2[,i],2),log(Tmx_HFD_T2[,i],2), var.equal = TRUE,alternative = "two.sided")
Tmx_T2_CD_HFD.fc <- mean(Tmx_HFD_T2[,i])/mean(Tmx_CD_T2[,i])

#Control vs treatment
CD_T1_Ctrl_Tgt.p <- my.t.test.p.value(log(Ctrl_CD_T1[,i],2),log(Tgt_CD_T1[,i],2), var.equal = TRUE,alternative = "two.sided")
CD_T1_Ctrl_Tgt.fc <- mean(Tgt_CD_T1[,i])/mean(Ctrl_CD_T1[,i])
CD_T2_Ctrl_Tgt.p <- my.t.test.p.value(log(Ctrl_CD_T2[,i],2),log(Tgt_CD_T2[,i],2), var.equal = TRUE,alternative = "two.sided")
CD_T2_Ctrl_Tgt.fc <- mean(Tgt_CD_T2[,i])/mean(Ctrl_CD_T2[,i])

HFD_T1_Ctrl_Tgt.p <- my.t.test.p.value(log(Ctrl_HFD_T1[,i],2),log(Tgt_HFD_T1[,i],2), var.equal = TRUE,alternative = "two.sided")
HFD_T1_Ctrl_Tgt.fc <- mean(Tgt_HFD_T1[,i])/mean(Ctrl_HFD_T1[,i])
HFD_T2_Ctrl_Tgt.p <- my.t.test.p.value(log(Ctrl_HFD_T2[,i],2),log(Tgt_HFD_T2[,i],2), var.equal = TRUE,alternative = "two.sided")
HFD_T2_Ctrl_Tgt.fc <- mean(Tgt_HFD_T2[,i])/mean(Ctrl_HFD_T2[,i])

CD_T1_Ctrl_Tmx.p <- my.t.test.p.value(log(Ctrl_CD_T1[,i],2),log(Tmx_CD_T1[,i],2), var.equal = TRUE,alternative = "two.sided")
CD_T1_Ctrl_Tmx.fc <- mean(Tmx_CD_T1[,i])/mean(Ctrl_CD_T1[,i])
CD_T2_Ctrl_Tmx.p <- my.t.test.p.value(log(Ctrl_CD_T2[,i],2),log(Tmx_CD_T2[,i],2), var.equal = TRUE,alternative = "two.sided")
CD_T2_Ctrl_Tmx.fc <- mean(Tmx_CD_T2[,i])/mean(Ctrl_CD_T2[,i])

HFD_T1_Ctrl_Tmx.p <- my.t.test.p.value(log(Ctrl_HFD_T1[,i],2),log(Tmx_HFD_T1[,i],2), var.equal = TRUE,alternative = "two.sided")
HFD_T1_Ctrl_Tmx.fc <- mean(Tmx_HFD_T1[,i])/mean(Ctrl_HFD_T1[,i])
HFD_T2_Ctrl_Tmx.p <- my.t.test.p.value(log(Ctrl_HFD_T2[,i],2),log(Tmx_HFD_T2[,i],2), var.equal = TRUE,alternative = "two.sided")
HFD_T2_Ctrl_Tmx.fc <- mean(Tmx_HFD_T2[,i])/mean(Ctrl_HFD_T2[,i])


re[[i-15]] <- c(metabolite=names(Ctrl_CD_T1)[i],
                Ctrl_T1_CD_HFD.p=Ctrl_T1_CD_HFD.p, Ctrl_T1_CD_HFD.fc=Ctrl_T1_CD_HFD.fc,
                Ctrl_T2_CD_HFD.p=Ctrl_T2_CD_HFD.p, Ctrl_T2_CD_HFD.fc=Ctrl_T2_CD_HFD.fc,
                
                Tgt_T1_CD_HFD.p=Tgt_T1_CD_HFD.p, Tgt_T1_CD_HFD.fc=Tgt_T1_CD_HFD.fc,
                Tgt_T2_CD_HFD.p=Tgt_T2_CD_HFD.p, Tgt_T2_CD_HFD.fc=Tgt_T2_CD_HFD.fc,
                Tmx_T1_CD_HFD.p=Tmx_T1_CD_HFD.p, Tmx_T1_CD_HFD.fc=Tmx_T1_CD_HFD.fc,
                Tmx_T2_CD_HFD.p=Tmx_T2_CD_HFD.p, Tmx_T2_CD_HFD.fc=Tmx_T2_CD_HFD.fc,
                
                CD_T1_Ctrl_Tgt.p=CD_T1_Ctrl_Tgt.p, CD_T1_Ctrl_Tgt.fc=CD_T1_Ctrl_Tgt.fc,
                CD_T2_Ctrl_Tgt.p=CD_T2_Ctrl_Tgt.p, CD_T2_Ctrl_Tgt.fc=CD_T2_Ctrl_Tgt.fc,
                HFD_T1_Ctrl_Tgt.p=HFD_T1_Ctrl_Tgt.p, HFD_T1_Ctrl_Tgt.fc=HFD_T1_Ctrl_Tgt.fc,
                HFD_T2_Ctrl_Tgt.p=HFD_T2_Ctrl_Tgt.p, HFD_T2_Ctrl_Tgt.fc=HFD_T2_Ctrl_Tgt.fc,
                CD_T1_Ctrl_Tmx.p=CD_T1_Ctrl_Tmx.p, CD_T1_Ctrl_Tmx.fc=CD_T1_Ctrl_Tmx.fc,
                CD_T2_Ctrl_Tmx.p=CD_T2_Ctrl_Tmx.p, CD_T2_Ctrl_Tmx.fc=CD_T2_Ctrl_Tmx.fc,
                HFD_T1_Ctrl_Tmx.p=HFD_T1_Ctrl_Tmx.p, HFD_T1_Ctrl_Tmx.fc=HFD_T1_Ctrl_Tmx.fc,
                HFD_T2_Ctrl_Tmx.p=HFD_T2_Ctrl_Tmx.p, HFD_T2_Ctrl_Tmx.fc=HFD_T2_Ctrl_Tmx.fc)
cat(i)
}
re[[1]]
result1 <- data.frame(do.call(rbind, re), stringsAsFactors = FALSE)
dim(result1)
write.csv(result1, "two_Drug_statistic_p_fc_568metabolite_020619.csv")

result1 <- read.csv("two_Drug_statistic_p_fc_568metabolite_020619.csv", stringsAsFactors = FALSE)
str(result1)
str
str(re)
result1[,"Ctrl_T1_CD_HFD.FDR"] <- p.adjust(result1$Ctrl_T1_CD_HFD.p, method = "BH")
result1[,"Ctrl_T2_CD_HFD.FDR"] <- p.adjust(result1$Ctrl_T2_CD_HFD.p, method = "BH")

result1[,"Tgt_T1_CD_HFD.FDR"] <- p.adjust(result1$Tgt_T1_CD_HFD.p, method = "BH")
result1[,"Tgt_T2_CD_HFD.FDR"] <- p.adjust(result1$Tgt_T2_CD_HFD.p, method = "BH")
result1[,"Tmx_T1_CD_HFD.FDR"] <- p.adjust(result1$Tmx_T1_CD_HFD.p, method = "BH")
result1[,"Tmx_T2_CD_HFD.FDR"] <- p.adjust(result1$Tmx_T2_CD_HFD.p, method = "BH")


result1[,"CD_T1_Ctrl_Tgt.FDR"] <- p.adjust(result1$CD_T1_Ctrl_Tgt.p, method = "BH")
result1[,"CD_T2_Ctrl_Tgt.FDR"] <- p.adjust(result1$CD_T2_Ctrl_Tgt.p, method = "BH")
result1[,"HFD_T1_Ctrl_Tgt.FDR"] <- p.adjust(result1$HFD_T1_Ctrl_Tgt.p, method = "BH")
result1[,"HFD_T2_Ctrl_Tgt.FDR"] <- p.adjust(result1$HFD_T2_Ctrl_Tgt.p, method = "BH")

result1[,"CD_T1_Ctrl_Tmx.FDR"] <- p.adjust(result1$CD_T1_Ctrl_Tmx.p, method = "BH")
result1[,"CD_T2_Ctrl_Tmx.FDR"] <- p.adjust(result1$CD_T2_Ctrl_Tmx.p, method = "BH")
result1[,"HFD_T1_Ctrl_Tmx.FDR"] <- p.adjust(result1$HFD_T1_Ctrl_Tmx.p, method = "BH")
result1[,"HFD_T2_Ctrl_Tmx.FDR"] <- p.adjust(result1$HFD_T2_Ctrl_Tmx.p, method = "BH")


result2 <- result1[,c("metabolite",
                      "Ctrl_T1_CD_HFD.p", "Ctrl_T1_CD_HFD.fc", "Ctrl_T1_CD_HFD.FDR", "Ctrl_T2_CD_HFD.p", "Ctrl_T2_CD_HFD.fc","Ctrl_T2_CD_HFD.FDR",  
                      
                      "Tgt_T1_CD_HFD.p", "Tgt_T1_CD_HFD.fc", "Tgt_T1_CD_HFD.FDR","Tgt_T2_CD_HFD.p", "Tgt_T2_CD_HFD.fc", "Tgt_T2_CD_HFD.FDR",
                      "Tmx_T1_CD_HFD.p", "Tmx_T1_CD_HFD.fc", "Tmx_T1_CD_HFD.FDR", "Tmx_T2_CD_HFD.p", "Tmx_T2_CD_HFD.fc", "Tmx_T2_CD_HFD.FDR",
                      
                      "CD_T1_Ctrl_Tgt.p", "CD_T1_Ctrl_Tgt.fc", "CD_T1_Ctrl_Tgt.FDR", "CD_T2_Ctrl_Tgt.p", "CD_T2_Ctrl_Tgt.fc", "CD_T2_Ctrl_Tgt.FDR", 
                      "HFD_T1_Ctrl_Tgt.p", "HFD_T1_Ctrl_Tgt.fc","HFD_T1_Ctrl_Tgt.FDR", "HFD_T2_Ctrl_Tgt.p", "HFD_T2_Ctrl_Tgt.fc", "HFD_T2_Ctrl_Tgt.FDR",
                      "CD_T1_Ctrl_Tmx.p",  "CD_T1_Ctrl_Tmx.fc", "CD_T1_Ctrl_Tmx.FDR", "CD_T2_Ctrl_Tmx.p", "CD_T2_Ctrl_Tmx.fc", "CD_T2_Ctrl_Tmx.FDR",
                      "HFD_T1_Ctrl_Tmx.p", "HFD_T1_Ctrl_Tmx.fc", "HFD_T1_Ctrl_Tmx.FDR", "HFD_T2_Ctrl_Tmx.p", "HFD_T2_Ctrl_Tmx.fc","HFD_T2_Ctrl_Tmx.FDR")]


result3 <- inner_join(metaid,result2, by=c("Name"="metabolite"))
write.xlsx(result3,"two_Drug_p_fc_fdr_568metabolites.xlsx")
tem <- result3
result3<- read.xlsx("two_Drug_all_p_fc_fdr_sig_01label.xlsx")
result3<- read.xlsx("statistic_p_fc_fdr.xlsx")
names(result3)
pvalue <- 0.05
fc <- 1.2
fdr <- 0.2

str(result3)

Ctrl_T1_CD_HFD.sig <- result3 %>% filter((Ctrl_T1_CD_HFD.p<pvalue&(Ctrl_T1_CD_HFD.fc>=fc|Ctrl_T1_CD_HFD.fc<=1/fc)&Ctrl_T1_CD_HFD.FDR<fdr))
Ctrl_T2_CD_HFD.sig <- result3 %>% filter((Ctrl_T2_CD_HFD.p<pvalue&(Ctrl_T2_CD_HFD.fc>=fc|Ctrl_T2_CD_HFD.fc<=1/fc)&Ctrl_T2_CD_HFD.FDR<fdr))

Tgt_T1_CD_HFD.sig <- result3 %>% filter((Tgt_T1_CD_HFD.p<pvalue&(Tgt_T1_CD_HFD.fc>=fc|Tgt_T1_CD_HFD.fc<=1/fc)&Tgt_T1_CD_HFD.FDR<fdr))
Tgt_T2_CD_HFD.sig <- result3 %>% filter((Tgt_T2_CD_HFD.p<pvalue&(Tgt_T2_CD_HFD.fc>=fc|Tgt_T2_CD_HFD.fc<=1/fc)&Tgt_T2_CD_HFD.FDR<fdr))
Tmx_T1_CD_HFD.sig <- result3 %>% filter((Tmx_T1_CD_HFD.p<pvalue&(Tmx_T1_CD_HFD.fc>=fc|Tmx_T1_CD_HFD.fc<=1/fc)&Tmx_T1_CD_HFD.FDR<fdr))
Tmx_T2_CD_HFD.sig <- result3 %>% filter((Tmx_T2_CD_HFD.p<pvalue&(Tmx_T2_CD_HFD.fc>=fc|Tmx_T2_CD_HFD.fc<=1/fc)&Tmx_T2_CD_HFD.FDR<fdr))

CD_T1_Ctrl_Tgt.sig <- result3 %>% filter((CD_T1_Ctrl_Tgt.p<pvalue&(CD_T1_Ctrl_Tgt.fc>=fc|CD_T1_Ctrl_Tgt.fc<=1/fc)&CD_T1_Ctrl_Tgt.FDR<fdr))
CD_T2_Ctrl_Tgt.sig <- result3 %>% filter((CD_T2_Ctrl_Tgt.p<pvalue&(CD_T2_Ctrl_Tgt.fc>=fc|CD_T2_Ctrl_Tgt.fc<=1/fc)&CD_T2_Ctrl_Tgt.FDR<fdr))
HFD_T1_Ctrl_Tgt.sig <- result3 %>% filter((HFD_T1_Ctrl_Tgt.p<pvalue&(HFD_T1_Ctrl_Tgt.fc>=fc|HFD_T1_Ctrl_Tgt.fc<=1/fc)&HFD_T1_Ctrl_Tgt.FDR<fdr))
HFD_T2_Ctrl_Tgt.sig <- result3 %>% filter((HFD_T2_Ctrl_Tgt.p<pvalue&(HFD_T2_Ctrl_Tgt.fc>=fc|HFD_T2_Ctrl_Tgt.fc<=1/fc)&HFD_T2_Ctrl_Tgt.FDR<fdr))

CD_T1_Ctrl_Tmx.sig <- result3 %>% filter((CD_T1_Ctrl_Tmx.p<pvalue&(CD_T1_Ctrl_Tmx.fc>=fc|CD_T1_Ctrl_Tmx.fc<=1/fc)&CD_T1_Ctrl_Tmx.FDR<fdr))
CD_T2_Ctrl_Tmx.sig <- result3 %>% filter((CD_T2_Ctrl_Tmx.p<pvalue&(CD_T2_Ctrl_Tmx.fc>=fc|CD_T2_Ctrl_Tmx.fc<=1/fc)&CD_T2_Ctrl_Tmx.FDR<fdr))
HFD_T1_Ctrl_Tmx.sig <- result3 %>% filter((HFD_T1_Ctrl_Tmx.p<pvalue&(HFD_T1_Ctrl_Tmx.fc>=fc|HFD_T1_Ctrl_Tmx.fc<=1/fc)&HFD_T1_Ctrl_Tmx.FDR<fdr))
HFD_T2_Ctrl_Tmx.sig <- result3 %>% filter((HFD_T2_Ctrl_Tmx.p<pvalue&(HFD_T2_Ctrl_Tmx.fc>=fc|HFD_T2_Ctrl_Tmx.fc<=1/fc)&HFD_T2_Ctrl_Tmx.FDR<fdr))

#two drugs
all_sig <- list(Ctrl_T1_CD_HFD.sig=Ctrl_T1_CD_HFD.sig, Ctrl_T2_CD_HFD.sig=Ctrl_T2_CD_HFD.sig,
            
                Tgt_T1_CD_HFD.sig=Tgt_T1_CD_HFD.sig, Tgt_T2_CD_HFD.sig=Tgt_T2_CD_HFD.sig,
                Tmx_T1_CD_HFD.sig=Tmx_T1_CD_HFD.sig, Tmx_T2_CD_HFD.sig=Tmx_T2_CD_HFD.sig,
                
                CD_T1_Ctrl_Tgt.sig=CD_T1_Ctrl_Tgt.sig, CD_T2_Ctrl_Tgt.sig=CD_T2_Ctrl_Tgt.sig,
                HFD_T1_Ctrl_Tgt.sig=HFD_T1_Ctrl_Tgt.sig, HFD_T2_Ctrl_Tgt.sig=HFD_T2_Ctrl_Tgt.sig,
                
                CD_T1_Ctrl_Tmx.sig=CD_T1_Ctrl_Tmx.sig, CD_T2_Ctrl_Tmx.sig=CD_T2_Ctrl_Tmx.sig,
                HFD_T1_Ctrl_Tmx.sig=HFD_T1_Ctrl_Tmx.sig, HFD_T2_Ctrl_Tmx.sig=HFD_T2_Ctrl_Tmx.sig)
sapply(all_sig, nrow)
write.xlsx(all_sig, "two_Drug_all_p05fc1d2fdr02_sig_01label.xlsx")

###venn diagram
sapply(all_sig, nrow)

#Tgt
vennplot(list(Ctrl_T1=Ctrl_T1_CD_HFD.sig$Name, T1_CD=CD_T1_Ctrl_Tgt.sig$Name, 
              Tgt_T1=Tgt_T1_CD_HFD.sig$Name, T1_HFD=HFD_T1_Ctrl_Tgt.sig$Name),"T1 Tgt", labelsize=1.5)

vennplot(list(Ctrl_T2=Ctrl_T2_CD_HFD.sig$Name, T2_CD=CD_T2_Ctrl_Tgt.sig$Name, 
              Tgt_T2=Tgt_T2_CD_HFD.sig$Name, T2_HFD=HFD_T2_Ctrl_Tgt.sig$Name),"T2 Tgt", labelsize=1.5)

#Tmx
vennplot(list(Ctrl_T1=Ctrl_T1_CD_HFD.sig$Name, T1_CD=CD_T1_Ctrl_Tmx.sig$Name, 
              Tmx_T1=Tmx_T1_CD_HFD.sig$Name, T1_HFD=HFD_T1_Ctrl_Tmx.sig$Name),"T1 Tmx", labelsize=1.5)

vennplot(list(Ctrl_T2=Ctrl_T2_CD_HFD.sig$Name, T2_CD=CD_T2_Ctrl_Tmx.sig$Name, 
              Tmx_T2=Tmx_T2_CD_HFD.sig$Name, T2_HFD=HFD_T2_Ctrl_Tmx.sig$Name),"T2 Tmx", labelsize=1.5)



vennplot(list(Ctrl_T1=Ctrl_T1_CD_HFD.sig$Name, Ctrl_T2=Ctrl_T2_CD_HFD.sig$Name, 
              Tgt_T1=Tgt_T1_CD_HFD.sig$Name, Tgt_T2=Tgt_T2_CD_HFD.sig$Name),"CD_HFD Tgt", labelsize=1.5)

vennplot(list(Ctrl_T1=Ctrl_T1_CD_HFD.sig$Name, Ctrl_T2=Ctrl_T2_CD_HFD.sig$Name, 
              Tmx_T1=Tmx_T1_CD_HFD.sig$Name, Tmx_T2=Tmx_T2_CD_HFD.sig$Name),"CD_HFD Tmx", labelsize=1.5)

vennplot(list(T1_CD=CD_T1_Ctrl_Tgt.sig$Name, T2_CD=CD_T2_Ctrl_Tgt.sig$Name, 
              T1_HFD=HFD_T1_Ctrl_Tgt.sig$Name, T2_HFD=HFD_T2_Ctrl_Tgt.sig$Name),"Tgt_Ctrl", labelsize=1.5)

vennplot(list(T1_CD=CD_T1_Ctrl_Tmx.sig$Name, T2_CD=CD_T2_Ctrl_Tmx.sig$Name, 
              T1_HFD=HFD_T1_Ctrl_Tmx.sig$Name, T2_HFD=HFD_T2_Ctrl_Tmx.sig$Name),"Tmx_Ctrl", labelsize=1.5)




##label significant metabolite with 0 or 1

Ctrl_T1_CD_HFD.sig.01 <- as.numeric(result3$Name %in% Ctrl_T1_CD_HFD.sig$Name)
Ctrl_T2_CD_HFD.sig.01 <- as.numeric(result3$Name %in% Ctrl_T2_CD_HFD.sig$Name)

Tgt_T1_CD_HFD.sig.01 <- as.numeric(result3$Name %in% Tgt_T1_CD_HFD.sig$Name)
Tgt_T2_CD_HFD.sig.01 <- as.numeric(result3$Name %in% Tgt_T2_CD_HFD.sig$Name)
Tmx_T1_CD_HFD.sig.01 <- as.numeric(result3$Name %in% Tmx_T1_CD_HFD.sig$Name)
Tmx_T2_CD_HFD.sig.01 <- as.numeric(result3$Name %in% Tmx_T2_CD_HFD.sig$Name)

CD_T1_Ctrl_Tgt.sig.01 <- as.numeric(result3$Name %in% CD_T1_Ctrl_Tgt.sig$Name)
CD_T2_Ctrl_Tgt.sig.01 <- as.numeric(result3$Name %in% CD_T2_Ctrl_Tgt.sig$Name)
HFD_T1_Ctrl_Tgt.sig.01 <- as.numeric(result3$Name %in% HFD_T1_Ctrl_Tgt.sig$Name)
HFD_T2_Ctrl_Tgt.sig.01 <- as.numeric(result3$Name %in% HFD_T2_Ctrl_Tgt.sig$Name)

CD_T1_Ctrl_Tmx.sig.01 <- as.numeric(result3$Name %in% CD_T1_Ctrl_Tmx.sig$Name)
CD_T2_Ctrl_Tmx.sig.01 <- as.numeric(result3$Name %in% CD_T2_Ctrl_Tmx.sig$Name)
HFD_T1_Ctrl_Tmx.sig.01 <- as.numeric(result3$Name %in% HFD_T1_Ctrl_Tmx.sig$Name)
HFD_T2_Ctrl_Tmx.sig.01 <- as.numeric(result3$Name %in% HFD_T2_Ctrl_Tmx.sig$Name)


label01 <- paste(Ctrl_T1_CD_HFD.sig.01, Ctrl_T2_CD_HFD.sig.01,
              
                 Tgt_T1_CD_HFD.sig.01, Tgt_T2_CD_HFD.sig.01,
                 Tmx_T1_CD_HFD.sig.01, Tmx_T2_CD_HFD.sig.01,
                
                 CD_T1_Ctrl_Tgt.sig.01, CD_T2_Ctrl_Tgt.sig.01,
                 HFD_T1_Ctrl_Tgt.sig.01, HFD_T2_Ctrl_Tgt.sig.01,
                 
                 CD_T1_Ctrl_Tmx.sig.01, CD_T2_Ctrl_Tmx.sig.01,
                 HFD_T1_Ctrl_Tmx.sig.01, HFD_T2_Ctrl_Tmx.sig.01,
                 sep="")

label01CDHFD <- paste(Ctrl_T1_CD_HFD.sig.01, Ctrl_T2_CD_HFD.sig.01,
                 Tgt_T1_CD_HFD.sig.01, Tgt_T2_CD_HFD.sig.01,
                 Tmx_T1_CD_HFD.sig.01, Tmx_T2_CD_HFD.sig.01,
                  sep="")

label01Tgt <- paste(CD_T1_Ctrl_Tgt.sig.01, CD_T2_Ctrl_Tgt.sig.01,
                 HFD_T1_Ctrl_Tgt.sig.01, HFD_T2_Ctrl_Tgt.sig.01,
                 sep="")

label01Tmx <- paste(CD_T1_Ctrl_Tmx.sig.01, CD_T2_Ctrl_Tmx.sig.01,
                 HFD_T1_Ctrl_Tmx.sig.01, HFD_T2_Ctrl_Tmx.sig.01,
                 sep="")


result3[,"sig"] <- label01
result3[,"CD_HFD"] <- label01CDHFD
result3[,"Tgt"] <- label01Tgt
result3[,"Tmx"] <- label01Tmx



write.xlsx(result3, "two_Drug_all_p_fc_fdr_sig_01label_020619.xlsx")

result3<- read.xlsx("two_Drug_all_p_fc_fdr_sig_01label_020619.xlsx")
names(result3)

names(result3)
Ctrl
Lpt
Mfm
Tgt
Tmx

Irs

names(result3)



#bar chart
all_sig_num <- sapply(all_sig, nrow)
write.xlsx(data.frame(all_sig_num), "two_drug_p05fc1d2fdr02_sig_summary_020619.xlsx",row.names=TRUE)

sig_num <- read.xlsx("C:/zhijuncao/sun/diet-drug/01282019/two_drug_p05fc1d2fdr02_sig_summary_020619.xlsx")
IncDec <- sig_num %>% gather(key="change", value="num", -c(1:8))
IncDec[,"ypos"] <- c(IncDec$num[IncDec$change=="Inc"], IncDec$all_sig_num[IncDec$change=="Dec"])

sig_num$Group <- factor(sig_num$Group, levels=c("Ctrl_T1", "Tgt_T1", "Tmx_T1", "Ctrl_T2", 
                                                "Tgt_T2", "Tmx_T2",  "CD_T1", "HFD_T1",  "CD_T2",   
                                                  "HFD_T2"))
tem <- IncDec

pdf("significant change metabolites increasedecrease_bar.pdf",width = 5, height =5)
pdf("significant change metabolites nummber_diet_and_treatment_bar_stack_2drug_meta568.pdf",width = 5, height =3)
pdf("significant change metabolites nummber_diet_and_treatment_bar_stack_2drug_meta565b.pdf",width = 5, height =3)
ggplot(IncDec[IncDec$Type=="treatment",], aes(Group,num, fill=change))+
  theme_bw()+
  geom_bar(stat = "identity",position="stack")+
  geom_text(aes(y=ypos, label=num), vjust=2, color="black", size=2.5,fontface = "bold")+
  labs(fill="HFD vs CD")+
  theme(axis.title=element_text(size=8),axis.text=element_text(size=8,color="black"))+
  theme(axis.text.x=element_text(angle=45,  hjust = 1))
#dev.off()

#pdf("significant change metabolites increasedecrease_treatment_bar_stack_2grug.pdf",width = 5, height =3)
ggplot(IncDec[IncDec$Type=="TimePoint",], aes(Group,num, fill=change))+
  theme_bw()+
  facet_grid(.~Compare)+
  geom_bar(stat = "identity",position="stack")+
  geom_text(aes(y=ypos, label=num), vjust=2, color="black", size=2.5,fontface = "bold")+
  labs(fill="Drug vs Ctrl")+
  theme(axis.title=element_text(size=8),axis.text=element_text(size=8,color="black"))+
  theme(axis.text.x=element_text(angle=45,  hjust = 1))
dev.off()



result3 <- read.xlsx("C:/zhijuncao/sun/diet-drug/01282019/two_drug_all_p_fc_fdr_01lable_02132019.xlsx")
###vocanoplot: c("analytes", "pvalue","ratio", "fdr")
result3[is.na(result3)] <- 1 #replace NA with 1

pdf("volcanoplot_020619_meta565.pdf", width = 8, height = 8)


#########
pfc <- result3[,c("BIOCHEMICAL", "Ctrl_T1_CD_HFD.p", "Ctrl_T1_CD_HFD.fc", "Ctrl_T1_CD_HFD.FDR")]
title <- "Ctrl_T1_CD_HFD"
volcanoplotfdr(pfc,title=title,fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)

pfc <- result3[,c("BIOCHEMICAL", "Ctrl_T2_CD_HFD.p", "Ctrl_T2_CD_HFD.fc", "Ctrl_T2_CD_HFD.FDR")]
title <- "Ctrl_T2_CD_HFD"
volcanoplotfdr(pfc,title=title,fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)

#########


#########
pfc <- result3[,c("BIOCHEMICAL", "Tgt_T1_CD_HFD.p", "Tgt_T1_CD_HFD.fc", "Tgt_T1_CD_HFD.FDR")]
title <- "Tgt_T1_CD_HFD"
volcanoplotfdr(pfc,title=title,fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)

pfc <- result3[,c("BIOCHEMICAL", "Tgt_T2_CD_HFD.p", "Tgt_T2_CD_HFD.fc", "Tgt_T2_CD_HFD.FDR")]
title <- "Tgt_T2_CD_HFD"
volcanoplotfdr(pfc,title=title,fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)

#########
pfc <- result3[,c("BIOCHEMICAL", "Tmx_T1_CD_HFD.p", "Tmx_T1_CD_HFD.fc", "Tmx_T1_CD_HFD.FDR")]
title <- "Tmx_T1_CD_HFD"
volcanoplotfdr(pfc,title=title,fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)

pfc <- result3[,c("BIOCHEMICAL", "Tmx_T2_CD_HFD.p", "Tmx_T2_CD_HFD.fc", "Tmx_T2_CD_HFD.FDR")]
title <- "Tmx_T2_CD_HFD"
volcanoplotfdr(pfc,title=title,fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)

##################################

##################################
pfc <- result3[,c("BIOCHEMICAL", "CD_T1_Ctrl_Tgt.p", "CD_T1_Ctrl_Tgt.fc", "CD_T1_Ctrl_Tgt.FDR")]
title <- "CD_T1_Ctrl_Tgt"
volcanoplotfdr(pfc,title=title,fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)

pfc <- result3[,c("BIOCHEMICAL", "CD_T2_Ctrl_Tgt.p", "CD_T2_Ctrl_Tgt.fc", "CD_T2_Ctrl_Tgt.FDR")]
title <- "CD_T2_Ctrl_Tgt"
volcanoplotfdr(pfc,title=title,fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)

pfc <- result3[,c("BIOCHEMICAL", "HFD_T1_Ctrl_Tgt.p", "HFD_T1_Ctrl_Tgt.fc", "HFD_T1_Ctrl_Tgt.FDR")]
title <- "HFD_T1_Ctrl_Tgt"
volcanoplotfdr(pfc,title=title,fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)

pfc <- result3[,c("BIOCHEMICAL", "HFD_T2_Ctrl_Tgt.p", "HFD_T2_Ctrl_Tgt.fc", "HFD_T2_Ctrl_Tgt.FDR")]
title <- "HFD_T2_Ctrl_Tgt"
volcanoplotfdr(pfc,title=title,fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)

##################################
pfc <- result3[,c("BIOCHEMICAL", "CD_T1_Ctrl_Tmx.p", "CD_T1_Ctrl_Tmx.fc", "CD_T1_Ctrl_Tmx.FDR")]
title <- "CD_T1_Ctrl_Tmx"
volcanoplotfdr(pfc,title=title,fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)

pfc <- result3[,c("BIOCHEMICAL", "CD_T2_Ctrl_Tmx.p", "CD_T2_Ctrl_Tmx.fc", "CD_T2_Ctrl_Tmx.FDR")]
title <- "CD_T2_Ctrl_Tmx"
volcanoplotfdr(pfc,title=title,fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)

pfc <- result3[,c("BIOCHEMICAL", "HFD_T1_Ctrl_Tmx.p", "HFD_T1_Ctrl_Tmx.fc", "HFD_T1_Ctrl_Tmx.FDR")]
title <- "HFD_T1_Ctrl_Tmx"
volcanoplotfdr(pfc,title=title,fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)

pfc <- result3[,c("BIOCHEMICAL", "HFD_T2_Ctrl_Tmx.p", "HFD_T2_Ctrl_Tmx.fc", "HFD_T2_Ctrl_Tmx.FDR")]
title <- "HFD_T2_Ctrl_Tmx"
volcanoplotfdr(pfc,title=title,fccutoff=1.2, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)

dev.off()
names(result3)

Ctrl
Lpt
Mfm
Tgt
Tmx

Irs



#PCA 
sample <- read.csv("C:/zhijuncao/sun/diet-drug/twodrug_scaledata_1_remove2samples_565metabolites.csv")
metaid <- read.csv("C:/zhijuncao/sun/diet-drug/metaID_1.csv")
drugs <- paste(metaid$Name[metaid$SUB_PATHWAY=="Drug"])

re <- read.xlsx("C:/zhijuncao/sun/diet-drug/01282019/two_drug_all_p_fc_fdr_01lable_02132019.xlsx")

sample<- sample[, !(names(sample)%in% drugs)]
dim(sample1)
dim(sample)
names(sample)[1:20]
dim(sample)
dim(tem)
logmeta <- log(sample[,-c(1:15)],2)
dim(logmeta)

length(names(logmeta))

library(mixOmics)
pca_analysis <- pca(logmeta, ncomp = 10)
 y1 <- paste(sample$DIET, sample$TIME2,sep="-")
 y2 <- sample$TREATMENT
plotIndiv(pca_analysis,ind.names = FALSE,group=y1, comp=c(1,2), title="pca_565",
          ellipse =TRUE, ellipse.level = 0.95,cex=3,legend=T,
          legend.title.pch = "Treatmet", pch = as.numeric(factor(y2)), pch.levels = factor(y2),size.legend=0.2)

plotIndiv(pca_analysis,ind.names = FALSE,group=y2, comp=c(1,3), title="pca_565",
          ellipse =F, ellipse.level = 0.95,cex=3,legend=T,
          legend.title.pch = "Treatmet", pch = as.numeric(factor(y2)), pch.levels = factor(y2),size.legend=0.2)





pcaplot(pca_analysis, cex=3, y1=paste(sample$DIET, sample$TIME2,sep="-"),y=sample$TREATMENT,ellipse=F, title ="PCA_twodrug_diet-time-treatment_565")

dev.off()
tem <- re$Name[re$sig!="00000000000000"]
pca_analysis <- pca(logmeta[,tem], ncomp = 10)
pcaplot(pca_analysis, cex=3, y1=paste(sample$DIET, sample$TIME2,sep="-"),y=sample$TREATMENT, title ="PCA_twodrug_diet-time-treatment_sig531")

tem <- re$Name[re$Tmx!="0000"]
pca_analysis <- pca(logmeta[,tem], ncomp = 10)
pcaplot(pca_analysis, cex=3, y1=paste(sample$DIET, sample$TIME2,sep="-"),y=sample$TREATMENT, title ="PCA_twodrug_diet-time-treatment_tmx314")

tem <- re$Name[re$Tgt!="0000"]
pca_analysis <- pca(logmeta[,tem], ncomp = 10)
pcaplot(pca_analysis, cex=3, y1=paste(sample$DIET, sample$TIME2,sep="-"),y=sample$TREATMENT, title ="PCA_twodrug_diet-time-treatment_tgt353")

tem <- re$Name[re$Tgt!="0000"|re$Tmx!="0000"]
pca_analysis <- pca(logmeta[,tem], ncomp = 10)
pcaplot(pca_analysis, cex=3, y1=paste(sample$DIET, sample$TIME2,sep="-"),y=sample$TREATMENT, title ="PCA_twodrug_diet-time-treatment_tgttmx457")

#splsda analysis
sample_tgt <- sample %>% filter(TREATMENT!="Tmx") %>% droplevels()
sample_tmx <- sample %>% filter(TREATMENT!="Tgt") %>% droplevels()

logmeta <- log(sample_tgt[,-c(1:15)],2)
dim(logmeta)
splsda_analysis <- splsda(logmeta ,sample_tgt$TREATMENT,ncomp =10, keepX = rep(565,10))
pcaplot(splsda_analysis, cex=3, y1=paste(sample_tgt$DIET, sample_tgt$TIME2,sep="-"),y=sample_tgt$TREATMENT, title ="splsda10x565_tgt")
#tgt_loading <- splsda_analysis$loadings
write.xlsx(tgt_loading,"tgt_loading.xlsx", row.names=TRUE)

dim(metaid)
pdf("splsda10x565_twodrug_diet-time-tgt_loading_top120.pdf",width = 9, height = 12)
plotLoadings(splsda_analysis, comp =1, size.name=0.5,contrib = 'max', 
             method = 'median', name.var = metaid$BIOCHEMICAL[metaid$Name %in% names(logmeta)], ndisplay=120, title=paste("comp:",1))
plotLoadings(splsda_analysis, comp =2, size.name=0.5,contrib = 'max', 
             method = 'median', name.var = metaid$BIOCHEMICAL[metaid$Name %in% names(logmeta)], ndisplay=120, title=paste("comp:",2))
dev.off()

###tmx
logmeta <- log(sample_tmx[,-c(1:15)],2)
logmeta1 <-logmeta %>% select(-targretin) 
dim(logmeta1)
names(logmeta)
splsda_analysis <- splsda(logmeta1 ,sample_tmx$TREATMENT,ncomp =10, keepX = rep(564,10))
pcaplot(splsda_analysis, cex=3, y1=paste(sample_tmx$DIET, sample_tmx$TIME2,sep="-"),y=sample_tmx$TREATMENT, title ="splsda10x564_tmx")

#tmx_loading <- splsda_analysis$loadings
write.xlsx(tmx_loading,"tmx_loading.xlsx", row.names=TRUE) 
dim(metaid)
pdf("splsda10x565_twodrug_diet-time-tmx_loading_top120.pdf",width = 9, height = 12)
plotLoadings(splsda_analysis, comp =1, size.name=0.5,contrib = 'max', 
             method = 'median', name.var = metaid$BIOCHEMICAL[metaid$Name %in% names(logmeta1)], ndisplay=120, title=paste("comp:",1))
plotLoadings(splsda_analysis, comp =2, size.name=0.5,contrib = 'max', 
             method = 'median', name.var = metaid$BIOCHEMICAL[metaid$Name %in% names(logmeta1)], ndisplay=120, title=paste("comp:",2))
dev.off()


tgt_loading_df<- tgt_loading$X %>% data.frame()%>% rownames_to_column(var = "Name")
tmx_loading_df<- tmx_loading$X %>% data.frame()%>% rownames_to_column(var = "Name")
loading_df <- left_join(tgt_loading_df[,1:2],tmx_loading_df[,1:2], by="Name", suffix = c(".Tgt", ".Tmx"))
str(loading_df)
dim(loading_df)
re_loadint <- left_join(re,loading_df) 
write.xlsx(re_loadint, "C:/zhijuncao/sun/diet-drug/03182019/twodrug_diet-time-treatment_loading_with_stat_result.xlsx", rowNames=TRUE)

dev.off()
?plotLoadings()
re$Name



pcaplot(pca_analysis, cex=3, y1=paste(sample$DIET, sample$TIME2,sep="-"),y=sample$TREATMENT, title ="PCA_twodrug_diet-time-treatment_565")
plotVar(pca_analysis, comp = c(1,2), cutoff = 0.9, cex=3)
plotVar(pca_analysis, comp = c(1,3), cutoff = 0.8, cex=1)
pcaplot(pca_analysis, cex=3, y1=paste(sample$DIET, sample$TIME2,sep="-"),y=sample$TREATMENT, title ="PCA_twodrug_diet-time-treatment_nodrug")
plotIndiv(pca_analysis, cex=2, group=sample$TUMOR, ind.names = sample$SAMPLE_NAME, title ="name_time-diet-tumor_pca")
pcaplot(pca_analysis, cex=2, y1=sample$GROUP_DESCRIPTION,y=sample$TUMOR, title ="twodrug_group-tumor_pca")
dev.off()
?plotLoadings()
t1_cd <- cbind(sample[,1:15],log(sample[,-c(1:15)],2)) %>% filter(DIET=="control diet",TIME2=="T1") %>% droplevels()
t1_hfd <- cbind(sample[,1:15],log(sample[,-c(1:15)],2)) %>% filter(DIET=="high fat", TIME2=="T1")%>% droplevels()

t2_cd <- cbind(sample[,1:15],log(sample[,-c(1:15)],2)) %>% filter(DIET=="control diet",TIME2=="T2")%>% droplevels()
t2_hfd <- cbind(sample[,1:15],log(sample[,-c(1:15)],2)) %>% filter(DIET=="high fat",TIME2=="T2")%>% droplevels()

tem <-t1_cd
title1 <- "twodrug_t1_cd"

tem <-t1_hfd
title1 <- "twodrug_t1_hfd"

tem <-t2_cd
title1 <- "twodrug_t2_cd"

tem <-t2_hfd
title1 <- "twodrug_t2_hfd"

pca_analysis <- pca(tem[,-c(1:15)], ncomp = 10)
#pcaplot(pca_analysis, cex=2, y1=paste(tem$TIME2,tem$DIET,sep="-"),y=tem$TUMOR, title =paste(title1,"time2-diet-tumor_pca"))
pcaplot(pca_analysis, cex=2, y1=tem$GROUP_DESCRIPTION,y=tem$TUMOR, title =paste(title1,"group-tumor_pca"))



#boxplot
sample[,"Diet_Time"] <- NULL
twodrug <- read.csv("C:/zhijuncao/sun/diet-drug/twodrug_scaledata_1_remove2samples_568metabolites.csv")
metaid <- read.csv("C:/zhijuncao/sun/diet-drug/metaID_1.csv")
sample <- twodrug
dim(sample)
pdf("C:/zhijuncao/sun/diet-drug/01282019/Diet_tumor_scaledImp_time_diet_group_boxplot_all_pathway.pdf", width = 9, height = 3)
for (name in metaid$Name){
  #name <- "pinitol"
  xtitle <- paste('Name:',metaid$BIOCHEMICAL[metaid$Name==name],'\n','Super Pathway:',metaid$SUPER_PATHWAY[metaid$Name==name],'\n','Sub Pathway:',metaid$SUB_PATHWAY[metaid$Name==name])
p <- ggplot(sample, aes(TIME2, log(eval(parse(text = name)),2), fill=DIET))+
  theme_bw()+
  facet_wrap(GROUP_DESCRIPTION~.,nrow=1)+
  geom_boxplot()+
  labs(y="log2 ScaledImp", x=xtitle)+
  theme(axis.title=element_text(size=8),axis.text=element_text(size=8,color="black"))+
  theme(axis.text.x=element_text(angle=0,  hjust = 1))+
  theme(strip.background = element_rect(fill = NULL, colour = NA))+
  theme(strip.text = element_text(colour = "black", size = rel(0.8)))

print (p)
  
}
dev.off()
i <- 1
dim(twodrug)
### 2 drug boxplots
pdf("C:/zhijuncao/sun/diet-drug/12142018/Boxplot_twodrug_scaledImp_time_diet_treatment_1.pd.pdf", width = 6, height = 3)
for (name in metaid$Name){
  #name <- "pinitol"
  xtitle <- paste('Name:',metaid$BIOCHEMICAL[metaid$Name==name],'\n','Super Pathway:',metaid$SUPER_PATHWAY[metaid$Name==name],'\n','Sub Pathway:',metaid$SUB_PATHWAY[metaid$Name==name])
  p <- ggplot(twodrug, aes(TIME2, log(eval(parse(text = name)),2), fill=DIET))+
    theme_bw()+
    facet_wrap(TREATMENT~.,nrow=1)+
    geom_boxplot()+
    labs(y="log2 ScaledImp", x=xtitle)+
    theme(axis.title=element_text(size=8),axis.text=element_text(size=8,color="black"))+
    theme(axis.text.x=element_text(angle=0,  hjust = 1))+
    theme(strip.background = element_rect(fill = NULL, colour = NA))+
    theme(strip.text = element_text(colour = "black", size = rel(0.8)))
  
  print (p)
  cat(i)
  i <- i+1
}
dev.off()

i <- 1
pdf("C:/zhijuncao/sun/diet-drug/112018/twodrug_Diet_tumor_scaledImp_diet_group_time_boxplot.pdf", width = 6, height = 3)
pdf("C:/zhijuncao/sun/diet-drug/12142018/Boxplot_twodrug_scaledImp_diet_treatment_time_1.pdf", width = 6, height = 3)
for (name in metaid$Name){
  #name <- "pinitol"
  xtitle <- paste('Name:',metaid$BIOCHEMICAL[metaid$Name==name],'\n','Super Pathway:',metaid$SUPER_PATHWAY[metaid$Name==name],'\n','Sub Pathway:',metaid$SUB_PATHWAY[metaid$Name==name])
  p <- ggplot(twodrug, aes(DIET, log(eval(parse(text = name)),2), fill=TREATMENT))+
    theme_bw()+
    facet_wrap(TIME2~.,nrow=1)+
    geom_boxplot()+
    labs(y="log2 ScaledImp", x=xtitle, fill="Treatment")+
    theme(axis.title=element_text(size=8),axis.text=element_text(size=8,color="black"))+
    theme(axis.text.x=element_text(angle=0,  hjust = 1))+
    theme(strip.background = element_rect(fill = NULL, colour = NA))+
    theme(strip.text = element_text(colour = "black", size = rel(0.8)))
  
  print (p)
  cat(i)
  i <- i+1
}
dev.off()





i <- 0
pdf("C:/zhijuncao/sun/diet-drug/01282019/Diet_tumor_scaledImp_diet_group_time_boxplot_all_pathway_6x3.pdf", width = 6, height = 3)
for (name in metaid$Name){
  #name <- "pinitol"
  xtitle <- paste('Name:',metaid$BIOCHEMICAL[metaid$Name==name],'\n','Super Pathway:',metaid$SUPER_PATHWAY[metaid$Name==name],'\n','Sub Pathway:',metaid$SUB_PATHWAY[metaid$Name==name])
  p <- ggplot(sample, aes(DIET, log(eval(parse(text = name)),2), fill=TREATMENT))+
    theme_bw()+
    facet_wrap(TIME2~.,nrow=1)+
    geom_boxplot()+
    labs(y="log2 ScaledImp", x=xtitle, fill="Treatment")+
    theme(axis.title=element_text(size=8),axis.text=element_text(size=8,color="black"))+
    theme(axis.text.x=element_text(angle=0,  hjust = 1))+
    theme(strip.background = element_rect(fill = NULL, colour = NA))+
    theme(strip.text = element_text(colour = "black", size = rel(0.8)))
  
  print (p)
  cat(i)
  i <- i+1
  
}
dev.off()


i <- 0
pdf("C:/zhijuncao/sun/diet-drug/01282019/Diet_tumor_scaledImp_group_diet_time_boxplot_all_pathway_6x3.pdf", width = 6, height = 3)
for (name in metaid$Name){
  #name <- "pinitol"
  xtitle <- paste('Name:',metaid$BIOCHEMICAL[metaid$Name==name],'\n','Super Pathway:',metaid$SUPER_PATHWAY[metaid$Name==name],'\n','Sub Pathway:',metaid$SUB_PATHWAY[metaid$Name==name])
  p <- ggplot(sample, aes(TREATMENT, log(eval(parse(text = name)),2), fill=DIET))+
    theme_bw()+
    facet_wrap(TIME2~.,nrow=1)+
    geom_boxplot()+
    labs(y="log2 ScaledImp", x=xtitle, fill="Diet")+
    theme(axis.title=element_text(size=8),axis.text=element_text(size=8,color="black"))+
    theme(axis.text.x=element_text(angle=0,  hjust = 1))+
    theme(strip.background = element_rect(fill = NULL, colour = NA))+
    theme(strip.text = element_text(colour = "black", size = rel(0.8)))
  
  print (p)
  cat(i)
  i <- i+1
  
}
dev.off()





#tumor

pdf("C:/zhijuncao/sun/diet-drug/Diet_tumor_scaledImp_tumor_boxplot_all_pathway.pdf", width = 5, height = 3)
for (name in metaid$Name){
  #name <- "pinitol"
  xtitle <- paste('Name:',metaid$BIOCHEMICAL[metaid$Name==name],'\n','Super Pathway:',metaid$SUPER_PATHWAY[metaid$Name==name],'\n','Sub Pathway:',metaid$SUB_PATHWAY[metaid$Name==name])
  p <- ggplot(sample, aes(DIET, log(eval(parse(text = name)),2), fill=TUMOR))+
    theme_bw()+
    facet_wrap(TIME2~.)+
    geom_boxplot()+
    labs(y="log2 ScaledImp", x=xtitle)+
    theme(axis.title=element_text(size=8),axis.text=element_text(size=8,color="black"))+
    theme(axis.text.x=element_text(angle=45,  hjust = 1))+
    theme(strip.background = element_rect(fill = NULL, colour = NA))+
    theme(strip.text = element_text(colour = "black", size = rel(0.8)))
  
  print (p)
  
}
dev.off()



names(metaid)












