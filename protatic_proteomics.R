
library(openxlsx)
library(mixOmics)
library(tidyr)
library(dplyr)
library(magrittr)
library(tibble)
library(gridExtra)
library(ggplot2)
library(plotly)
library(purrr)
library(multcomp)
library(stringr)


re <- read.xlsx("C:/zhijuncao/nori/proteomics/12072020/Supplemental Table 3_12032020_for manuscript_R.xlsx")
names(re)
re1 <- re %>% select(protein_name,`stat_status(90EB,100EB,100EBTE,145EB,145EBTE,200EB,200EBTE)`) %>% 
  separate(2, into = c("pnd90EB","pnd100EB","pnd100EBTE","pnd145EB","pnd145EBTE","pnd200EB","pnd200EBTE"), sep = c(1,2,3,4,5,6))

pnd100EBTE <- re1 %>% filter(pnd100EBTE==1)
pnd145EBTE <- re1 %>% filter(pnd145EBTE==1)
pnd200EBTE <- re1 %>% filter(pnd200EBTE==1)

vennplot(list(PND100EBTE=pnd100EBTE$protein_name, PND145EBTE=pnd145EBTE$protein_name, PND200EBTE=pnd200EBTE$protein_name),
         "C:/zhijuncao/nori/proteomics/12072020/venn_EBTE.tiff", labelsize=0.7)

names(re1)

rm(list = ls())
protein <- read.xlsx("C:/zhijuncao/nori/proteomics/result/protein_groups_iBAQ_2000_user.xlsx")
proteinid <- read.xlsx("C:/zhijuncao/nori/proteomics/result/protein_groups_iBAQ_2000_user.xlsx", sheet = 'proteinid')
dim(proteinid)
proteinid <- proteinid %>% filter(contaminant==0, nest_under==0) #remove contaminant and nest under proteins

protein1 <- cbind(protein[, c(1:4)], protein[,proteinid$protein_name])

rmNA <- function(pnd="PND 90", nonNA=4){
 tem <- protein1 %>% filter(PND==pnd)
 sumNA <- tem %>% map_df(~sum(is.na(.)))
 return (tem[,(nrow(tem)-sumNA)>=nonNA] ) 
}


pnd90 <- rmNA(pnd="PND 90", nonNA=4)
pnd100 <- rmNA(pnd="PND 100", nonNA=6)
pnd145 <- rmNA(pnd="PND 145", nonNA=6)
pnd200<- rmNA(pnd="PND 200", nonNA=6)

imputationNA <- function(data=protein){
  logdata <- log(data[,-c(1:4)],2)
  
  # logdata <- logdata[!is.na(apply(logdata , 1, function(x) var(x, na.rm = TRUE))),
  #                                !is.na(apply(logdata , 2, function(x) var(x, na.rm = TRUE)))]
  
  nipals.data = nipals(logdata, reconst = TRUE, ncomp = 8)$rec
  
  # only replace the imputation for the missing values
  id.na <- is.na(logdata)
  nipals.data[!id.na] <- logdata[!id.na]
  df <- cbind(data[,c(1:4)],nipals.data)
  return (df)
}

repNAwithmin <- function(x){
  if (sum(is.na(x))!=length(x))
  x[is.na(x)] <- min(x,na.rm = TRUE)
  return (x)
}

repNAwithmin05 <- function(x){
  if (sum(is.na(x))!=length(x))
    x[is.na(x)] <- min(x,na.rm = TRUE)*0.5
  return (x)
}

pnd4more <- protein1[, unique(c(names(pnd90),names(pnd100),names(pnd145),names(pnd200)))]

write.xlsx(list(pnd=pnd4more, pnd90=pnd90,pnd100=pnd100,pnd145=pnd145,pnd200=pnd200), "half or less missing for each pnd.xlsx")

names(pnd4more)

pnd90_NAmin05 <- pnd90 %>% group_by(Treatment)%>% mutate_at(-c(1:4), repNAwithmin) %>% ungroup() %>% mutate_at(-c(1:4), repNAwithmin05)%>% data.frame()
pnd100_NAmin05 <- pnd100 %>% group_by(Treatment)%>% mutate_at(-c(1:4), repNAwithmin) %>% ungroup() %>% mutate_at(-c(1:4), repNAwithmin05)%>%  data.frame()
pnd145_NAmin05 <- pnd145 %>% group_by(Treatment)%>% mutate_at(-c(1:4), repNAwithmin) %>% ungroup() %>% mutate_at(-c(1:4), repNAwithmin05)%>% data.frame()
pnd200_NAmin05 <- pnd200 %>% group_by(Treatment)%>% mutate_at(-c(1:4), repNAwithmin) %>% ungroup() %>% mutate_at(-c(1:4), repNAwithmin05)%>% data.frame()

imp1 <- merge(pnd90_NAmin05,pnd100_NAmin05, all = TRUE)
imp2 <- merge(pnd145_NAmin05 ,pnd200_NAmin05, all = TRUE)
imp <- merge(imp1,imp2, all=TRUE)

dim(imp)
# protein_median <- apply(imp[,-c(1:4)],1, function(x)median(x, na.rm = TRUE))
# normalizeValue <- median(protein_median)/protein_median

protein_sum <- apply(imp[,-c(1:4)],1, function(x)sum(x, na.rm = TRUE))
sumnormalizeValue <- median(protein_sum)/protein_sum
normalize <- function(x)x*sumnormalizeValue

imp_sumnormalized_log2 <- imp %>% mutate_at(-c(1:4), normalize) %>% mutate_at(-c(1:4), log2)
write.xlsx(imp_sumnormalized_log2, "half or less missing for each pnd NAminmin05_normalized_log2.xlsx")


imp_sumnormalized_log2a<- read.xlsx("C:/zhijuncao/nori/proteomics/result/nest0repminmin05normlog/half or less missing for each pnd NAminmin05_normalized_log2a.xlsx")

library(pheatmap)
anno_col <- imp_sumnormalized_log2a[, c("Treatment", "PND")]
anno_col$PND <- factor(anno_col$PND, level=c(90,100, 145, 200))
tem <- t(imp_sumnormalized_log2a[,-c(1:4)])
colnames(tem) <- rownames(anno_col)
pdf("C:/zhijuncao/nori/proteomics/result/nest0repminmin05normlog/heatmap6x4.pdf", width = 6, height = 4)
pheatmap(tem,scale="row",
         annotation_col=anno_col,
         show_rownames = T,fontsize_row=1,
         show_colnames = F,gaps_col = c(8,20,32 ),
         cluster_rows=F,cluster_cols=F)
dev.off()

pnd90_NAmin05log2_normalized <- imp_sumnormalized_log2 %>% filter(PND=="PND 90") %>% select(all_of(names(pnd90_NAmin05)))
pnd100_NAmin05log2_normalized <- imp_sumnormalized_log2 %>% filter(PND=="PND 100") %>% select(all_of(names(pnd100_NAmin05)))
pnd145_NAmin05log2_normalized <- imp_sumnormalized_log2 %>% filter(PND=="PND 145") %>% select(all_of(names(pnd145_NAmin05)))
pnd200_NAmin05log2_normalized <- imp_sumnormalized_log2 %>% filter(PND=="PND 200") %>% select(all_of(names(pnd200_NAmin05)))
pnd90_NAmin05log2_normalized$Treatment

pnd90_log2 <- pnd90 %>% mutate_at((names(pnd90)[-c(1:4)]), log2)
pnd100_log2 <- pnd100 %>% mutate_at((names(pnd100)[-c(1:4)]), log2)
pnd145_log2 <- pnd145 %>% mutate_at((names(pnd145)[-c(1:4)]), log2)
pnd200_log2 <- pnd200 %>% mutate_at((names(pnd200)[-c(1:4)]), log2)

pnd_log2a <- merge(pnd90_log2,pnd100_log2, all = TRUE)
pnd_log2b <- merge(pnd145_log2, pnd200_log2, all = TRUE)
pnd_log2 <- merge(pnd_log2a,pnd_log2b, all=TRUE)

write.xlsx(list(pnd_log2=pnd_log2,
                pnd90_log2=pnd90_log2,
                pnd100_log2=pnd100_log2,
                pnd145_log2=pnd145_log2,
                pnd200_log2=pnd200_log2), "half missing for each pndlog2.xlsx")


###pca function
pca_analysis <- function(title="PND200", tem=pnd200, colors=colors){
  
  tem_pca <- pca(tem[,-(1:4)], scale = TRUE, ncomp = 3)
  loading <- data.frame(tem_pca$loadings$X)
  absloading <- loading %>% mutate_all(abs)
  names(absloading) <- paste("abs_",names(absloading),sep="")
  loading_abs <- cbind(loading, absloading) %>% rownames_to_column ("protein_name")
  loading_abs <- inner_join(proteinid,loading_abs) 
  write.xlsx(list(explained=data.frame(tem_pca$explained_variance), loading=loading_abs), paste(title," pca results.xlsx"), rowNames=T)  

  explained <-paste(format(100*tem_pca$explained_variance, 3,2),"%", sep="")
  #3d pca
  
  p <- plot_ly(data.frame(tem_pca$variates$X), x = ~PC1, y = ~PC2, z = ~PC3, color = ~tem$Treatment, colors = colors) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = paste('PC 1:', explained[1])),
                        yaxis = list(title = paste('PC 2:', explained[2])),
                        zaxis = list(title = paste('PC 3:',  explained[3]))))
  
  
  htmlwidgets::saveWidget(as_widget(p), paste(title," pca_3d.html", sep=""))
  
  pdf(paste(title, "pca_2d.pdf"), width=8, height=6)
  plotIndiv(tem_pca, comp=c(1,2), group=tem$Treatment, col.per.group=colors, legend = TRUE, ind.names=F, title = title)
  plotIndiv(tem_pca, comp=c(1,3), group=tem$Treatment, col.per.group=colors, legend = TRUE, ind.names=F, title = title)
  plotIndiv(tem_pca, comp=c(2,3), group=tem$Treatment, col.per.group=colors, legend = TRUE, ind.names=F, title = title)
  dev.off()
  
  #variable correlation and loading
  pdf(paste(title, "pca_correlation and loading.pdf"),width=12, height=12)
  plotVar(tem_pca, comp=c(1,2), cex=4, title = title)
  plotVar(tem_pca, comp=c(1,3), cex=4, title = title)
  plotVar(tem_pca, comp=c(2,3), cex=4, title = title)
  
  plotLoadings(tem_pca,comp=1, size.name = 0.3, ndisplay = 100, title=paste(title,": contribution on comp 1"))
  plotLoadings(tem_pca,comp=2, size.name = 0.3, ndisplay = 100, title=paste(title,": contribution on comp 2"))
  plotLoadings(tem_pca,comp=3, size.name = 0.3, ndisplay = 100, title=paste(title,": contribution on comp 3"))
  
  dev.off()
  
}

###plsda function
plsda_analysis <- function(title="PND30", tem=pnd30, colors=colors){
  
  tem_plsda <- plsda(tem[,-(1:4)],tem$Treatment, scale = TRUE, ncomp = 3)
  perf_plsda <- perf(tem_plsda, validation = "loo", progressBar = FALSE, auc = TRUE)
  
  loading <- data.frame(tem_plsda$loadings$X)
  absloading <- loading %>% mutate_all(abs)
  names(absloading) <- paste("abs_",names(absloading),sep="")
  loading_abs <- cbind(loading, absloading) %>% rownames_to_column ("protein_name")
  loading_abs <- inner_join(proteinid,loading_abs) 
  
  write.xlsx(list(explained=data.frame(tem_plsda$explained_variance),
                  validatonerror=perf_plsda$error.rate$BER,
                  loading=loading_abs), paste(title," plsda results.xlsx"), rowNames=T)  
  
  explained <-paste(format(100*tem_plsda$explained_variance$X, 3,2),"%", sep="")
  #3d plsda
  
  p <- plot_ly(data.frame(tem_plsda$variates$X), x = ~comp1, y = ~comp2, z = ~comp3, color = ~tem$Treatment, colors = colors) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = paste('Xâvariate 1:', explained[1])),
                        yaxis = list(title = paste('Xâvariate 2:', explained[2])),
                        zaxis = list(title = paste('Xâvariate 3:',  explained[3]))))
  
  htmlwidgets::saveWidget(as_widget(p), paste(title," plsda_3d.html", sep=""))
  
  #variable correlation and loading
  pdf(paste(title, "plsda_2d.pdf"), width=8, height=6)
  plotIndiv(tem_plsda, comp=c(1,2), group=tem$Treatment,col.per.group=colors, legend = TRUE, ind.names=F, title = title)
  plotIndiv(tem_plsda, comp=c(1,3), group=tem$Treatment,col.per.group=colors, legend = TRUE, ind.names=F, title = title)
  plotIndiv(tem_plsda, comp=c(2,3), group=tem$Treatment,col.per.group=colors, legend = TRUE, ind.names=F, title = title)
  dev.off()
  
  pdf(paste(title, "plsda_correlation and loading.pdf"), width=12, height=12)
  plotVar(tem_plsda, comp=c(1,2), cex=4, title = title)
  plotVar(tem_plsda, comp=c(1,3), cex=4, title = title)
  plotVar(tem_plsda, comp=c(2,3), cex=4, title = title)
  
  plotLoadings(tem_plsda,comp=1,method='median', contrib="max", size.legend = 0.8, size.name = 0.3, ndisplay = 100, title=paste(title,": contribution on comp 1"))
  plotLoadings(tem_plsda,comp=2,method='median', contrib="max", size.legend = 0.8, size.name = 0.3, ndisplay = 100, title=paste(title,": contribution on comp 2") )
  plotLoadings(tem_plsda,comp=3,method='median', contrib="max", size.legend = 0.8, size.name = 0.3, ndisplay = 100, title=paste(title,": contribution on comp 3") )
  dev.off()
}

colors <- c('black', 'red')
plsda_analysis(title="PND90", tem=pnd90_NAmin05log2_normalized, colors=colors)
pca_analysis(title="PND90", tem=pnd90_NAmin05log2_normalized, colors=colors)


colors <- c('black', 'red', 'purple')
plsda_analysis(title="PND100", tem=pnd100_NAmin05log2_normalized, colors=colors)
pca_analysis(title="PND100", tem=pnd100_NAmin05log2_normalized, colors=colors)


colors <- c('black', 'red', 'purple')
plsda_analysis(title="PND145", tem=pnd145_NAmin05log2_normalized, colors=colors)
pca_analysis(title="PND145", tem=pnd145_NAmin05log2_normalized, colors=colors)

colors <- c('black', 'red', 'purple')
plsda_analysis(title="PND200", tem=pnd200_NAmin05log2_normalized, colors=colors)
pca_analysis(title="PND200", tem=pnd200_NAmin05log2_normalized, colors=colors)

#dev.off()

###wrap text
wrapper <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}
###point errorbar
protein1$PND <- factor(protein1$PND, levels = c('PND 90', 'PND 100', 'PND 145', 'PND 200'))
  #analytes <- unique(c(names(pnd90)[-c(1:4)], names(pnd100)[-c(1:4)],names(pnd145)[-c(1:4)],names(pnd200)[-c(1:4)]))
  analytes <- names(imp)[-c(1:4)]
  #analyte <- analytes[1:2]
  pdf("protein points errorbar.pdf", width = 4, height = 3)
  for (analyte in analytes){
    p <- ggplot(protein1, aes(PND, log(get(analyte),2),group=Treatment, color=Treatment))+
      theme_bw()+ 
      #facet_wrap(.~PND, nrow = 1)+
      #geom_boxplot()+
      geom_point(position=position_dodge(0.6),size=0.7,alpha=0.4)+
      scale_color_manual(values=c('black', 'red', 'purple'))+
      labs(y=analyte, x="", title = wrapper(paste(proteinid$description[proteinid$protein_name==analyte],
                                                  ':',
                                            proteinid$gene_name[proteinid$protein_name==analyte]),40))+
      theme(axis.title=element_text(size=10),axis.text=element_text(size=10),axis.text.x=element_text(angle=90))+
      theme(plot.title=element_text(size=10))+
      #theme(plot.margin=unit(c(1,1,1,1),"line"))+

      stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                     geom="errorbar",  width=0.4, position=position_dodge(0.6)) +
      stat_summary(fun=mean, geom="point", shape=4, position=position_dodge(0.6))
    
    print(p)
  }
  dev.off()

  imp_sumnormalized_log2$PND <- factor(imp_sumnormalized_log2$PND, levels = c('PND 90', 'PND 100', 'PND 145', 'PND 200'))
  ###point errorbar
  analytes <- names(imp_sumnormalized_log2)[-c(1:4)]
  pdf("protein imputation normalized_log2 points errorbar.pdf", width = 4, height = 3)
  for (analyte in analytes){
    p <- ggplot(imp_sumnormalized_log2, aes(PND, get(analyte),group=Treatment, color=Treatment))+
      theme_bw()+ 
      #facet_wrap(.~PND, nrow = 1)+
      #geom_boxplot()+
      geom_point(position=position_dodge(0.6),size=0.7,alpha=0.4)+
      scale_color_manual(values=c('black', 'red', 'purple'))+
      labs(y=analyte, x="", title = wrapper(paste(proteinid$description[proteinid$protein_name==analyte],
                                                  ':',
                                                  proteinid$gene_name[proteinid$protein_name==analyte]),40))+
      theme(axis.title=element_text(size=10),axis.text=element_text(size=10),axis.text.x=element_text(angle=90))+
      theme(plot.title=element_text(size=10))+
      #theme(plot.margin=unit(c(1,1,1,1),"line"))+
      
      stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                   geom="errorbar",  width=0.4, position=position_dodge(0.6)) +
      stat_summary(fun.y=mean, geom="point", shape=4, position=position_dodge(0.6))
    
    print(p)
  }
  dev.off()
  
  ###Before after NA inputation and normalization
  
  analytes <- names(imp_sumnormalized_log2)[-c(1:4)]
  pdf("Before and afterprotein imputation normalized_log2 points errorbar.pdf", width = 8, height = 3)
  pdf("Before and afterprotein imputation normalized_log2 significantly changed points errorbar.pdf", width = 8, height = 3)
  for (analyte in analytes){
    p1 <- ggplot(protein1, aes(PND, log(get(analyte),2),group=Treatment, color=Treatment))+
      theme_bw()+ 
      #facet_wrap(.~PND, nrow = 1)+
      #geom_boxplot()+
      geom_point(position=position_dodge(0.6),size=0.7,alpha=0.4)+
      scale_color_manual(values=c('black', 'red', 'purple'))+
      labs(y=analyte, x="", title = wrapper(paste(proteinid$description[proteinid$protein_name==analyte],
                                                  ':',
                                                  proteinid$gene_name[proteinid$protein_name==analyte]),40))+
      theme(axis.title=element_text(size=10),axis.text=element_text(size=10),axis.text.x=element_text(angle=90))+
      theme(plot.title=element_text(size=10))+
      #theme(plot.margin=unit(c(1,1,1,1),"line"))+
      
      stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                   geom="errorbar",  width=0.4, position=position_dodge(0.6)) +
      stat_summary(fun=mean, geom="point", shape=4, position=position_dodge(0.6))
    
    p2 <- ggplot(imp_sumnormalized_log2, aes(PND, get(analyte),group=Treatment, color=Treatment))+
      theme_bw()+ 
      #facet_wrap(.~PND, nrow = 1)+
      #geom_boxplot()+
      geom_point(position=position_dodge(0.6),size=0.7,alpha=0.4)+
      scale_color_manual(values=c('black', 'red', 'purple'))+
      labs(y=analyte, x="", title = wrapper(paste(proteinid$description[proteinid$protein_name==analyte],
                                                  ':',
                                                  proteinid$gene_name[proteinid$protein_name==analyte]),40))+
      theme(axis.title=element_text(size=10),axis.text=element_text(size=10),axis.text.x=element_text(angle=90))+
      theme(plot.title=element_text(size=10))+
      #theme(plot.margin=unit(c(1,1,1,1),"line"))+
      
      stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                   geom="errorbar",  width=0.4, position=position_dodge(0.6)) +
      stat_summary(fun.y=mean, geom="point", shape=4, position=position_dodge(0.6))
    
    grid.arrange(p1, p2, nrow=1)
  }
  dev.off()
  
  ## Statistical analysis

  aov_dunn <- function(data=pnd100_log2imp){
    
    Mean <- function(x)mean(x, na.rm=TRUE)
    exp2 <- function(x)2^x
    ratio <- function(x)(x/x[1])
    
    data$Treatment <- factor(data$Treatment)
    analytes <- names(data)[-c(1:4)]
 
      ratio_re <- data %>% select(-animalID,)%>%mutate_if(is.numeric, exp2) %>% group_by(Treatment) %>% summarise_if(is.numeric,Mean) %>% 
        mutate_if(is.numeric, ratio)
      ratio_re_df <- data.frame(t(as.matrix(ratio_re[,-1]))) %>% select(-1)
      
      pvalues <- list()
      
      for (x in analytes){
        m1 <- aov(data[,x] ~ Treatment, data=data)
        re <- glht(m1, linfct=mcp(Treatment="Dunnett"))
        sum_re_pvalue <- summary(re)$test$pvalues
        names(sum_re_pvalue) <- names(coefficients(re))
        pvalues[[x]] <- sum_re_pvalue
      }
      
      df_re <- data.frame(do.call(rbind, pvalues))
      fdr_re <- map_df(df_re, ~p.adjust(.x, method="BH"))
      names(fdr_re) <- paste(names(df_re), " FDR")
      names(ratio_re_df) <- paste(names(df_re), " Ratio")
      ratio_re_df <- ratio_re_df %>% rownames_to_column(var="analytes")
      names(df_re) <- paste(names(df_re), "p value")
      pnd_re <- cbind(ratio_re_df, df_re, fdr_re)
      pnd_re1 <- inner_join(proteinid, pnd_re, by=c('protein_name'='analytes'))
      return(pnd_re1)
  }
 

  pnd90_re <- aov_dunn(data=pnd90_NAmin05log2_normalized)
  pnd100_re <- aov_dunn(data=pnd100_NAmin05log2_normalized)
  pnd145_re <- aov_dunn(data=pnd145_NAmin05log2_normalized)
  pnd200_re <- aov_dunn(data=pnd200_NAmin05log2_normalized)
  
  class(pnd90_NAmin05log2)
  pnd90_rep005fc1d5fdr02 <-  pnd90_re %>% filter(`EB...CTRL p value`<0.05, `EB...CTRL  FDR`<0.2, `EB...CTRL  Ratio`>=1.5|`EB...CTRL  Ratio`<=1/1.5)
 
  
  pnd100_rep005fc1d5fdr02 <-  pnd100_re %>% filter((`EB...CTRL p value`<0.05 & `EB...CTRL  FDR`<0.2 & (`EB...CTRL  Ratio`>=1.5|`EB...CTRL  Ratio`<=1/1.5))
                                                  |(`EB.T.E...CTRL p value`<0.05 & `EB.T.E...CTRL  FDR`<0.2 & (`EB.T.E...CTRL  Ratio`>=1.5|`EB.T.E...CTRL  Ratio`<=1/1.5)))
  
  pnd145_rep005fc1d5fdr02 <-  pnd145_re %>% filter((`EB...CTRL p value`<0.05 & `EB...CTRL  FDR`<0.2 & (`EB...CTRL  Ratio`>=1.5|`EB...CTRL  Ratio`<=1/1.5))
                                                   |(`EB.T.E...CTRL p value`<0.05 & `EB.T.E...CTRL  FDR`<0.2 & (`EB.T.E...CTRL  Ratio`>=1.5|`EB.T.E...CTRL  Ratio`<=1/1.5)))
  
  pnd200_rep005fc1d5fdr02 <-  pnd200_re %>% filter((`EB...CTRL p value`<0.05 & `EB...CTRL  FDR`<0.2 & (`EB...CTRL  Ratio`>=1.5|`EB...CTRL  Ratio`<=1/1.5))
                                                   |(`EB.T.E...CTRL p value`<0.05 & `EB.T.E...CTRL  FDR`<0.2 & (`EB.T.E...CTRL  Ratio`>=1.5|`EB.T.E...CTRL  Ratio`<=1/1.5)))

  pnd90_EB <- pnd90_re %>% mutate(EB.sig=if_else(`EB...CTRL p value`<0.05&`EB...CTRL  FDR`<0.2& (`EB...CTRL  Ratio`>=1.5|`EB...CTRL  Ratio`<=1/1.5), 1,0)) %>% 
    select(protein_name,EB.sig, `EB...CTRL  Ratio`,`EB...CTRL p value`,`EB...CTRL  FDR`)
  pnd100_EB <- pnd100_re %>% mutate(EB.sig=if_else(`EB...CTRL p value`<0.05&`EB...CTRL  FDR`<0.2& (`EB...CTRL  Ratio`>=1.5|`EB...CTRL  Ratio`<=1/1.5), 1,0)) %>% 
    select(protein_name,EB.sig, `EB...CTRL  Ratio`,`EB...CTRL p value`,`EB...CTRL  FDR`)
  pnd145_EB <- pnd145_re %>% mutate(EB.sig=if_else(`EB...CTRL p value`<0.05&`EB...CTRL  FDR`<0.2& (`EB...CTRL  Ratio`>=1.5|`EB...CTRL  Ratio`<=1/1.5), 1,0)) %>% 
    select(protein_name,EB.sig, `EB...CTRL  Ratio`,`EB...CTRL p value`,`EB...CTRL  FDR`)
  pnd200_EB <- pnd200_re %>% mutate(EB.sig=if_else(`EB...CTRL p value`<0.05&`EB...CTRL  FDR`<0.2& (`EB...CTRL  Ratio`>=1.5|`EB...CTRL  Ratio`<=1/1.5), 1,0)) %>% 
    select(protein_name,EB.sig, `EB...CTRL  Ratio`,`EB...CTRL p value`,`EB...CTRL  FDR`)
  
  pnd100_EBTE <- pnd100_re %>% mutate(EB.T.E.sig=if_else(`EB.T.E...CTRL p value`<0.05&`EB.T.E...CTRL  FDR`<0.2& (`EB.T.E...CTRL  Ratio`>=1.5|`EB.T.E...CTRL  Ratio`<=1/1.5), 1,0)) %>%
    select(protein_name,EB.T.E.sig, `EB.T.E...CTRL  Ratio`,`EB.T.E...CTRL p value`,`EB.T.E...CTRL  FDR`)
  pnd145_EBTE <- pnd145_re %>% mutate(EB.T.E.sig=if_else(`EB.T.E...CTRL p value`<0.05&`EB.T.E...CTRL  FDR`<0.2& (`EB.T.E...CTRL  Ratio`>=1.5|`EB.T.E...CTRL  Ratio`<=1/1.5), 1,0)) %>% 
    select(protein_name,EB.T.E.sig, `EB.T.E...CTRL  Ratio`,`EB.T.E...CTRL p value`,`EB.T.E...CTRL  FDR`)
  pnd200_EBTE <- pnd200_re %>% mutate(EB.T.E.sig=if_else(`EB.T.E...CTRL p value`<0.05&`EB.T.E...CTRL  FDR`<0.2& (`EB.T.E...CTRL  Ratio`>=1.5|`EB.T.E...CTRL  Ratio`<=1/1.5), 1,0)) %>% 
    select(protein_name,EB.T.E.sig, `EB.T.E...CTRL  Ratio`,`EB.T.E...CTRL p value`,`EB.T.E...CTRL  FDR`)
  
  unique(c(pnd90_EB$protein_name, pnd100_EB$protein_name, pnd145_EB$protein_name,pnd200_EB$protein_name,pnd100_EBTE$protein_name,pnd145_EBTE$protein_name,pnd200_EBTE$protein_name))
  
  names(pnd90_EB)[-1] <- paste("PND90",names(pnd90_EB)[-1])
  names(pnd100_EB)[-1] <- paste("PND100",names(pnd100_EB)[-1])
  names(pnd145_EB)[-1] <- paste("PND145",names(pnd145_EB)[-1])
  names(pnd200_EB)[-1] <- paste("PND200",names(pnd200_EB)[-1])
  
  names(pnd100_EBTE)[-1] <- paste("PND100",names(pnd100_EBTE)[-1])
  names(pnd145_EBTE)[-1] <- paste("PND145",names(pnd145_EBTE)[-1])
  names(pnd200_EBTE)[-1] <- paste("PND200",names(pnd200_EBTE)[-1])
  
  pnd90_EB1 <- inner_join(proteinid,pnd90_EB)
  pnd100_EB1 <- inner_join(proteinid,pnd100_EB)
  pnd145_EB1 <- inner_join(proteinid,pnd145_EB)
  pnd200_EB1 <- inner_join(proteinid,pnd200_EB)
  
  pnd100_EBTE1 <- inner_join(proteinid,pnd100_EBTE)
  pnd145_EBTE1 <- inner_join(proteinid,pnd145_EBTE)
  pnd200_EBTE1 <- inner_join(proteinid,pnd200_EBTE)
  
  pnd100 <- inner_join(pnd100_EB,pnd100_EBTE )
  pnd145 <- inner_join(pnd145_EB,pnd145_EBTE )
  pnd200 <- inner_join(pnd200_EB,pnd200_EBTE )
  
  pnd100a <- inner_join(proteinid,pnd100)
  pnd145a <- inner_join(proteinid,pnd145)
  pnd200a <- inner_join(proteinid,pnd200)
  
  # sum(pnd90_EB$`PND90 EB.sig`,na.rm = T)
  # sum(pnd100_EB$`PND100 EB.sig`,na.rm = T)
  # sum(pnd145_EB$`PND145 EB.sig`,na.rm = T)
  # sum(pnd200_EB$`PND200 EB.sig`,na.rm = T)
  # 
  # sum(pnd100_EBTE$`PND100 EB.T.E.sig`,na.rm = T)
  # sum(pnd145_EBTE$`PND145 EB.T.E.sig`,na.rm = T)
  # sum(pnd200_EBTE$`PND200 EB.T.E.sig`,na.rm = T)
  
pndsig <- plyr::join_all(list(pnd90_EB,pnd100_EB,pnd100_EBTE,pnd145_EB,pnd145_EBTE,pnd200_EB,pnd200_EBTE), by="protein_name", type="full", match="all")
sig <- paste(pndsig$`PND90 EB.sig`,pndsig$`PND100 EB.sig`,pndsig$`PND100 EB.T.E.sig`, pndsig$`PND145 EB.sig`,pndsig$`PND145 EB.T.E.sig`,pndsig$`PND200 EB.sig`,pndsig$`PND200 EB.T.E.sig`)

pndsig <- cbind(sigstatus=sig, pndsig)
pndsig1<- inner_join(proteinid, pndsig)
write.xlsx(list(pnd=pndsig1,pnd90=pnd90_EB1,
                pnd100=pnd100a,
                pnd145=pnd145a,
                pnd200=pnd200a), "statistic result_NAminmin05_normalizedlog2_10052020_.xlsx")

  
  write.xlsx(list(pnd90=pnd90_re,pnd100=pnd100_re,pnd145=pnd145_re, pnd200=pnd200_re), "statistic result_NAminmin05_normalizedlog2.xlsx" )
  
  write.xlsx(list(pnd90=pnd90_rep005fc1d5fdr02,pnd100=pnd100_rep005fc1d5fdr02,pnd145=pnd145_rep005fc1d5fdr02, pnd200=pnd200_rep005fc1d5fdr02), "statistic result_NAminmin05_normalizedlog2 p005fc1d5fdr02.xlsx" )
  
  ### overlap between groups
  rm(list=ls())
  pnd <- read.xlsx("C:/zhijuncao/nori/proteomics/result/nest0repminmin05normlog/statistic result_NAminmin05_normalizedlog2_10052020.xlsx", sheet = "pnd")
 
  pnd1 <- pnd %>% unite(sig, c(PND90.EB.sig, PND100.EB.sig, PND100.EB.T.E.sig, PND145.EB.sig, PND145.EB.T.E.sig,PND200.EB.sig, PND200.EB.T.E.sig,),sep = "")
  write.xlsx(pnd1, "C:/zhijuncao/nori/proteomics/result/nest0repminmin05normlog/statistic result_NAminmin05_normalizedlog2_12032020.xlsx")
 
  ###significantly changed
  length(analytes)
  analytes <- unique(c(pnd90_rep005fc1d5fdr02$protein_name,pnd100_rep005fc1d5fdr02$protein_name,pnd145_rep005fc1d5fdr02$protein_name,pnd200_rep005fc1d5fdr02$protein_name))
  pdf("protein imputation normalized_log2_significantly changed points errorbar.pdf", width = 4, height = 3)
  for (analyte in analytes){
    p <- ggplot(imp_sumnormalized_log2, aes(PND, get(analyte),group=Treatment, color=Treatment))+
      theme_bw()+ 
      #facet_wrap(.~PND, nrow = 1)+
      #geom_boxplot()+
      geom_point(position=position_dodge(0.6),size=0.7,alpha=0.4)+
      scale_color_manual(values=c('black', 'red', 'purple'))+
      labs(y=analyte, x="", title = wrapper(paste(proteinid$description[proteinid$protein_name==analyte],
                                                  ':',
                                                  proteinid$gene_name[proteinid$protein_name==analyte]),40))+
      theme(axis.title=element_text(size=10),axis.text=element_text(size=10),axis.text.x=element_text(angle=90))+
      theme(plot.title=element_text(size=10))+
      #theme(plot.margin=unit(c(1,1,1,1),"line"))+
      
      stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                   geom="errorbar",  width=0.4, position=position_dodge(0.6)) +
      stat_summary(fun.y=mean, geom="point", shape=4, position=position_dodge(0.6))
    
    print(p)
  }
  dev.off()

##volcanoplot  
  names(pnd90_re)
  
  names(pnd100_re)
  
  protein <- "protein_name"
  protein <- "description"
  
  pdf("volcanoplotp05fc1d5fdr0d2 label description nest0nocontam NAminmin05_normalizedlog2_1.pdf", 6,6)
  pnd90_re <- pnd90_re %>% filter(contaminant==0, nest_under==0)
   pfc <- pnd90_re[,c(protein, "EB...CTRL p value",  "EB...CTRL  Ratio", "EB...CTRL  FDR")]
   volcanoplotfdr(pfc,title="pnd90 EB",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)
  
   pnd100_re <- pnd100_re%>% filter(contaminant==0, nest_under==0)
  pfc <- pnd100_re[,c(protein, "EB...CTRL p value",  "EB...CTRL  Ratio", "EB...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd100 EB",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)
  
  pfc <- pnd100_re[,c(protein, "EB.T.E...CTRL p value",  "EB.T.E...CTRL  Ratio", "EB.T.E...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd100 EB.T.E",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)
  
  pnd145_re <-  pnd145_re%>% filter(contaminant==0, nest_under==0)
  pfc <- pnd145_re[,c(protein, "EB...CTRL p value",  "EB...CTRL  Ratio", "EB...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd145 EB",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)
  
  pfc <- pnd145_re[,c(protein, "EB.T.E...CTRL p value",  "EB.T.E...CTRL  Ratio", "EB.T.E...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd145 EB.T.E",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)
  
  pnd200_re <- pnd200_re%>% filter(contaminant==0, nest_under==0)
  pfc <- pnd200_re[,c(protein, "EB...CTRL p value",  "EB...CTRL  Ratio", "EB...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd200 EB",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)
  
  pfc <- pnd200_re[,c(protein, "EB.T.E...CTRL p value",  "EB.T.E...CTRL  Ratio", "EB.T.E...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd200 EB.T.E",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=0.2)
  
  dev.off()
  
  
  
  pdf("volcanoplotp05fc1d5fdr1 label description.pdf", 6,6)
  
  pfc <- pnd90_re[,c(protein, "EB...CTRL p value",  "EB...CTRL  Ratio", "EB...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd90 EB",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=1)
  
  pfc <- pnd100_re[,c(protein, "EB...CTRL p value",  "EB...CTRL  Ratio", "EB...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd100 EB",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=1)
  
  pfc <- pnd100_re[,c(protein, "EB.T.E...CTRL p value",  "EB.T.E...CTRL  Ratio", "EB.T.E...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd100 EB.T.E",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=1)
  
  pfc <- pnd145_re[,c(protein, "EB...CTRL p value",  "EB...CTRL  Ratio", "EB...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd145 EB",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=1)
  
  pfc <- pnd145_re[,c(protein, "EB.T.E...CTRL p value",  "EB.T.E...CTRL  Ratio", "EB.T.E...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd145 EB.T.E",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=1)
  
  pfc <- pnd200_re[,c(protein, "EB...CTRL p value",  "EB...CTRL  Ratio", "EB...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd200 EB",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=1)
  
  pfc <- pnd200_re[,c(protein, "EB.T.E...CTRL p value",  "EB.T.E...CTRL  Ratio", "EB.T.E...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd200 EB.T.E",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff=1)
  
  dev.off()
  
  pnd90_re1 <- aov_dunn(data=pnd90_log2) 
  pnd100_re1 <- aov_dunn(data=pnd100_log2)
  pnd145_re1 <- aov_dunn(data=pnd145_log2)
  pnd200_re1 <- aov_dunn(data=pnd200_log2)
  
  fdrcutoff=0.2
  pdf("volcanoplotp05fc1d5fdr0d2 label description nest0nocontam lessthanhalfmissing noimputation.pdf", 6,6)

  pfc <- pnd90_re1[,c(protein, "EB...CTRL p value",  "EB...CTRL  Ratio", "EB...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd90 EB",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff= fdrcutoff)
  
  pfc <- pnd100_re1[,c(protein, "EB...CTRL p value",  "EB...CTRL  Ratio", "EB...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd100 EB",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff= fdrcutoff)
  
  pfc <- pnd100_re1[,c(protein, "EB.T.E...CTRL p value",  "EB.T.E...CTRL  Ratio", "EB.T.E...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd100 EB.T.E",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff= fdrcutoff)
  
  pfc <- pnd145_re1[,c(protein, "EB...CTRL p value",  "EB...CTRL  Ratio", "EB...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd145 EB",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff= fdrcutoff)
  
  pfc <- pnd145_re1[,c(protein, "EB.T.E...CTRL p value",  "EB.T.E...CTRL  Ratio", "EB.T.E...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd145 EB.T.E",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff= fdrcutoff)
  
  pfc <- pnd200_re1[,c(protein, "EB...CTRL p value",  "EB...CTRL  Ratio", "EB...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd200 EB",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff= fdrcutoff)
  
  pfc <- pnd200_re1[,c(protein, "EB.T.E...CTRL p value",  "EB.T.E...CTRL  Ratio", "EB.T.E...CTRL  FDR")]
  volcanoplotfdr(pfc,title="pnd200 EB.T.E",fccutoff=1.5, label="labelsig",labelsize=3, pcutoff=0.05, fdrcutoff= fdrcutoff)
  
  dev.off()
