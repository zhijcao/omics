library(ggplot2)
#library(mixOmics)
library(magrittr)
#library(gridExtra)
library(dplyr)
library(tidyr)
library(tibble)
#library(gplots)
library(openxlsx)
library(purrr)

library(multcomp)
library(nlme)
library(lme4)
library(MASS)
library(broom)

library(ggpubr)
#library(ggpmisc)
library(cowplot)

rm(list=ls())
##############################################
####### removed BA analysis ###############
##############################################

#intba <- read.xlsx("C:/zhijuncao/sun/stool/Copy of Copy of BA individual data (002).xlsx", sheet = "Intetine BA long")
#cecum <- read.xlsx("C:/zhijuncao/sun/stool/Copy of Copy of BA individual data (002).xlsx", sheet = "Cecum BA long")
#stool <- read.xlsx("C:/zhijuncao/sun/stool/Copy of Copy of BA individual data (002).xlsx", sheet = "Stool BA long")

tem <- read.xlsx("C:/zhijuncao/sun/stool/mr1/Stool Heat map.xlsx", sheet = "Sheet1")
class(tem)
names(tem)
class(sample_meta)

lipid_report_clean <- function(tem=tem, sample_info=c(1:3), lipid_info=c(1:6)){
  
  sample_meta <- t(tem[sample_info,-lipid_info]) %>% data.frame()
  names(sample_meta) <- c("LCMS", "Samples", "LCMS1")
 
  lipids_meta <- tem[-sample_info,lipid_info]
  is_class <- is.na(lipids_meta) %>% rowSums(.)
  
  #create lipid class column
  Class <- c()
  for (i in 1:length(is_class)){
    if (is_class[i]) 
    {tem_class <- lipids_meta[[1]][i]}
    Class[i] <- tem_class
  }
  lipid_meta1 <- lipids_meta %>% add_column(uniID=paste0("X", 1:nrow(.)), Class=Class,.before=1)
  names(lipid_meta1)[-c(1:2)] <- tem[max(sample_info),lipid_info]
  
  #transpose feature to column
  data_df <- t(tem[-sample_info,-lipid_info]) %>% data.frame() %>% mutate_all(as.numeric)
  names(data_df) <- lipid_meta1$uniID
  
  #combined sample meta data with data
  data <- cbind(sample_meta, data_df)
  
  return (list(data=data, ids=lipid_meta1))
  
}

stool <- lipid_report_clean(tem=tem)

stool$data
stool$ids
as.numeric("10")
dim(feature_meta )
names(class_name)
class_name[1]
names(Class)
lipidids <- 
tem_t <- t(tem)
names(tem_t)


intest <- read.xlsx("C:/zhijuncao/sun/stool/mr1/Intestine heat map.xlsx", sheet = "clean")
cecum <- read.xlsx("C:/zhijuncao/sun/stool/mr1/Cecum heat map.xlsx", sheet = "clean")
stool <- read.xlsx("C:/zhijuncao/sun/stool/mr1/Stool Heat map.xlsx", sheet = "clean")
blood <- read.xlsx("C:/zhijuncao/sun/stool/mr1/Blood heat map.xlsx", sheet = "clean")

names(intest)[1:10]
names(cecum)[1:10]
names(stool)[1:10]
names(blood)[1:10]

###change lipid in row to in column

featuretocolunm <- function(data=intest){
  df <- t(data[, -c(1:7)]) %>% data.frame()
  samples <- names(data)[-c(1:7)]
  Type <- substr(samples, 1,2)
  Day <- substr(samples, 3,4)
  df1 <- cbind(data.frame(Sample=samples, Type=Type, Day=Day), df)
  lipidid <- cbind(uniID=names(df), data[,1:2])
  return (list(df1=df1, lipidid=lipidid))
}


intestl<- featuretocolunm(data=intest)
cecuml<- featuretocolunm(data=cecum)
stooll<- featuretocolunm(data=stool)
bloodl<- featuretocolunm(data=blood)

write.xlsx(intestl,"C:/zhijuncao/sun/stool/mr1/intestine_R_09242021.xlsx" )
write.xlsx(cecuml,"C:/zhijuncao/sun/stool/mr1/cecum_R_09242021.xlsx" )
write.xlsx(stooll,"C:/zhijuncao/sun/stool/mr1/stool_R_09242021.xlsx" )
write.xlsx(bloodl,"C:/zhijuncao/sun/stool/mr1/blood_R.xlsx" )

stoolid <- read.xlsx("C:/zhijuncao/sun/stool/mr1/stool_R.xlsx", sheet="df1" )


### remove lipid with # of missing 

rmNA <- function(data=stooll$df1, nonNA=0.5){
  sumNA <- data %>% map_df(~sum(is.na(.)))
  return (data[,(sumNA/nrow(data))<=nonNA]) 
}

### replace missing value with 0.5 minimum
repNA_05min <- function(x){
  #x is vector or list
  if (sum(is.na(x))!=length(x))
    x[is.na(x)] <- min(x,na.rm = TRUE)*0.5
  return (x)
}



bloodwithNA <- bloodl$df1 %>% rmNA
cecumwithNA <-cecuml$df1 %>% rmNA
stoolwithNA <- stooll$df1  %>% rmNA
intestwithNA <- intestl$df1 %>% rmNA

blood05min <- bloodl$df1 %>% rmNA %>% mutate_at(-c(1:3), repNA_05min)
cecum05min <-cecuml$df1 %>% rmNA %>% mutate_at(-c(1:3), ~repNA_05min(.))

#stool05min <- stooll$df1 %>% filter(X101!="nodrug")  %>% dplyr::select(-X101) %>% rmNA %>% mutate_at(-c(1:3), ~repNA_05min(as.numeric(.)))
stool05min_aminobutyrate <-stooll$df1 %>% dplyr::select(-X101) %>% rmNA %>% mutate_at(-c(1:3), ~repNA_05min(as.numeric(.)))
stool05min <-stooll$df1 %>% dplyr::select(-X101) %>% rmNA %>% mutate_at(-c(1:3), ~repNA_05min(as.numeric(.)))
names(stooll$df1)
intest05min <- intestl$df1 %>% rmNA %>% mutate_at(-c(1:3), repNA_05min)


dim(bloodl$df1)
dim(cecuml$df1)
dim(stooll$df1)
dim(intestl$df1)
#
dim(blood05min)
dim(cecum05min)
dim(stool05min)
dim(intest05min)



###correlation

### return coefficient and p value for correlation test
cor_p <- function(...){
  re <- cor.test(...)
  cor_p <- c(re$estimate, p=re$p.value)
  return (cor_p)
}

### return coefficients and p values as dataframe
cor_p_df <- function(data=data, x="X184",Y=c("X182", "X184","X183"), ids=lipidid){
  df <- map2_df(x, Y, ~cor_p(data[[.x]], data[[.y]], method="pearson")) %>% 
  add_column(uniID=Y, .before=1) %>% inner_join(ids, .)
  return (df)
}

cor_p_df(data=intest05min, x="X184",Y=c("X182", "X184","X183"),  ids=lipidid)

###all, WT, KO correlation
wt_ko_cor_p <- function (data=intest05min, x="X184",Y=c("X182", "X184","X183") ,  ids=intestl$lipidid){
  data1 <- data %>% filter(Type=="WT") 
  data2 <- data %>% filter(Type=="KO") 
  all <- cor_p_df(data=data, x=x, Y=Y, ids=ids)
  WT <- cor_p_df(data=data1, x=x, Y=Y, ids=ids)
  KO <- cor_p_df(data=data2, x=x, Y=Y, ids=ids)
  df <- full_join(WT,KO, by=c("uniID", "Class","ID"), suffix = c("_WT", "_KO")) %>% 
    full_join(all, by=c("uniID", "Class","ID"))
  return (df)
}

names(cecum05min)
names(stool05min)
names(intest05min)
ids <- intestl$lipidid
x <- ids$uniID[ids$ID=="Cefoperazone"]
Y <- names(intest05min)[-c(1:3)]



###Riboflavin and IMIDAZOLE_PROPIONATE in stool
stool_Riboflavin <- wt_ko_cor_p(data=stool05min, x=stooll$lipidid$uniID[stooll$lipidid$ID=="Riboflavin"],Y=names(stool05min)[-c(1:3)],  ids=stooll$lipidid)
stool_IMIDAZOLE_PROPIONATE <- wt_ko_cor_p(data=stool05min, x=stooll$lipidid$uniID[stooll$lipidid$ID=="IMIDAZOLE_PROPIONATE"],Y=names(stool05min)[-c(1:3)],  ids=stooll$lipidid)

###Riboflavin in cecal
cecum_Riboflavin <- wt_ko_cor_p(data=cecum05min, x=cecuml$lipidid$uniID[cecuml$lipidid$ID=="Riboflavin"],Y=names(cecum05min)[-c(1:3)],  ids=cecuml$lipidid)


###Riboflavin in intestine
intest_Riboflavin <- wt_ko_cor_p(data=intest05min, x=intestl$lipidid$uniID[intestl$lipidid$ID=="Riboflavin"],Y=names(intest05min)[-c(1:3)],  ids=intestl$lipidid)

write.xlsx(list(stool_Riboflavin=stool_Riboflavin,
                stool_IMIDAZOLE_PROPIONATE=stool_IMIDAZOLE_PROPIONATE,
                cecum_Riboflavin=cecum_Riboflavin,
                intest_Riboflavin=intest_Riboflavin), "C:/zhijuncao/sun/stool/mr1/imidazole_riboflavin_correlations.xlsx")

stool_Riboflavin_sig <- stool_Riboflavin %>% filter(p_WT<0.05|p_KO<0.05|p<0.05)
stool_IMIDAZOLE_PROPIONATE_sig <- stool_IMIDAZOLE_PROPIONATE %>% filter(p_WT<0.05|p_KO<0.05|p<0.05)
cecum_Riboflavin_sig <- cecum_Riboflavin %>% filter(p_WT<0.05|p_KO<0.05|p<0.05)
intest_Riboflavin_sig <- intest_Riboflavin %>% filter(p_WT<0.05|p_KO<0.05|p<0.05)

write.xlsx(list(stool_Riboflavin=stool_Riboflavin_sig,
                stool_IMIDAZOLE_PROPIONATE=stool_IMIDAZOLE_PROPIONATE_sig,
                cecum_Riboflavin=cecum_Riboflavin_sig,
                intest_Riboflavin=intest_Riboflavin_sig), "C:/zhijuncao/sun/stool/mr1/imidazole_riboflavin_correlations_sig.xlsx")

pdf("C:/zhijuncao/sun/stool/mr1/stool_Riboflavin_corplot_sig.pdf", width = 8, height = 4)

cor_plot(data=stool05min, xvariable="Riboflavin",Y=stool_Riboflavin_sig$uniID, ids=stooll$lipidid)

dev.off()

pdf("C:/zhijuncao/sun/stool/mr1/stool_IMIDAZOLE_PROPIONATE_corplot_sig.pdf", width = 8, height = 4)

cor_plot(data=stool05min, xvariable="IMIDAZOLE_PROPIONATE",Y=stool_IMIDAZOLE_PROPIONATE_sig$uniID, ids=stooll$lipidid)

dev.off()

###
pdf("C:/zhijuncao/sun/stool/mr1/cecum_Riboflavin_corplot_sig.pdf", width = 8, height = 4)

cor_plot(data=cecum05min, xvariable="Riboflavin",Y=cecum_Riboflavin_sig$uniID, ids=cecuml$lipidid)

dev.off()

###
pdf("C:/zhijuncao/sun/stool/mr1/intest_Riboflavin_corplot_sig.pdf", width = 8, height = 4)

cor_plot(data=intest05min, xvariable="Riboflavin",Y=intest_Riboflavin_sig$uniID, ids=intestl$lipidid)

dev.off()




##cefoperazion
intest_cefoperazone <- wt_ko_cor_p(data=intest05min, x=intestl$lipidid$uniID[intestl$lipidid$ID=="Cefoperazone"],Y=names(intest05min)[-c(1:3)],  ids=intestl$lipidid)

stool_cefoperazone <- wt_ko_cor_p(data=stool05min, x=stooll$lipidid$uniID[stooll$lipidid$ID=="Cefoperazone"],Y=names(stool05min)[-c(1:3)],  ids=stooll$lipidid)

cecum_cefoperazone <- wt_ko_cor_p(data=cecum05min, x=cecuml$lipidid$uniID[cecuml$lipidid$ID=="Cefoperazone"],Y=names(cecum05min)[-c(1:3)],  ids=cecuml$lipidid)

#aminobutyrate
intest_aminobutyrate <- wt_ko_cor_p(data=intest05min, x=intestl$lipidid$uniID[intestl$lipidid$ID=="AMINOBUTYRATE"],Y=names(intest05min)[-c(1:3)],  ids=intestl$lipidid)

stool_aminobutyrate <- wt_ko_cor_p(data=stool05min, x=stooll$lipidid$uniID[stooll$lipidid$ID=="_2_AMINOBUTYRATE"],Y=names(stool05min)[-c(1:3)],  ids=stooll$lipidid)

stool_aminobutyrate_1 <- wt_ko_cor_p(data=stool05min_aminobutyrate, x=stooll$lipidid$uniID[stooll$lipidid$ID=="_2_AMINOBUTYRATE"],Y=names(stool05min_aminobutyrate)[-c(1:3)],  ids=stooll$lipidid)


cecum_aminobutyrate <- wt_ko_cor_p(data=cecum05min, x=cecuml$lipidid$uniID[cecuml$lipidid$ID=="_2_AMINOBUTYRATE"],Y=names(cecum05min)[-c(1:3)],  ids=cecuml$lipidid)

write.xlsx(list(intest_cefoperazone=intest_cefoperazone,
                stool_cefoperazone=stool_cefoperazone,
                cecum_cefoperazone=cecum_cefoperazone,
                intest_aminobutyrate=intest_aminobutyrate,
                stool_aminobutyrate=stool_aminobutyrate,
                cecum_aminobutyrate=cecum_aminobutyrate,
                stool_aminobutyrate_1=stool_aminobutyrate_1
                ),
           "C:/zhijuncao/sun/stool/mr1/cefoperazon or aminobutyrate correlation_1.xlsx")

intest_cefoperazone_sig <- intest_cefoperazone %>% filter(p_WT<0.05|p_KO<0.05|p<0.05)
stool_cefoperazone_sig <- stool_cefoperazone %>% filter(p_WT<0.05|p_KO<0.05|p<0.05)
cecum_cefoperazone_sig <- cecum_cefoperazone %>% filter(p_WT<0.05|p_KO<0.05|p<0.05)

intest_aminobutyrate_sig <- intest_aminobutyrate %>% filter(p_WT<0.05|p_KO<0.05|p<0.05)
stool_aminobutyrate_sig <- stool_aminobutyrate %>% filter(p_WT<0.05|p_KO<0.05|p<0.05)
cecum_aminobutyrate_sig <- cecum_aminobutyrate %>% filter(p_WT<0.05|p_KO<0.05|p<0.05)

stool_aminobutyrate_1_sig <- stool_aminobutyrate_1 %>% filter(p_WT<0.05|p_KO<0.05|p<0.05)

write.xlsx(list(intest_cefoperazone=intest_cefoperazone_sig,
                stool_cefoperazone=stool_cefoperazone_sig,
                cecum_cefoperazone=cecum_cefoperazone_sig,
                intest_aminobutyrate=intest_aminobutyrate_sig,
                stool_aminobutyrate=stool_aminobutyrate_sig,
                cecum_aminobutyrate=cecum_aminobutyrate_sig,
                stool_aminobutyrate_1_sig=stool_aminobutyrate_1_sig
),
"C:/zhijuncao/sun/stool/mr1/cefoperazon or aminobutyrate correlation_1_sig.xlsx")


###correlation plot
stat_regline_equation
cor_plot <- function(data=intest05min, xvariable="Cefoperazone",Y=intest_cefoperazone_sig$uniID, ids=intestl$lipidid){
  x <- ids$uniID[ids$ID==xvariable]
  for (y in Y){
    ylab <- ids$ID[ids$uniID==y]
    cat(y, " ")
    topy <- max(data[[y]],na.rm = T)
    p1 <- ggscatter(data, x = x, y = y, add = "reg.line") +
      facet_wrap(~Type, ncol = 2, scales="fixed")+
      stat_cor(label.y=topy*0.82) +
      stat_regline_equation(label.y= topy*0.9)+
      labs(x= xvariable, y=ylab)
    
    p2 <- ggscatter(data, x = x, y = y, add = "reg.line") +
      facet_wrap(~"KO+WT")+
      stat_cor(label.y=topy*0.82) +
      stat_regline_equation(label.y= topy*0.9)+
      labs(x= xvariable, y=ylab)
    
    #grid.arrange(p2, p1,ncol=2)
    print(plot_grid(p2, p1, align = "h", nrow = 1, rel_widths = c(1.1/3,  1.9/3)))
    }
}

#intest Cefoperazone
pdf("C:/zhijuncao/sun/stool/mr1/intestine_Cefoperazone_corplot_sig.pdf", width = 8, height = 4)

cor_plot(data=intest05min, xvariable="Cefoperazone",Y=intest_cefoperazone_sig$uniID, ids=intestl$lipidid)

dev.off()

#stool Cefoperazone
pdf("C:/zhijuncao/sun/stool/mr1/stool_Cefoperazone_corplot_sig.pdf", width = 8, height = 4)

cor_plot(data=stool05min, xvariable="Cefoperazone",Y=stool_cefoperazone_sig$uniID, ids=stooll$lipidid)

dev.off()

#cecum Cefoperazone
pdf("C:/zhijuncao/sun/stool/mr1/cecum_Cefoperazone_corplot_sig.pdf", width = 8, height = 4)

cor_plot(data=cecum05min, xvariable="Cefoperazone",Y=cecum_cefoperazone_sig$uniID, ids=cecuml$lipidid)

dev.off()
#################################
#intest AMINOBUTYRATE
pdf("C:/zhijuncao/sun/stool/mr1/intestine_AMINOBUTYRATE_corplot_sig.pdf", width = 8, height = 4)

cor_plot(data=intest05min, xvariable="AMINOBUTYRATE",Y=intest_aminobutyrate_sig$uniID, ids=intestl$lipidid)

dev.off()

#stool _2_AMINOBUTYRATE
pdf("C:/zhijuncao/sun/stool/mr1/stool_2_AMINOBUTYRATE_corplot_sig.pdf", width = 8, height = 4)

cor_plot(data=stool05min, xvariable="_2_AMINOBUTYRATE",Y=stool_aminobutyrate_sig$uniID, ids=stooll$lipidid)

dev.off()

#stool _2_AMINOBUTYRATE
pdf("C:/zhijuncao/sun/stool/mr1/stool_2_AMINOBUTYRATE_corplot_1_sig.pdf", width = 8, height = 4)

cor_plot(data=stool05min_aminobutyrate, xvariable="_2_AMINOBUTYRATE",Y=stool_aminobutyrate_1_sig$uniID, ids=stooll$lipidid)

dev.off()

#cecum _2_AMINOBUTYRATE
pdf("C:/zhijuncao/sun/stool/mr1/cecum_2_AMINOBUTYRATE_corplot_sig.pdf", width = 8, height = 4)

cor_plot(data=cecum05min, xvariable="_2_AMINOBUTYRATE",Y=cecum_aminobutyrate_sig$uniID, ids=cecuml$lipidid)

dev.off()



str(lm_sum)

sqrt(0.3176)
write.xlsx(intestl,"C:/zhijuncao/sun/stool/mr1/intestine_R.xlsx" )
write.xlsx(cecuml,"C:/zhijuncao/sun/stool/mr1/cecum_R.xlsx" )
write.xlsx(stooll,"C:/zhijuncao/sun/stool/mr1/stool_R.xlsx" )
write.xlsx(bloodl,"C:/zhijuncao/sun/stool/mr1/blood_R.xlsx" )




stoolid <- read.xlsx("C:/zhijuncao/sun/stool/mr1/stool_R.xlsx", sheet="df1" )# add animal id
stool05min <- stoolid  %>% rmNA %>% mutate_at(-c(1:4), repNA_05min)

se <- function(x)sd(x)/sqrt(length(x))

data_summary <- function(data= blood05min, lipidid=bloodl$lipidid){
  
  summary <- data %>% group_by(Type,Day) %>% summarise_if(is.numeric,list(mean=mean, sd=sd)) %>% 
  pivot_longer(-c(1:2), names_to="summary", values_to="value") %>% 
  separate(summary, into = c("uniID","meansd"), sep="_") %>% 
  pivot_wider(id_cols = c(1:3),names_from = meansd, values_from = value) %>% 
  mutate(meansd=paste(round(mean,2),round(sd,2), sep = " ? "))%>% unite("Type_Day", Type, Day) %>% 
  pivot_wider(id_cols = uniID, names_from = Type_Day, values_from = meansd) %>% inner_join(lipidid, .)
  
  return (summary)
}

blood05min_summary <- data_summary(data= blood05min, lipidid=bloodl$lipidid)
cecum05min_summary <- data_summary(data= cecum05min, lipidid=cecuml$lipidid)
stool05min_summary <- data_summary(data= stool05min, lipidid=stooll$lipidid)
intest05min_summary <- data_summary(data=intest05min, lipidid=intestl$lipidid)

write.xlsx(list(blood05min=blood05min_summary,
                cecum05min=cecum05min_summary,
                stool05min=stool05min_summary,
                intest05min=intest05min_summary), "C:/zhijuncao/sun/stool/mr1/mean_sd_table NAreplace05minium.xlsx")

####Wilcox and Welch t test
noerrorttest <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(list(p.value=NA)) else return(obj)
}

wilcox_t <- function(intba1 = intba1) {
  WT_d0 <- intba1 %>% filter(Type == "WT", Day == "D0")
  WT_d1 <- intba1 %>% filter(Type == "WT", Day == "D1")
  WT_d3 <- intba1 %>% filter(Type == "WT", Day == "D3")
  WT_d5 <- intba1 %>% filter(Type == "WT", Day == "D5")
  
  KO_d0 <- intba1 %>% filter(Type == "KO", Day == "D0")
  KO_d1 <- intba1 %>% filter(Type == "KO", Day == "D1")
  KO_d3 <- intba1 %>% filter(Type == "KO", Day == "D3")
  KO_d5 <- intba1 %>% filter(Type == "KO", Day == "D5")
  
  
  ba <- names(intba1)[-c(1:3)]
  #i <- ba[15]
  
  re <- list()
  for (i in ba) {
    #Wilcoxon rank sum test
    WT_d1_d0_wil <-  wilcox.test(WT_d1[, i], WT_d0[, i])$p.value
    WT_d3_d0_wil <-  wilcox.test(WT_d3[, i], WT_d0[, i])$p.value
    WT_d5_d0_wil <-  wilcox.test(WT_d5[, i], WT_d0[, i])$p.value
    
    KO_d1_d0_wil <-  wilcox.test(KO_d1[, i], KO_d0[, i])$p.value
    KO_d3_d0_wil <-  wilcox.test(KO_d3[, i], KO_d0[, i])$p.value
    KO_d5_d0_wil <-  wilcox.test(KO_d5[, i], KO_d0[, i])$p.value
    
    d0_KO_WT_wil <-  wilcox.test(WT_d0[, i], KO_d0[, i])$p.value
    d1_KO_WT_wil <-  wilcox.test(WT_d1[, i], KO_d1[, i])$p.value
    d3_KO_WT_wil <-  wilcox.test(WT_d3[, i], KO_d3[, i])$p.value
    d5_KO_WT_wil <-  wilcox.test(WT_d5[, i], KO_d5[, i])$p.value
    
    #fold change
    WT_d1_d0_fc <- mean(WT_d1[, i]) / mean(WT_d0[, i])
    WT_d3_d0_fc <- mean(WT_d3[, i]) / mean(WT_d0[, i])
    WT_d5_d0_fc <- mean(WT_d5[, i]) / mean(WT_d0[, i])
    
    KO_d1_d0_fc <- mean(KO_d1[, i]) / mean(KO_d0[, i])
    KO_d3_d0_fc <- mean(KO_d3[, i]) / mean(KO_d0[, i])
    KO_d5_d0_fc <- mean(KO_d5[, i]) / mean(KO_d0[, i])
    
    d0_KO_WT_fc <- mean(KO_d0[, i]) / mean(WT_d0[, i])
    d1_KO_WT_fc <- mean(KO_d1[, i]) / mean(WT_d1[, i])
    d3_KO_WT_fc <- mean(KO_d3[, i]) / mean(WT_d3[, i])
    d5_KO_WT_fc <- mean(KO_d5[, i]) / mean(WT_d5[, i])
    
    #t test
    WT_d1_d0_t <-  noerrorttest(log(WT_d1[, i]), log(WT_d0[, i]))$p.value
    WT_d3_d0_t <-  noerrorttest(log(WT_d3[, i]), log(WT_d0[, i]))$p.value
    WT_d5_d0_t <-  noerrorttest(log(WT_d5[, i]), log(WT_d0[, i]))$p.value
    
    KO_d1_d0_t <-  noerrorttest(log(KO_d1[, i]), log(KO_d0[, i]))$p.value
    KO_d3_d0_t <-  noerrorttest(log(KO_d3[, i]), log(KO_d0[, i]))$p.value
    KO_d5_d0_t <-  noerrorttest(log(KO_d5[, i]), log(KO_d0[, i]))$p.value
    
    d0_KO_WT_t <-  noerrorttest(log(WT_d0[, i]), log(KO_d0[, i]))$p.value
    d1_KO_WT_t <-  noerrorttest(log(WT_d1[, i]), log(KO_d1[, i]))$p.value
    d3_KO_WT_t <-  noerrorttest(log(WT_d3[, i]), log(KO_d3[, i]))$p.value
    d5_KO_WT_t <-  noerrorttest(log(WT_d5[, i]), log(KO_d5[, i]))$p.value
    
    re[[i]] <- c(
      WT_d1_d0_fc,
      WT_d1_d0_wil,
      WT_d1_d0_t,
      WT_d3_d0_fc,
      WT_d3_d0_wil,
      WT_d3_d0_t,
      WT_d5_d0_fc,
      WT_d5_d0_wil,
      WT_d5_d0_t,
      KO_d1_d0_fc,
      KO_d1_d0_wil,
      KO_d1_d0_t,
      KO_d3_d0_fc,
      KO_d3_d0_wil,
      KO_d3_d0_t,
      KO_d5_d0_fc,
      KO_d5_d0_wil,
      KO_d5_d0_t,
      d0_KO_WT_fc,
      d0_KO_WT_wil,
      d0_KO_WT_t,
      d1_KO_WT_fc,
      d1_KO_WT_wil,
      d1_KO_WT_t,
      d3_KO_WT_fc,
      d3_KO_WT_wil,
      d3_KO_WT_t,
      d5_KO_WT_fc,
      d5_KO_WT_wil,
      d5_KO_WT_t
    )
    
  }
  re_name <- c(
    "WT_d1_d0_fc",
    "WT_d1_d0_wil",
    "WT_d1_d0_t",
    "WT_d3_d0_fc",
    "WT_d3_d0_wil",
    "WT_d3_d0_t",
    "WT_d5_d0_fc",
    "WT_d5_d0_wil",
    "WT_d5_d0_t",
    "KO_d1_d0_fc",
    "KO_d1_d0_wil",
    "KO_d1_d0_t",
    "KO_d3_d0_fc",
    "KO_d3_d0_wil",
    "KO_d3_d0_t",
    "KO_d5_d0_fc",
    "KO_d5_d0_wil",
    "KO_d5_d0_t",
    "d0_KO_WT_fc",
    "d0_KO_WT_wil",
    "d0_KO_WT_t",
    "d1_KO_WT_fc",
    "d1_KO_WT_wil",
    "d1_KO_WT_t",
    "d3_KO_WT_fc",
    "d3_KO_WT_wil",
    "d3_KO_WT_t",
    "d5_KO_WT_fc",
    "d5_KO_WT_wil",
    "d5_KO_WT_t"
  )
  re_df <- do.call(rbind, re) %>% data.frame()
  names(re_df) <- re_name
  
  re_df$WT_d1_d0_wil_FDR <- p.adjust(re_df$WT_d1_d0_wil, method = "BH")
  re_df$WT_d3_d0_wil_FDR <- p.adjust(re_df$WT_d3_d0_wil, method = "BH")
  re_df$WT_d5_d0_wil_FDR <- p.adjust(re_df$WT_d5_d0_wil, method = "BH")
  
  re_df$KO_d1_d0_wil_FDR <- p.adjust(re_df$KO_d1_d0_wil, method = "BH")
  re_df$KO_d3_d0_wil_FDR <- p.adjust(re_df$KO_d3_d0_wil, method = "BH")
  re_df$KO_d5_d0_wil_FDR <- p.adjust(re_df$KO_d5_d0_wil, method = "BH")
  
  re_df$d0_KO_WT_wil_FDR <- p.adjust(re_df$d0_KO_WT_wil, method = "BH")
  re_df$d1_KO_WT_wil_FDR <- p.adjust(re_df$d1_KO_WT_wil, method = "BH")
  re_df$d3_KO_WT_wil_FDR <- p.adjust(re_df$d3_KO_WT_wil, method = "BH")
  re_df$d5_KO_WT_wil_FDR <- p.adjust(re_df$d5_KO_WT_wil, method = "BH")
  
  re_df$WT_d1_d0_t_FDR <- p.adjust(re_df$WT_d1_d0_t, method = "BH")
  re_df$WT_d3_d0_t_FDR <- p.adjust(re_df$WT_d3_d0_t, method = "BH")
  re_df$WT_d5_d0_t_FDR <- p.adjust(re_df$WT_d5_d0_t, method = "BH")
  
  re_df$KO_d1_d0_t_FDR <- p.adjust(re_df$KO_d1_d0_t, method = "BH")
  re_df$KO_d3_d0_t_FDR <- p.adjust(re_df$KO_d3_d0_t, method = "BH")
  re_df$KO_d5_d0_t_FDR <- p.adjust(re_df$KO_d5_d0_t, method = "BH")
  
  re_df$d0_KO_WT_t_FDR <- p.adjust(re_df$d0_KO_WT_t, method = "BH")
  re_df$d1_KO_WT_t_FDR <- p.adjust(re_df$d1_KO_WT_t, method = "BH")
  re_df$d3_KO_WT_t_FDR <- p.adjust(re_df$d3_KO_WT_t, method = "BH")
  re_df$d5_KO_WT_t_FDR <- p.adjust(re_df$d5_KO_WT_t, method = "BH")
  re_df <- re_df %>% tibble::rownames_to_column("uniID")
  return(re_df)
  
}

wilcox_t_stool <- function(intba1 = intba1) {
  WT_d0 <- intba1 %>% filter(Type == "WT", Day == "D0")
  WT_d1 <- intba1 %>% filter(Type == "WT", Day == "D1")
  WT_d2 <- intba1 %>% filter(Type == "WT", Day == "D2")
  WT_d3 <- intba1 %>% filter(Type == "WT", Day == "D3")
  WT_d4 <- intba1 %>% filter(Type == "WT", Day == "D4")
  WT_d5 <- intba1 %>% filter(Type == "WT", Day == "D5")
  
  KO_d0 <- intba1 %>% filter(Type == "KO", Day == "D0")
  KO_d1 <- intba1 %>% filter(Type == "KO", Day == "D1")
  KO_d2 <- intba1 %>% filter(Type == "KO", Day == "D2")
  KO_d3 <- intba1 %>% filter(Type == "KO", Day == "D3")
  KO_d4 <- intba1 %>% filter(Type == "KO", Day == "D4")
  KO_d5 <- intba1 %>% filter(Type == "KO", Day == "D5")
  
  
  ba <- names(intba1)[-c(1:3)]
  #i <- ba[15]
  
  re <- list()
  for (i in ba) {
    #Wilcoxon rank sum test
    WT_d1_d0_wil <-  wilcox.test(WT_d1[, i], WT_d0[, i])$p.value
    WT_d2_d0_wil <-  wilcox.test(WT_d2[, i], WT_d0[, i])$p.value
    WT_d3_d0_wil <-  wilcox.test(WT_d3[, i], WT_d0[, i])$p.value
    WT_d4_d0_wil <-  wilcox.test(WT_d4[, i], WT_d0[, i])$p.value
    WT_d5_d0_wil <-  wilcox.test(WT_d5[, i], WT_d0[, i])$p.value
    
    KO_d1_d0_wil <-  wilcox.test(KO_d1[, i], KO_d0[, i])$p.value
    KO_d2_d0_wil <-  wilcox.test(KO_d2[, i], KO_d0[, i])$p.value
    KO_d3_d0_wil <-  wilcox.test(KO_d3[, i], KO_d0[, i])$p.value
    KO_d4_d0_wil <-  wilcox.test(KO_d4[, i], KO_d0[, i])$p.value
    KO_d5_d0_wil <-  wilcox.test(KO_d5[, i], KO_d0[, i])$p.value
    
    d0_KO_WT_wil <-  wilcox.test(WT_d0[, i], KO_d0[, i])$p.value
    d1_KO_WT_wil <-  wilcox.test(WT_d1[, i], KO_d1[, i])$p.value
    d2_KO_WT_wil <-  wilcox.test(WT_d2[, i], KO_d2[, i])$p.value
    d3_KO_WT_wil <-  wilcox.test(WT_d3[, i], KO_d3[, i])$p.value
    d4_KO_WT_wil <-  wilcox.test(WT_d4[, i], KO_d4[, i])$p.value
    d5_KO_WT_wil <-  wilcox.test(WT_d5[, i], KO_d5[, i])$p.value
    
    #fold change
    WT_d1_d0_fc <- mean(WT_d1[, i]) / mean(WT_d0[, i])
    WT_d2_d0_fc <- mean(WT_d2[, i]) / mean(WT_d0[, i])
    WT_d3_d0_fc <- mean(WT_d3[, i]) / mean(WT_d0[, i])
    WT_d4_d0_fc <- mean(WT_d4[, i]) / mean(WT_d0[, i])
    WT_d5_d0_fc <- mean(WT_d5[, i]) / mean(WT_d0[, i])
    
    KO_d1_d0_fc <- mean(KO_d1[, i]) / mean(KO_d0[, i])
    KO_d2_d0_fc <- mean(KO_d2[, i]) / mean(KO_d0[, i])
    KO_d3_d0_fc <- mean(KO_d3[, i]) / mean(KO_d0[, i])
    KO_d4_d0_fc <- mean(KO_d4[, i]) / mean(KO_d0[, i])
    KO_d5_d0_fc <- mean(KO_d5[, i]) / mean(KO_d0[, i])
    
    d0_KO_WT_fc <- mean(KO_d0[, i]) / mean(WT_d0[, i])
    d1_KO_WT_fc <- mean(KO_d1[, i]) / mean(WT_d1[, i])
    d2_KO_WT_fc <- mean(KO_d2[, i]) / mean(WT_d2[, i])
    d3_KO_WT_fc <- mean(KO_d3[, i]) / mean(WT_d3[, i])
    d4_KO_WT_fc <- mean(KO_d4[, i]) / mean(WT_d4[, i])
    d5_KO_WT_fc <- mean(KO_d5[, i]) / mean(WT_d5[, i])
    
    #t test
    WT_d1_d0_t <-  noerrorttest(log(WT_d1[, i]), log(WT_d0[, i]))$p.value
    WT_d2_d0_t <-  noerrorttest(log(WT_d2[, i]), log(WT_d0[, i]))$p.value
    WT_d3_d0_t <-  noerrorttest(log(WT_d3[, i]), log(WT_d0[, i]))$p.value
    WT_d4_d0_t <-  noerrorttest(log(WT_d4[, i]), log(WT_d0[, i]))$p.value
    WT_d5_d0_t <-  noerrorttest(log(WT_d5[, i]), log(WT_d0[, i]))$p.value
    
    KO_d1_d0_t <-  noerrorttest(log(KO_d1[, i]), log(KO_d0[, i]))$p.value
    KO_d2_d0_t <-  noerrorttest(log(KO_d2[, i]), log(KO_d0[, i]))$p.value
    KO_d3_d0_t <-  noerrorttest(log(KO_d3[, i]), log(KO_d0[, i]))$p.value
    KO_d4_d0_t <-  noerrorttest(log(KO_d4[, i]), log(KO_d0[, i]))$p.value
    KO_d5_d0_t <-  noerrorttest(log(KO_d5[, i]), log(KO_d0[, i]))$p.value
    
    d0_KO_WT_t <-  noerrorttest(log(WT_d0[, i]), log(KO_d0[, i]))$p.value
    d1_KO_WT_t <-  noerrorttest(log(WT_d1[, i]), log(KO_d1[, i]))$p.value
    d2_KO_WT_t <-  noerrorttest(log(WT_d2[, i]), log(KO_d2[, i]))$p.value
    d3_KO_WT_t <-  noerrorttest(log(WT_d3[, i]), log(KO_d3[, i]))$p.value
    d4_KO_WT_t <-  noerrorttest(log(WT_d4[, i]), log(KO_d4[, i]))$p.value
    d5_KO_WT_t <-  noerrorttest(log(WT_d5[, i]), log(KO_d5[, i]))$p.value
    
    re[[i]] <- c(
      WT_d1_d0_fc,
      WT_d1_d0_wil,
      WT_d1_d0_t,
      WT_d2_d0_fc,
      WT_d2_d0_wil,
      WT_d2_d0_t,
      WT_d3_d0_fc,
      WT_d3_d0_wil,
      WT_d3_d0_t,
      WT_d4_d0_fc,
      WT_d4_d0_wil,
      WT_d4_d0_t,
      WT_d5_d0_fc,
      WT_d5_d0_wil,
      WT_d5_d0_t,
      KO_d1_d0_fc,
      KO_d1_d0_wil,
      KO_d1_d0_t,
      KO_d2_d0_fc,
      KO_d2_d0_wil,
      KO_d2_d0_t,
      KO_d3_d0_fc,
      KO_d3_d0_wil,
      KO_d3_d0_t,
      KO_d4_d0_fc,
      KO_d4_d0_wil,
      KO_d4_d0_t,
      KO_d5_d0_fc,
      KO_d5_d0_wil,
      KO_d5_d0_t,
      d0_KO_WT_fc,
      d0_KO_WT_wil,
      d0_KO_WT_t,
      d1_KO_WT_fc,
      d1_KO_WT_wil,
      d1_KO_WT_t,
      d2_KO_WT_fc,
      d2_KO_WT_wil,
      d2_KO_WT_t,
      d3_KO_WT_fc,
      d3_KO_WT_wil,
      d3_KO_WT_t,
      d4_KO_WT_fc,
      d4_KO_WT_wil,
      d4_KO_WT_t,
      d5_KO_WT_fc,
      d5_KO_WT_wil,
      d5_KO_WT_t
    )
    
  }
  re_name <- c(
    "WT_d1_d0_fc",
    "WT_d1_d0_wil",
    "WT_d1_d0_t",
    "WT_d2_d0_fc",
    "WT_d2_d0_wil",
    "WT_d2_d0_t",
    "WT_d3_d0_fc",
    "WT_d3_d0_wil",
    "WT_d3_d0_t",
    "WT_d4_d0_fc",
    "WT_d4_d0_wil",
    "WT_d4_d0_t",
    "WT_d5_d0_fc",
    "WT_d5_d0_wil",
    "WT_d5_d0_t",
    "KO_d1_d0_fc",
    "KO_d1_d0_wil",
    "KO_d1_d0_t",
    "KO_d2_d0_fc",
    "KO_d2_d0_wil",
    "KO_d2_d0_t",
    "KO_d3_d0_fc",
    "KO_d3_d0_wil",
    "KO_d3_d0_t",
    "KO_d4_d0_fc",
    "KO_d4_d0_wil",
    "KO_d4_d0_t",
    "KO_d5_d0_fc",
    "KO_d5_d0_wil",
    "KO_d5_d0_t",
    "d0_KO_WT_fc",
    "d0_KO_WT_wil",
    "d0_KO_WT_t",
    "d1_KO_WT_fc",
    "d1_KO_WT_wil",
    "d1_KO_WT_t",
    "d2_KO_WT_fc",
    "d2_KO_WT_wil",
    "d2_KO_WT_t",
    "d3_KO_WT_fc",
    "d3_KO_WT_wil",
    "d3_KO_WT_t",
    "d4_KO_WT_fc",
    "d4_KO_WT_wil",
    "d4_KO_WT_t",
    "d5_KO_WT_fc",
    "d5_KO_WT_wil",
    "d5_KO_WT_t"
  )
  re_df <- do.call(rbind, re) %>% data.frame()
  names(re_df) <- re_name
  
  re_df$WT_d1_d0_wil_FDR <- p.adjust(re_df$WT_d1_d0_wil, method = "BH")
  re_df$WT_d2_d0_wil_FDR <- p.adjust(re_df$WT_d2_d0_wil, method = "BH")
  re_df$WT_d3_d0_wil_FDR <- p.adjust(re_df$WT_d3_d0_wil, method = "BH")
  re_df$WT_d4_d0_wil_FDR <- p.adjust(re_df$WT_d4_d0_wil, method = "BH")
  re_df$WT_d5_d0_wil_FDR <- p.adjust(re_df$WT_d5_d0_wil, method = "BH")
  
  re_df$KO_d1_d0_wil_FDR <- p.adjust(re_df$KO_d1_d0_wil, method = "BH")
  re_df$KO_d2_d0_wil_FDR <- p.adjust(re_df$KO_d2_d0_wil, method = "BH")
  re_df$KO_d3_d0_wil_FDR <- p.adjust(re_df$KO_d3_d0_wil, method = "BH")
  re_df$KO_d4_d0_wil_FDR <- p.adjust(re_df$KO_d4_d0_wil, method = "BH")
  re_df$KO_d5_d0_wil_FDR <- p.adjust(re_df$KO_d5_d0_wil, method = "BH")
  
  re_df$d0_KO_WT_wil_FDR <- p.adjust(re_df$d0_KO_WT_wil, method = "BH")
  re_df$d1_KO_WT_wil_FDR <- p.adjust(re_df$d1_KO_WT_wil, method = "BH")
  re_df$d2_KO_WT_wil_FDR <- p.adjust(re_df$d2_KO_WT_wil, method = "BH")
  re_df$d3_KO_WT_wil_FDR <- p.adjust(re_df$d3_KO_WT_wil, method = "BH")
  re_df$d4_KO_WT_wil_FDR <- p.adjust(re_df$d4_KO_WT_wil, method = "BH")
  re_df$d5_KO_WT_wil_FDR <- p.adjust(re_df$d5_KO_WT_wil, method = "BH")
  
  re_df$WT_d1_d0_t_FDR <- p.adjust(re_df$WT_d1_d0_t, method = "BH")
  re_df$WT_d2_d0_t_FDR <- p.adjust(re_df$WT_d2_d0_t, method = "BH")
  re_df$WT_d3_d0_t_FDR <- p.adjust(re_df$WT_d3_d0_t, method = "BH")
  re_df$WT_d4_d0_t_FDR <- p.adjust(re_df$WT_d4_d0_t, method = "BH")
  re_df$WT_d5_d0_t_FDR <- p.adjust(re_df$WT_d5_d0_t, method = "BH")
  
  re_df$KO_d1_d0_t_FDR <- p.adjust(re_df$KO_d1_d0_t, method = "BH")
  re_df$KO_d2_d0_t_FDR <- p.adjust(re_df$KO_d2_d0_t, method = "BH")
  re_df$KO_d3_d0_t_FDR <- p.adjust(re_df$KO_d3_d0_t, method = "BH")
  re_df$KO_d4_d0_t_FDR <- p.adjust(re_df$KO_d4_d0_t, method = "BH")
  re_df$KO_d5_d0_t_FDR <- p.adjust(re_df$KO_d5_d0_t, method = "BH")
  
  re_df$d0_KO_WT_t_FDR <- p.adjust(re_df$d0_KO_WT_t, method = "BH")
  re_df$d1_KO_WT_t_FDR <- p.adjust(re_df$d1_KO_WT_t, method = "BH")
  re_df$d2_KO_WT_t_FDR <- p.adjust(re_df$d2_KO_WT_t, method = "BH")
  re_df$d3_KO_WT_t_FDR <- p.adjust(re_df$d3_KO_WT_t, method = "BH")
  re_df$d4_KO_WT_t_FDR <- p.adjust(re_df$d4_KO_WT_t, method = "BH")
  re_df$d5_KO_WT_t_FDR <- p.adjust(re_df$d5_KO_WT_t, method = "BH")
  re_df <- re_df %>% tibble::rownames_to_column("uniID")
  return(re_df)
  
}


fc <- c("KO_d1_d0_fc", "KO_d3_d0_fc", "KO_d5_d0_fc", 
        "WT_d1_d0_fc","WT_d3_d0_fc","WT_d5_d0_fc", 
        "d0_KO_WT_fc","d1_KO_WT_fc", "d3_KO_WT_fc" , "d5_KO_WT_fc")     

blood_re_df<- wilcox_t(intba1=blood05min)
blood_re_df <- inner_join(bloodl$lipidid,blood_re_df)
blood_fc <- blood_re_df[,fc]  
names(blood_fc) <- paste(contrasts1,"Ratio", sep = ": ")
blood_fc[,"uniID"] <- blood_re_df$uniID


cecum_re_df<- wilcox_t(intba1=cecum05min)
cecum_re_df <- inner_join(cecuml$lipidid,cecum_re_df)
cecum_fc <- cecum_re_df[,fc]  
names(cecum_fc) <- paste(contrasts1,"Ratio", sep = ": ")
cecum_fc[,"uniID"] <- cecum_re_df$uniID

intest_re_df<- wilcox_t(intba1=intest05min)
intest_re_df <- inner_join(intestl$lipidid,intest_re_df)
intest_fc <- intest_re_df[,fc]  
names(intest_fc) <- paste(contrasts1,"Ratio", sep = ": ")
intest_fc[,"uniID"] <- intest_re_df$uniID


st_fc <- c("KO_d1_d0_fc", "KO_d2_d0_fc", "KO_d3_d0_fc","KO_d4_d0_fc","KO_d5_d0_fc", 
           "WT_d1_d0_fc", "WT_d2_d0_fc", "WT_d3_d0_fc","WT_d4_d0_fc","WT_d5_d0_fc", 
           "d0_KO_WT_fc", "d1_KO_WT_fc", "d2_KO_WT_fc" , "d3_KO_WT_fc", "d4_KO_WT_fc", "d5_KO_WT_fc") 

stool_re_df<- wilcox_t_stool(intba1=stool05min)
stool_re_df <- inner_join(stooll$lipidid,stool_re_df)
stool_fc <- stool_re_df[,st_fc]  
names(stool_fc) <- paste(contrasts,"Ratio", sep = ": ")
stool_fc[,"uniID"] <- stool_re_df$uniID

write.xlsx(list(blood=blood_re_df,intestine=intest_re_df, cecum=cecum_re_df, stool=stool_re_df),
           "C:/zhijuncao/sun/stool/mr1/wilcox_t_ratiopfdr.xlsx")



############linear regression###############
library(lmerTest)
library(psycho)
library(mutoss)
library(modelbased)

#devtools::install_version("psycho", version = "0.4.91", repos = "http://cran.us.r-project.org")

contrasts1 <- c("KO D0 - KO D1",  "KO D0 - KO D3", "KO D0 - KO D5", 
                "WT D0 - WT D1",  "WT D0 - WT D3", "WT D0 - WT D5",
                "WT D0 - KO D0", "WT D1 - KO D1", "WT D3 - KO D3", "WT D5 - KO D5")


contrasts <- c("KO D0 - KO D1", "KO D0 - KO D2", "KO D0 - KO D3", "KO D0 - KO D4", "KO D0 - KO D5", 
               "WT D0 - WT D1", "WT D0 - WT D2", "WT D0 - WT D3", "WT D0 - WT D4", "WT D0 - WT D5",
               "WT D0 - KO D0", "WT D1 - KO D1", "WT D2 - KO D2","WT D3 - KO D3", "WT D4 - KO D4",  "WT D5 - KO D5")
allcontrast <- paste(results$Level1, results$Level2, sep = " - ")

blood05minlog<- blood05min %>% mutate_if(is.numeric, log) 
cecum05minlog<- cecum05min %>% mutate_if(is.numeric, log)
stool05minlog<-  stool05min %>% mutate_if(is.numeric, log)
intest05minlog<- intest05min %>% mutate_if(is.numeric, log)

blood05minlog$Type <- factor(blood05minlog$Type, levels = c("WT", "KO"))
cecum05minlog$Type <- factor(cecum05minlog$Type, levels = c("WT", "KO"))
stool05minlog$Type <- factor(stool05minlog$Type, levels = c("WT", "KO"))
intest05minlog$Type <- factor(intest05minlog$Type, levels = c("WT", "KO"))


x <- "X2"
fdr <- function(x)p.adjust(x, method="BH")

#function of linear regression
ba_lm <- function(data=blood05minlog,lipidid=bloodl$lipidid, contrasts1=contrasts1){

re <- list()
ratio <- list()
for (x in names(data)[-c(1:3)]){
  cat(x,"  ")
  fml <- paste(x,"~ Type*Day")
  fml <-  as.formula(fml)
  fit <- lm(fml, data=data)
  
  results <- estimate_contrasts(fit, adjust = "none") %>% data.frame()
  row.names(results) <- paste(results$Level1, results$Level2, sep = " - ")
  results1 <- results[contrasts1, ]
  re[[x]]  <- p.adjust(results1$p, method="BH")
  ratio[[x]] <- 1/exp(results1$Difference)
}

re1 <- do.call(rbind, re) %>% data.frame()
ratio1 <- do.call(rbind, ratio) %>% data.frame()
names(re1) <-contrasts1
names(ratio1) <- paste(contrasts1, "ratio", sep = ": ")
re_fdr <- re1 %>% mutate_all(fdr)

names(re_fdr) <- paste(names(re_fdr), "FDR", sep = ": ")
names(re1) <- paste(names(re1), "p", sep = ": ")
pfdr <- cbind(re1,re_fdr, ratio1) %>% tibble::rownames_to_column("uniID")

pfdr1 <- inner_join(lipidid, pfdr)

return(pfdr1)
}
ba_lmmix <- function(data=stool05minlog,lipidid=stooll$lipidid, contrasts1=contrasts1){
  
  re <- list()
  ratio <- list()
  for (x in names(data)[-c(1:4)]){
    cat(x,"  ")
    fml <- paste(x,"~ Type*Day + (1|AnimalID)")
    fml <-  as.formula(fml)
    fit <- glmer(fml, data=data)
    
    results <- estimate_contrasts(fit, adjust = "none") %>% data.frame()
    row.names(results) <- paste(results$Level1, results$Level2, sep = " - ")
    results1 <- results[contrasts1, ]
    re[[x]]  <- p.adjust(results1$p, method="BH")
    ratio[[x]] <- 1/exp(results1$Difference)
  }
  
  re1 <- do.call(rbind, re) %>% data.frame()
  ratio1 <- do.call(rbind, ratio) %>% data.frame()
  names(re1) <-contrasts1
  names(ratio1) <- paste(contrasts1, "ratio", sep = ": ")
  re_fdr <- re1 %>% mutate_all(fdr)
  
  names(re_fdr) <- paste(names(re_fdr), "FDR", sep = ": ")
  names(re1) <- paste(names(re1), "p", sep = ": ")
  pfdr <- cbind(re1,re_fdr, ratio1) %>% tibble::rownames_to_column("uniID")
  
  pfdr1 <- inner_join(lipidid, pfdr)
  
  return(pfdr1)
}

blood_re <- ba_lm(data=blood05minlog,lipidid=bloodl$lipidid, contrasts1=contrasts1)
cecum_re <- ba_lm(data=cecum05minlog,lipidid=cecuml$lipidid, contrasts1=contrasts1)
stool_re <- ba_lm(data=stool05minlog,lipidid=stooll$lipidid, contrasts1=contrasts)
intest_re <- ba_lm(data=intest05minlog,lipidid=intestl$lipidid, contrasts1=contrasts1)

stool_re_mixmode <- ba_lmmix(data=stool05minlog,lipidid=stooll$lipidid, contrasts1=contrasts)
write.xlsx(stool_re_mixmode,
           "C:/zhijuncao/sun/stool/mr1/stool_re_mixmode_ratio_by_exp.xlsx")
write.xlsx(list(bloodpfdr=blood_re, cecumpfdr=cecum_re, stoolpfdr=stool_re, intestpfdr=intest_re),
           "C:/zhijuncao/sun/stool/mr1/lm_p_fdr_ratio_by_exp.xlsx")



library(ggthemes)
library(gridExtra)


################################################
####   Manuscript figures
#################################################

####manuscript figure 1
lipidid <- bloodl$lipidid
analyte <- "X1"

point_errorbar_plot <- function(pdfname="c:/zhijuncao/sun/stool/mr1/test.pdf",data=blood05min,analytes=names(blood05min)[c(5:7)], lipidid=bloodl$lipidid){
  pdf(pdfname,width=3.5,height=2.6)
  for (analyte in analytes){
  title <- lipidid$ID[lipidid$uniID==analyte]
  p <- ggplot(data, aes(Day, get(analyte)))+
    facet_wrap(~Type)+
    theme_classic()+ 
    geom_point(alpha=0.5)+
    stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                 geom="errorbar",  width=0.4, position=position_dodge(0.6)) +
    stat_summary(fun.y=mean, geom="point", shape=4, position=position_dodge(0.6))+
    labs(y="Intensity", x="Day", title = title)+
    theme(axis.title=element_text(size=10),axis.text=element_text(size=10),axis.text.x=element_text(angle=0))+
    theme(plot.title=element_text(size=7))
  
  print(p)
}
dev.off()
}
point_errorbar_plot(pdfname="c:/zhijuncao/sun/stool/mr1/point_errorbar_plot_blood.pdf",data=blood05min,analytes=names(blood05min)[-c(1:3)], lipidid=bloodl$lipidid) 
point_errorbar_plot(pdfname="c:/zhijuncao/sun/stool/mr1/point_errorbar_plot_cecum.pdf",data=cecum05min,analytes=names(cecum05min)[-c(1:3)], lipidid=cecuml$lipidid) 
point_errorbar_plot(pdfname="c:/zhijuncao/sun/stool/mr1/point_errorbar_plot_stool.pdf",data=stool05min1,analytes=names(stool05min1)[-c(1:3)], lipidid=stooll$lipidid) 
point_errorbar_plot(pdfname="c:/zhijuncao/sun/stool/mr1/point_errorbar_plot_intest.pdf",data=intest05min,analytes=names(intest05min)[-c(1:3)], lipidid=intestl$lipidid) 



###heat map
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
colors <- brewer.pal(9, "Set1")
names(colors) <- unique(bloodl$lipidid$Class)
toMean <- function(x) x/mean(x, na.rm=TRUE)

data <- bloodwithNA
lipidid <- bloodl$lipidid
prepheatmap <- function(data=bloodwithNA, lipidid=bloodl$lipidid){
  
  labrow <- left_join(data.frame(uniID=names(data)[-c(1:3)]),lipidid)
  toMean1 <- data[data$Type=="KO" ,] %>% mutate_at(-c(1:3),list(toMean))
  toMean2 <- data[data$Type=="WT" ,] %>% mutate_at(-c(1:3),list(toMean))
  toMean12 <- rbind(toMean1,toMean2)
  
  return (list(data=toMean12, labrow=labrow))
}

data <- cecumwithNA
data <- stoolwithNA
data <- intestwithNA
blood_hp <- prepheatmap(data=bloodwithNA, lipidid=bloodl$lipidid)
cecum_hp <- prepheatmap(data=cecumwithNA, lipidid=cecuml$lipidid)
stool_hp <- prepheatmap(data=stoolwithNA, lipidid=stooll$lipidid)
intest_hp <- prepheatmap(data=intestwithNA, lipidid=intestl$lipidid)

length(names(toMean12))
unique(bloodl$lipidid$Class)
unique(cecuml$lipidid$Class)
unique(stooll$lipidid$Class)
unique(intestl$lipidid$Class)

###no separated normalization
length(stoolwithNA)
colors <- brewer.pal(9, "Set1")
colColors <- c("grey", "antiquewhite4")
names(colColors) <- c("KO", "WT")
labrow <- left_join(data.frame(uniID=names(stoolwithNA)[-c(1:3, 84)]),stooll$lipidid)
names(colors) <- unique(stooll$lipidid$Class)[1:9]
tem <- stoolwithNA[,-c(1:3,84)] %>% mutate_all(as.numeric) %>% data.frame()
pdf("C:/zhijuncao/sun/stool/mr1/stoolwithNA_heatmap_color.pdf", width = 10,height = 10)
heatmap.2(as.matrix(log(t(tem))),col=bluered, Rowv = F, Colv = FALSE, margins = c(4, 15),
          labCol=stoolwithNA$Sample,labRow = labrow$ID ,trace="none",
          dendrogram="none",
          scale = c("row"),
          lhei = c(1,6),
          #lwid=c(1,7),
          main = "stool",
          RowSideColors=colors[labrow$Class],
          ColSideColors=colColors[stoolwithNA$Type],
          na.color = "grey",srtCol=90, cexCol = 0.8, cexRow = 0.5, key.xlab = "log(Intensity)",
          density.info="none",breaks = seq(-1.5, 1.5, 0.1))
legend("left", legend = names(colors),fill = colors,title = "Class", cex = 0.7)
dev.off()


labrow <- left_join(data.frame(uniID=names(cecumwithNA)[-c(1:3, 132)]),cecuml$lipidid)
names(colors) <- unique(cecuml$lipidid$Class)[1:9]
pdf("C:/zhijuncao/sun/stool/mr1/cecumwithNA_heatmap_color.pdf", width = 10,height = 10)
heatmap.2(as.matrix(log(t(cecumwithNA[,-c(1:3,132)]))),col=bluered, Rowv = F, Colv = FALSE, margins = c(4, 15),
          labCol=cecumwithNA$Sample,labRow = labrow$ID ,trace="none",
          dendrogram="none",
          scale = c("row"),
          lhei = c(1,6),
          #lwid=c(1,7),
          main = "Cecum",
          RowSideColors=colors[labrow$Class], 
          ColSideColors=colColors[cecumwithNA$Type],
          na.color = "grey",srtCol=90, cexCol = 0.8, cexRow = 0.5, key.xlab = "log(Intensity)",
          density.info="none",breaks = seq(-1.5, 1.5, 0.1))
legend("left", legend = names(colors),fill = colors,title = "Class", cex = 0.7)
dev.off()

colors <- brewer.pal(9, "Set1")[1:8]
labrow <- left_join(data.frame(uniID=names(intestwithNA)[-c(1:3, 175)]),intestl$lipidid)
names(colors) <- unique(intestl$lipidid$Class)[1:8]
pdf("C:/zhijuncao/sun/stool/mr1/intestwithNA_heatmap_color.pdf", width = 10,height = 10)
heatmap.2(as.matrix(log(t(intestwithNA[,-c(1:3,175)]))),col=bluered, Rowv = F, Colv = FALSE, margins = c(4, 15),
          labCol=intestwithNA$Sample,labRow = labrow$ID ,trace="none",
          dendrogram="none",
          scale = c("row"),
          lhei = c(1,6),
          #lwid=c(1,7),
          main = "Intestine",
          RowSideColors=colors[labrow$Class], 
          ColSideColors=colColors[intestwithNA$Type],
          na.color = "grey",srtCol=90, cexCol = 0.8, cexRow = 0.5, key.xlab = "log(Intensity)",
          density.info="none",breaks = seq(-1.5, 1.5, 0.1))
legend("left", legend = names(colors),fill = colors,title = "Class", cex = 0.7)
dev.off()



str(cecumwithNA[,-c(1:3,132)])
str(intestwithNA[,-c(1:150)])

length(cecumwithNA)
length(intestwithNA)
cecuml$lipidid
intestl$lipidid
#####normazied
toMean12 <- blood_hp$data
labrow <- blood_hp$labrow
names(colors) <- unique(bloodl$lipidid$Class)
pdf("C:/zhijuncao/sun/stool/mr1/blood_KO(tomean)_WT(tomean)_heatmap_color1.pdf", width = 10,height = 10)
heatmap.2(as.matrix(t(log(toMean12[,-c(1:3)]))),col=bluered, Rowv = F, Colv = FALSE, margins = c(4, 15),
          labCol=toMean12$Sample,labRow = labrow$ID ,trace="none",
          dendrogram="none",
          lhei = c(1,6),
          #lwid=c(1,7),
          main = "Blood",
          RowSideColors=colors[labrow$Class], 
          na.color = "grey",srtCol=90, cexCol = 0.8, cexRow = 0.5, key.xlab = "log(toMean)",
          density.info="none",breaks = seq(-1, 1, 0.1))
legend("left", legend = names(colors),fill = colors,title = "Class", cex = 0.7)
dev.off()



toMean12 <- cecum_hp$data
labrow <- cecum_hp$labrow
names(colors) <- unique(cecuml$lipidid$Class)
pdf("C:/zhijuncao/sun/stool/mr1/cecum_KO(tomean)_WT(tomean)_heatmap_color1.pdf", width = 10,height = 10)
heatmap.2(as.matrix(t(log(toMean12[,-c(1:3)]))),col=bluered, Rowv = F, Colv = FALSE, margins = c(4, 15),
          labCol=toMean12$Sample,labRow = labrow$ID ,trace="none",
          dendrogram="none",
          lhei = c(1,6),
          #lwid=c(1,7),
          main = "Cecum",
          RowSideColors=colors[labrow$Class], 
          na.color = "grey",srtCol=90, cexCol = 0.8, cexRow = 0.5, key.xlab = "log(toMean)",
          density.info="none",breaks = seq(-1, 1, 0.1))
legend("left", legend = names(colors),fill = colors,title = "Class", cex = 0.55)
dev.off()

toMean12 <- stool_hp$data
labrow <- stool_hp$labrow
names(colors) <- unique(stooll$lipidid$Class)
pdf("C:/zhijuncao/sun/stool/mr1/stool_KO(tomean)_WT(tomean)_heatmap_color1.pdf", width = 12,height = 10)
heatmap.2(as.matrix(t(log(toMean12[,-c(1:3)]))),col=bluered, Rowv = F, Colv = FALSE, margins = c(4, 15),
          labCol=toMean12$Sample,labRow = labrow$ID ,trace="none",
          dendrogram="none",
          lhei = c(1,6),
          #lwid=c(1,7),
          main = "Stool",
          RowSideColors=colors[labrow$Class], 
          na.color = "grey",srtCol=90, cexCol = 0.8, cexRow = 0.5, key.xlab = "log(toMean)",
          density.info="none",breaks = seq(-1, 1, 0.1))
legend("left", legend = names(colors),fill = colors,title = "Class", cex = 0.7)
dev.off()

toMean12 <- intest_hp$data
labrow <- intest_hp$labrow
names(colors) <- unique(intestl$lipidid$Class)
pdf("C:/zhijuncao/sun/stool/mr1/intestine_KO(tomean)_WT(tomean)_heatmap_color1.pdf", width = 10,height = 10)
heatmap.2(as.matrix(t(log(toMean12[,-c(1:3)]))),col=bluered, Rowv = F, Colv = FALSE, margins = c(4, 15),
          labCol=toMean12$Sample,labRow = labrow$ID ,trace="none",
          dendrogram="none",
          lhei = c(1,6),
          #lwid=c(1,7),
          main = "Intestine",
          RowSideColors=colors[labrow$Class], 
          na.color = "grey",srtCol=90, cexCol = 0.8, cexRow = 0.5, key.xlab = "log(toMean)",
          density.info="none",breaks = seq(-1, 1, 0.1))
legend("left", legend = names(colors)[-9],fill = colors,title = "Class", cex = 0.7)
dev.off()

####volcanoplot

bloodpfdr <- read.xlsx("C:/zhijuncao/sun/stool/mr1/lm_stat_p_fdr_ratio_by_mean.xlsx", sheet="bloodpfdr")
cecumpfdr <- read.xlsx("C:/zhijuncao/sun/stool/mr1/lm_stat_p_fdr_ratio_by_mean.xlsx", sheet="cecumpfdr")
intestpfdr <- read.xlsx("C:/zhijuncao/sun/stool/mr1/lm_stat_p_fdr_ratio_by_mean.xlsx", sheet="intestpfdr")
stoolpfdr <- read.xlsx("C:/zhijuncao/sun/stool/mr1/lm_stat_p_fdr_ratio_by_mean.xlsx", sheet="stoolpfdr_mixmode")

stoolpfdr_sig <- stoolpfdr %>% dplyr::filter(`KO.D0.-.KO.D1:.p`<0.05|
                                               `KO.D0.-.KO.D3:.p`<0.05|
                                               `KO.D0.-.KO.D5:.p`<0.05|
                                               `WT.D0.-.WT.D1:.p`<0.05|
                                               `WT.D0.-.WT.D3:.p`<0.05|
                                               `WT.D0.-.WT.D5:.p`<0.05|
                                               `WT.D0.-.KO.D0:.p`<0.05|
                                               `WT.D1.-.KO.D1:.p`<0.05|
                                               `WT.D3.-.KO.D3:.p`<0.05|
                                               `WT.D5.-.KO.D5:.p`<0.05|
                                               `KO.D0.-.KO.D2:.p`<0.05|
                                               `KO.D0.-.KO.D4:.p`<0.05|
                                               `WT.D0.-.WT.D2:.p`<0.05|
                                               `WT.D0.-.WT.D4:.p`<0.05)
  
intestpfdr_sig <- intestpfdr %>% dplyr::filter(`KO.D0.-.KO.D1:.p`<0.05|
                                        `KO.D0.-.KO.D3:.p`<0.05|
                                        `KO.D0.-.KO.D5:.p`<0.05|
                                        `WT.D0.-.WT.D1:.p`<0.05|
                                        `WT.D0.-.WT.D3:.p`<0.05|
                                        `WT.D0.-.WT.D5:.p`<0.05|
                                        `WT.D0.-.KO.D0:.p`<0.05|
                                        `WT.D1.-.KO.D1:.p`<0.05|
                                        `WT.D3.-.KO.D3:.p`<0.05|
                                        `WT.D5.-.KO.D5:.p`<0.05)

cecumpfdr_sig <- cecumpfdr %>% dplyr::filter(`KO.D0.-.KO.D1:.p`<0.05|
                                                 `KO.D0.-.KO.D3:.p`<0.05|
                                                 `KO.D0.-.KO.D5:.p`<0.05|
                                                 `WT.D0.-.WT.D1:.p`<0.05|
                                                 `WT.D0.-.WT.D3:.p`<0.05|
                                                 `WT.D0.-.WT.D5:.p`<0.05|
                                                 `WT.D0.-.KO.D0:.p`<0.05|
                                                 `WT.D1.-.KO.D1:.p`<0.05|
                                                 `WT.D3.-.KO.D3:.p`<0.05|
                                                 `WT.D5.-.KO.D5:.p`<0.05)

write.xlsx(list(stoolpfdr_sig=stoolpfdr_sig, 
                intestpfdr_sig=intestpfdr_sig,
                cecumpfdr_sig=cecumpfdr_sig), 
           "C:/zhijuncao/sun/stool/mr1/lm_stat_p05_fdr_ratio_by_mean.xlsx")

####
#### get ids for pathway analysis
####
stool_compoundid <- read.xlsx("C:/zhijuncao/sun/stool/mr1/pathway/stool_cecum_intest_p05_id.xlsx", sheet = "stool_p05_id")
cecum_compoundid <- read.xlsx("C:/zhijuncao/sun/stool/mr1/pathway/stool_cecum_intest_p05_id.xlsx", sheet = "cecum_p05_id")
intest_compoundid <- read.xlsx("C:/zhijuncao/sun/stool/mr1/pathway/stool_cecum_intest_p05_id.xlsx", sheet = "intest_p05_id")

names(cecumpfdr)[c(4:13)]
getsig_ids <- function(data=cecumpfdr, compared=names(cecumpfdr)[c(4:13)], compoundid=cecum_compoundid){
  re <- list()
  
  for (x in compared){
     cat(x)
     matchedids <- data[data[[x]]<0.05,c(x, "ID")] %>% left_join(., compoundid) 
     re[[substr(x,1,13)]] <- matchedids
  }
    return (re)
}

str(cecum_ids)
cecum_ids <- getsig_ids(data=cecumpfdr, compared=names(cecumpfdr)[c(4:13)], compoundid=cecum_compoundid)
write.xlsx(cecum_ids, "C:/zhijuncao/sun/stool/mr1/cecum_ids_forpathway.xlsx") 

intest_ids <- getsig_ids(data=intestpfdr, compared=names(intestpfdr)[c(4:13)], compoundid=intest_compoundid)
write.xlsx(intest_ids, "C:/zhijuncao/sun/stool/mr1/intest_ids_forpathway.xlsx") 

stool_ids <- getsig_ids(data=stoolpfdr, compared=names(stoolpfdr)[c(4:19)], compoundid=stool_compoundid)
write.xlsx(stool_ids, "C:/zhijuncao/sun/stool/mr1/stool_ids_forpathway.xlsx") 

library(purrr)
library(openxlsx)
library(dplyr)


###pathway results
options(repos = getOption("repos")["CRAN"])
devtools::install_github("b2slab/FELLA",  force = TRUE)
library(FELLA)

bg <- read.xlsx("C:/zhijuncao/sun/stool/mr1/pathway/background_pathway.xlsx")



### build local database

g <- buildGraphFromKEGGREST(
  organism = "hsa", 
  filter.path = "hsa01100"
)
buildDataFromGraph(g)

g.mmu <- buildGraphFromKEGGREST(
  organism = "mmu", 
  filter.path = "mmu01100"
)
buildDataFromGraph(g.mmu)

stool_ids[[names(stool_ids)[1]]]$KEGG

### load mouse kegg database
database <- loadKEGGdata(databaseDir = tail(listInternalDatabases(), 1),
             internalDir = TRUE, loadMatrix = NULL)

###pathway enrichment analysis
###stool
stool_ids[["KO.D0.-.KO.D2"]]
x <- "KO.D0.-.KO.D2"
re <- list()
for (x in names(stool_ids)){
  obj.wrap <- enrich(
  compounds = stool_ids[[x]]$KEGG, 
  methods = FELLA::listMethods()[1], 
  #compoundsBackground=bg$KEGG,
  data = database)
  
  re[[x]] <-generateResultsTableName(data = database,
                                     object = obj.wrap,
                                     method="hypergeom",
                                     threshold = 1)
}

write.xlsx(re, "C:/zhijuncao/sun/stool/mr1/pathway/stool_pathway_result1.xlsx")

###cecum
re <- list()
for (x in names(cecum_ids)){
  cat(x, " ")
  if (is.na(cecum_ids[[x]]$KEGG)) next
  obj.wrap <- enrich(
    compounds = cecum_ids[[x]]$KEGG, 
    methods = FELLA::listMethods()[1], 
    #compoundsBackground=bg$KEGG,
    data = database)
  
  re[[x]] <-generateResultsTableName(data = data,
                                     object = obj.wrap,
                                     method="hypergeom",
                                     threshold = 1)
}

write.xlsx(re, "C:/zhijuncao/sun/stool/mr1/pathway/cecum_pathway_result1.xlsx")

###intest
re <- list()
x <- "KO.D0.-.KO.D3" 
for (x in names(intest_ids)){
  obj<-try(obj.wrap <- enrich(
    compounds = intest_ids[[x]]$KEGG, 
    methods = FELLA::listMethods()[1], 
    #compoundsBackground=bg$KEGG,
    data = database), silent=TRUE)
  if (is(obj, "try-error")) next
  
  
  re[[x]] <-generateResultsTableName(data = data,
                                     object = obj.wrap,
                                     method="hypergeom",
                                     threshold = 1)
}

write.xlsx(re, "C:/zhijuncao/sun/stool/mr1/pathway/intest_pathway_result1.xlsx")

intest_ids[[x]]
nrow(intest_ids[[x]])>0
length(intest_ids[["KO.D0.-.KO.D1"]]$KEGG)

#############add hits metabolites names
library(tibble)
library(dplr)
library(magrittr)
library(FELLA)

id_to_name <- function(list.hit) {
  ids_names <- getName(data, list.hit)
  ids_names <- map2(ids_names, names(ids_names), ~ paste(.y, paste(.x, collapse = ","), sep = ": ")) %>% unlist(.)
  return (list(Compound = paste(ids_names, collapse = " | ")))
}


generateResultsTableName <- function(data = database,
                         object = obj.wrap,
                         method="hypergeom",
                         threshold = 0.1) {
  
  tab.res <- generateResultsTable(
    method = method,
    threshold = threshold, 
    object = object, 
    data = data)
  
  if (is.null(tab.res)) return (NA)
  
  mat.hypergeom <- FELLA:::getMatrix(data, "hypergeom")
  list.cpd <- getInput(object)
  
  # all the metabolites (path) and those within the list (hits)
  list.path <- apply(mat.hypergeom, 2, function(x) names(which(x)))
  
  target.path <- purrr::map(tab.res$KEGG.id, ~ list.path[[.x]])
  names(target.path) <- tab.res$KEGG.id
  
  list.hits <- lapply(target.path, function(x) intersect(x, list.cpd))
 
  ids_names <- map(list.hits, ~ id_to_name(.x))

  df <- do.call(rbind, ids_names) %>% data.frame() %>% 
    rownames_to_column("KEGG.id") %>% left_join(tab.res, .)
  return (df)
}

df1 <- generateResultsTableName(data = database,
                         object = re[["KO.D0.-.KO.D2"]],
                         method="hypergeom",
                         threshold = 0.1)


data("input.sample")
input.full <- c(input.sample, paste0("intruder", 1:10))
show(input.full)


browseVignettes("FELLA")
nchar("cecum_WT.D0.-.KO.D0")
substr("cecum_WT.D0.-.KO.D0_pathway_results.csv", 1, 19)
####

stoolpfdr_fixed <- read.xlsx("C:/zhijuncao/sun/stool/mr1/lm_stat_p_fdr_ratio_by_mean.xlsx", sheet="stoolpfdr")

names(cecumpfdr)

tem <- cecumpfdr
pdf("C:/zhijuncao/sun/stool/mr1/cecum_volcanoplot.pdf", width = 7, height=6)

tem <- bloodpfdr
pdf("C:/zhijuncao/sun/stool/mr1/blood_volcanoplot.pdf", width = 7, height=6)

tem <- intestpfdr
pdf("C:/zhijuncao/sun/stool/mr1/intestine_volcanoplot.pdf", width = 7, height=6)

tem <- stoolpfdr
pdf("C:/zhijuncao/sun/stool/mr1/stool_volcanoplot.pdf", width = 7, height=6)

tem <- stoolpfdr_fixed
pdf("C:/zhijuncao/sun/stool/mr1/stool_fixed_volcanoplot.pdf", width = 7, height=6)

pfc <- tem[, c("ID", "KO.D0.-.KO.D1:.p", "KO.D0.-.KO.D1:.Ratio", "KO.D0.-.KO.D1:.FDR")]
volcanoplotfdr(pfc,title="KO.D0.-.KO.D1",fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2)

pfc <- tem[, c("ID", "KO.D0.-.KO.D3:.p", "KO.D0.-.KO.D3:.Ratio", "KO.D0.-.KO.D3:.FDR")]
volcanoplotfdr(pfc,title="KO.D0.-.KO.D3",fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2)

pfc <- tem[, c("ID", "KO.D0.-.KO.D5:.p", "KO.D0.-.KO.D5:.Ratio", "KO.D0.-.KO.D5:.FDR")]
volcanoplotfdr(pfc,title="KO.D0.-.KO.D5",fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2)
####
pfc <- tem[, c("ID", "WT.D0.-.WT.D1:.p", "WT.D0.-.WT.D1:.Ratio", "WT.D0.-.WT.D1:.FDR")]
volcanoplotfdr(pfc,title="WT.D0.-.WT.D1",fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2)

pfc <- tem[, c("ID", "WT.D0.-.WT.D3:.p", "WT.D0.-.WT.D3:.Ratio", "WT.D0.-.WT.D3:.FDR")]
volcanoplotfdr(pfc,title="WT.D0.-.WT.D3",fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2)

pfc <- tem[, c("ID", "WT.D0.-.WT.D5:.p", "WT.D0.-.WT.D5:.Ratio", "WT.D0.-.WT.D5:.FDR")]
volcanoplotfdr(pfc,title="WT.D0.-.WT.D5",fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2)

####
pfc <- tem[, c("ID", "WT.D0.-.KO.D0:.p", "WT.D0.-.KO.D0:.Ratio", "WT.D0.-.KO.D0:.FDR")]
volcanoplotfdr(pfc,title="WT.D0.-.KO.D0",fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2)

pfc <- tem[, c("ID", "WT.D1.-.KO.D1:.p", "WT.D1.-.KO.D1:.Ratio", "WT.D1.-.KO.D1:.FDR")]
volcanoplotfdr(pfc,title="WT.D1.-.KO.D1",fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2)

pfc <- tem[, c("ID", "WT.D3.-.KO.D3:.p", "WT.D3.-.KO.D3:.Ratio", "WT.D3.-.KO.D3:.FDR")]
volcanoplotfdr(pfc,title="WT.D3.-.KO.D3",fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2)

pfc <- tem[, c("ID", "WT.D5.-.KO.D5:.p", "WT.D5.-.KO.D5:.Ratio", "WT.D5.-.KO.D5:.FDR")]
volcanoplotfdr(pfc,title="WT.D5.-.KO.D5",fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2)

#############################
###stool add day 2 and 4
pfc <- tem[, c("ID", "KO.D0.-.KO.D2:.p", "KO.D0.-.KO.D2:.Ratio", "KO.D0.-.KO.D2:.FDR")]
volcanoplotfdr(pfc,title="KO.D0.-.KO.D2",fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2)

pfc <- tem[, c("ID", "KO.D0.-.KO.D4:.p", "KO.D0.-.KO.D4:.Ratio", "KO.D0.-.KO.D4:.FDR")]
volcanoplotfdr(pfc,title="KO.D0.-.KO.D4",fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2)
###
pfc <- tem[, c("ID", "WT.D0.-.WT.D2:.p", "WT.D0.-.WT.D2:.Ratio", "WT.D0.-.WT.D2:.FDR")]
volcanoplotfdr(pfc,title="WT.D0.-.WT.D2",fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2)

pfc <- tem[, c("ID", "WT.D0.-.WT.D4:.p", "WT.D0.-.WT.D4:.Ratio", "WT.D0.-.WT.D4:.FDR")]
volcanoplotfdr(pfc,title="WT.D0.-.WT.D4",fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2)
###
pfc <- tem[, c("ID", "WT.D2.-.KO.D2:.p", "WT.D2.-.KO.D2:.Ratio", "WT.D2.-.KO.D2:.FDR")]
volcanoplotfdr(pfc,title="WT.D2.-.KO.D2",fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2)

pfc <- tem[, c("ID", "WT.D4.-.KO.D4:.p", "WT.D4.-.KO.D4:.Ratio", "WT.D4.-.KO.D4:.FDR")]
volcanoplotfdr(pfc,title="WT.D4.-.KO.D4",fccutoff=1.2, label="labelsig",labelsize=3,pointsize=1, pcutoff=0.05, fdrcutoff=2)

dev.off()
