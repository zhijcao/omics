library(VennDiagram)
library(rJava)

source("C:/zhijuncao/skyline/QA/qadataanalysis/zcfunction/calculate_overlap_function.R")

# calculate and venndiagram demmo
# first<-c(1:30)
# second<-c(5:35)
# third<-c(10:40)
# four<-c(15:50)
# five<-c(20:55)
# 
# x2<-list(first=first, second=second)
# x3<-list(first=first, second=second, third=third)
# x4<-list(first=first, second=second, third=third, four=four)
# x5<-list(first=first, second=second, third=third, four=four,five=five)
# x<-x4
# savename<-"fourset"
# 
# overlap2<-col(x2)
# overlap3<-col(x3)
# overlap4<-col(x4)
# overlap4<-col(x5)
# vennplot(x2,"twosetvenn")
# vennplot(x3,"threesetvenn")
# vennplot(x4,"foursetvenn")
# vennplot(x5,"fivesetvenn")





# venn plot function definition

vennplot<-function(x,savename, labelsize=1.5){

  if (length(x)==2) {
# Complex two set Venn 
venn.plot <- venn.diagram(x ,
  filename = paste(savename,".tiff",collapse=""),
  col = "black",
  lty = 0,
  #lwd = 1,
  fill = c("cornflowerblue", "darkorchid1"),
  alpha = 0.75,
  label.col = "white",
  cex = 2,
  fontfamily = "serif",
  fontface = "bold",
  #cat.col = c("cornflowerblue", "darkorchid1"),
  cat.col = c("black", "black"),
  cat.cex = labelsize,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  cat.dist = c(0.05, 0.05),
  cat.pos = c(-20, 14))
}
  else if (length(x)==3) {

# Complex three set Venn
venn.plot <- venn.diagram(x ,
  category.names = names(x),
  filename = paste(savename,".tiff",collapse=""),
  col = "black",
  lty = 0,
  output = TRUE,
  height = 3000,
  width = 3000,
  resolution = 1000,
  compression = 'lzw',
  units = 'px',
  lwd = 1,
  
  fill = c('yellow', 'purple', 'green'),
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = labelsize,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "serif",
  rotation = 1,
  print.mode="raw")
  }
  else if (length(x)==4) {

# four set venn
venn.plot <- venn.diagram(x[1:4],
  filename = paste(savename,".tiff",collapse=""),
  col = "black",
  lty = 0,
  output = TRUE,
  height = 3000,
  width = 3000,
  resolution = 1000,
  compression = 'lzw',
  lwd = 1,
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  # label.col = c("orange", "white", "darkorchid4", "white", "white", "white",
  #               "white", "white", "darkblue", "white",
  #               "white", "white", "white", "darkgreen", "white"),
  label.col = c("black", "black", "black", "black", "black", "black",
                "black", "black", "black", "black",
                "black", "black", "black", "black", "black"),
  cex = 0.8,
  fontfamily = "sans",
  fontface = "bold",
  cat.pos=c(-5,5,0,0),
  cat.dist=c(0.2,0.2,0.1,0.1),
  #cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.col = c("black", "black", "black", "black"),
  cat.cex = labelsize,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  cat.default.pos = "outer",
  print.mode="raw")
  }
  
  else if (length(x)==5) {

# five-set Venn Diagram
    # x <- pepsequence_overlap
    # savename <- "tem"
venn.plot <- venn.diagram(x,
  filename = paste(savename,".tiff",collapse=""),
  lty = 0,
  col = "black",
  #fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  fill = c(rgb(246,138,50,maxColorValue=255), rgb(194,194,194,maxColorValue=255),rgb(0,158,115,maxColorValue=255), rgb(204,121,167,maxColorValue=255), rgb(239,227,66,maxColorValue=255)),
  alpha = 0.50,
  cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
          1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5)*1.7,
  #cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.col = c("black", "black", "black", "black", "black"),
  cat.cex = labelsize,
  cat.fontface = "bold",
  cat.dist=c(0.2,0.3,0.2, 0.25, 0.3),
   margin = 0.25)
  }
  else print ("no venn drawed")
}
