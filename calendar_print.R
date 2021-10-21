cal <- function(month, year) {
  
  if(!require(chron)) stop('Unable to load chron package')
  
  if(missing(year) && missing(month)) {
    tmp <- month.day.year(Sys.Date())
    year <- tmp$year
    month <- tmp$month
  }
  
  
  if(missing(year) || missing(month)){  # year calendar
    if(missing(year)) year <- month
    par(mfrow=c(4,3))
    tmp <- seq.dates( from=julian(1,1,year), to=julian(12,31,year) )
    tmp2 <- month.day.year(tmp)
    wd <- do.call(day.of.week, tmp2)
    par(mar=c(1.5,1.5,2.5,1.5))
    for(i in 1:12){
      w <- tmp2$month == i
      cs <- cumsum(wd[w]==0)
      if(cs[1] > 0) cs <- cs - 1
      nr <- max( cs ) + 1
      plot.new()
      plot.window( xlim=c(0,6), ylim=c(0,nr+1) )
      text( wd[w], nr - cs -0.5 , tmp2$day[w] )
      title( main=month.name[i] )
      text( 0:6, nr+0.5, c('S','M','T','W','T','F','S') )
    }
    
  } else {  # month calendar
    
    ld <- seq.dates( from=julian(month,1,year), length=2, by='months')[2]-1
    days <- seq.dates( from=julian(month,1,year), to=ld)
    tmp <- month.day.year(days)
    wd <- do.call(day.of.week, tmp)
    cs <- cumsum(wd == 0)
    if(cs[1] > 0) cs <- cs - 1
    nr <- max(cs) + 1
    par(oma=c(2,2,7,2))
    par(mfrow=c(nr,7))
    par(mar=c(0,0,0,0))
    for(i in seq_len(wd[1])){ 
      plot.new()
      #box()
    }
    day.name <- c('Sun','Mon','Tues','Wed','Thur','Fri','Sat')
    drivers <- c("-", "Hailin Tang",  "Zhijun Cao", "Alokita Kamakar", "Huanli Liu", "Seongwon Nho", "-")
    
    for(i in tmp$day){
      plot.new()
      box()
      text(0,1, i, adj=c(0,1))
      if(i < 8) mtext( day.name[wd[i]+1], line=0.5,
                       at=grconvertX(0.5,to='ndc'), outer=TRUE ) 
      #text(0,1, i, adj=c(2,1))
      if(i < 8) mtext( drivers[wd[i]+1], line=2,
                       at=grconvertX(0.5,to='ndc'), outer=TRUE ) 
    }
    mtext(paste(month.name[month],year,sep=", "), line=3.5, at=0.5, cex=1.5, outer=TRUE)
    #box('inner') #optional 
  }
}

pdf("C:/Users/zcao/Desktop/van drive record/van 1 driver record sheet_Jan_Mar2020.pdf", width = 9,height =  3)
cal(1,2020)
cal(2,2020)
cal(3,2020)
dev.off()


pdf("C:/Users/zcao/Desktop/van drive record/van 1 driver record sheet_Apr_Jun2020.pdf", width = 9,height =  3)
cal(4,2020)
cal(5,2020)
cal(6,2020)
dev.off()

pdf("C:/Users/zcao/Desktop/van drive record/van 1 driver record sheet_Jul_Sep2020.pdf", width = 9,height =  3)
cal(7,2020)
cal(8,2020)
cal(9,2020)
dev.off()

pdf("C:/Users/zcao/Desktop/van drive record/van 1 driver record sheet_Oct_Dec2020.pdf", width = 9,height =  3)
cal(10,2020)
cal(11,2020)
cal(12,2020)
dev.off()

cal(1, 2019)

TeachingDemos:: cal(1,2019)


par(mfg=c(3,2))  # monday oct 10
text(.5,.5, 'Some\nText', cex=2)

par(mfg=c(2,3)) #Tues oct 4
text(1,1, 'Top Right', adj=c(1,1))

par(mfg=c(2,4)) # Wed oct 5
text(0,0, 'Bottom Left', adj=c(0,0))

par(mfg=c(6,2)) # oct 31
tmp.x <- runif(25)
tmp.y <- rnorm(25,tmp.x,.1)
par(usr=c( range(tmp.x), range(tmp.y) ) )
points(tmp.x,tmp.y)
