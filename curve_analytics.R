library(minpack.lm)
library(nlmrt)
library(optimx)
library(ggplot2)
library(nlstools)
library(nlme)

source('generic.R')

# [64-bit] C:\Program Files\R\R-3.3.1

rhs1 <- function(b, mydata) {
  sum((mydata$y - ((b[1]*mydata$x*b[2]) / (b[1]*mydata$x+b[2])) - b[3])^2)
}


rhs2 <- function(b, mydata) {
  #eunsc = ydat ~ ((tdat * Pgmax)/(tdat + lo)) - RD
  sum((mydata$y - (((mydata$x*b[2]) / (mydata$x+b[1])) - b[3]))^2)
}

rhsfit1 <- function(b, mydata) {
  #ydat ~ P25*exp((ha/(8.314*298.15)) - ha/(8.314*tdat))
  sum((mydata$y - (b[1]*exp((b[2]/(8.314*298.15)) - b[2] /(8.314*mydata$x)) ) )^2)
}

rhsfit2 <- function(b, mydata) {
  #(P25*exp( (log(1 + exp((S*298.15-hd)/(8.314*298.15))) + (ha/(8.314*298.15))) 
   #         - (ha/(8.314*tdat)))) / (1 + exp((S*tdat-hd)/(8.312*tdat)))
  sum((mydata$y - (  (b[1]*exp(( log(1+ exp((b[3]*298.15-b[4])/(8.314*298.15))) 
                                 + (b[2]/(8.314*298.15)) 
                               - (b[2]/(8.314*mydata$x)))) / 
                        (1 + exp((b[3]*mydata$x-b[4])/(8.314*mydata$x)))
                      )) )^2)
}

optima <- function(start1, w, fc)
{
  ## Optim
  w.optx <- optimx(par=start1, fn=fc, mydata=w, 
                   control=list(all.methods=TRUE, save.failures=TRUE, 
                                maxit=10000))
  return(w.optx)
}
  
plotCurve <- function(x,y, fm, eq, yl, spc, tempr, repl)
{
  plot(x, y, ylab = yl,
       xlab = 'PHOTOSYNTHETIC PHOTON FLUX DENSITY', cex.axis=0.7, cex.lab=0.6)
  lines(x, fm, col = 2, lwd = 2)
  text(max(x)/2, max(y)/2, paste('specie: ', spc,
    '\ntemp:', tempr, '\n',  'rep:', repl, '\n', 
    'eq:', eq, sep=''), cex = 0.9)
}


searchOutlier <- function(data)
{
  library('outliers')
  data <- data[order(x),]
  row.names(data) <- 1:nrow(data)
  u <- chisq.out.test(data$y, variance=var(data$y), opposite = FALSE)
  i <- which(data$y == as.numeric(strsplit(u$alternative, ' ')[[1]][3]))
  data <- data[-i,]
  row.names(data) <- 1:nrow(data)
  s <- c()
  mnn <- round(data$y[1],3)
  mno<-mnn
  for(i in 2:nrow(data))
  {
    if(mno > mnn)
    {
      mno <- mno
    }
    else
    {
      mno <- mnn
    }

    if(mno > round(data$y[i],3) )
    {
  
      s <- append(s, i)
      mnn <- round(data$y[i-1],3)
    }
    else
    {
      mnn <- round(data$y[i],3)
    }
    
  }
  
  data <- data[-s,]
  row.names(data) <- 1:nrow(data)
  return(data)
}


eqFitcl1 <- function(tdat, ydat, subdata)
{
  o <- list()
  
  # Get starting values
  #cc <- coef(lm(log(ydat)~tdat,data=subdata)) 
  #cc <- c(cc)
  #names(cc) <- c("P25", "ha") 
  
  o$start<- c(P25 = 0, ha=38000)
  o$eunsc <- ydat ~ P25*exp((ha/(8.314*298.15)) - ha/(8.314*tdat))
  o$optx1 <- optima(o$start, subdata, rhsfit1)
  o1 <- nlsLM(o$eunsc, 
              start=list(P25 = o$optx1$P25[1], ha = o$optx1$ha[1]), trace = TRUE)
  
  return(o1)
}

eqFitcl2 <- function(tdat, ydat, weeddata)
{
  o <- list()
  o$start<- c(P25 = 10, ha=20000, S=800, hd=258000)
  o$eunsc <- ydat ~ (P25*exp( (log(1 + exp((S*298.15-hd)/(8.314*298.15))) + 
                                (ha/(8.314*298.15))) 
                             - (ha/(8.314*tdat)))) / 
    (1 + exp((S*tdat-hd)/(8.312*tdat)))
  
  #(b[1]*exp(( log(1+ exp((b[3]*298.15-b[4])/(8.314*298.15))) 
  #+ (b[2]/(8.314*298.15)) 
  #- (b[2]/(8.314*mydata$x)))) / 
  #(1 + exp((b[3]*mydata$x-b[4])/(8.314*mydata$x))))
  
  o$optx1 <- optima(o$start, weeddata, rhsfit2)
  o2 <- nlsLM(o$eunsc, 
               start=list(P25 = o$optx$P25[1], ha = o$optx$ha[1], S = o$optx$S[1], 
                          hd = o$optx$hd[1]), 
               trace = TRUE)
  return(o2)
}

# Remove outliers
getlogistic <- function(x, y)
{
  
  subdata <- data.frame(y=y, x=x)
  f <- as.formula(paste('y', '~', 'SSlogis(x, phi1, phi2, phi3)', sep=''))
  
  # Get starting values
  cc <- coef(lm(log(y)~x,data=subdata)) 
  cc <- c(cc, 1)
  names(cc) <- c("phi1", "phi2", "phi3") 
  #nls(log(y)~b1+(b2/b3)*(x^b3-1),data = subdata, start = cc,trace = TRUE) 
  
  gr <- nlsLM(f, data = subdata, start = cc)
  r <- nlsResiduals(gr)
  i <- union(which(r$resi2[1:nrow(subdata),2] > 1), which(r$resi2[1:nrow(subdata),2] < -1))
  # plot(weeddata[['x']], weeddata[['y']], xlab = "DAS", ylab = 'ta')
  # lines(weeddata[['x']], as.vector(predict(gr)), col="red",lty=2,lwd=3)
  # text(170, min(log(sub[[trait]]))+1.5, paste(ename, '\n', 'phi1: ', round(alpha[1],2), 
  #                                             '\n  phi2:', round(alpha[2],2), 
  #'\n',  ' phi3:', round(alpha[3],2), sep=''),cex = 0.8)
  # 
  return(i)
}


#' Fit ArrheniusFit Curve1 and Curve 2 for a given parameter
#' @param xp is the time/space relates parameter (temp)
#' @param yp the responde variable
#' @param param param for which curve if fitted
eqFit <- function(tdat, ydat, param)
{
  r <- list()
  ## Remove odd value BUT it results in ood values
  i <- getlogistic(tdat, ydat)
 
  subdata <- data.frame(y=ydat, x=tdat)
  if(param == 'lo' |  param == 'RD' ){
    eqt <- 1
    r$oo <- eqFitcl1(tdat, ydat, subdata)
  }
  else {
    eqt <- 2
    r$oo <- eqFitcl2(tdat, ydat, subdata)
  
  }
  print( r$oo)
  k <- min(ydat) + (min(ydat)/2)
  plot(tdat, ydat, main = paste('Arrhenius equations\' curves, param: ', param, 
                                sep=''),
       cex=0.8, cex.axis=0.8, cex.lab=0.6)
  lines(tdat, fitted(r$oo), col = 'red', lwd = 2, cex=0.8, cex.axis=0.8)
  #lines(tdat, fitted(r$o2), col = 'blue', lwd = 2, cex=0.8, cex.axis=0.8)
  legend(290,k, legend = c(paste('eq:', eqt, ' cor=',round(cor(tdat, fitted(r$oo)), 3), 
                                 sep='')), 
                           #paste('eq2, cor=',round(cor(tdat, fitted(r$o2)), 3), 
                            #     sep='')), 
         col=c('red'), lty=1, cex=0.8)
  r$ydat <- ydat
  return(r)
}

st <- function(x,y, opt)
{
  if(opt == '1')
  {
    mx <- max(y);
    mn<-min(y)
    mm<-0
    start1<- c(lo = mm, Pgmax = mx, RD=mn)
  }
  if(opt == '2')
  {  
    mx <- max(y) 
    mn<-min(y)
    mm<-5
    start1<- c(lo = mm, Pgmax = mx, RD=mn)
  }
  return(start1)
}

makeFunctions <- function(ydat, tdat)
{
  
  elist <- list()
  elist[['1']] <- c(eq = ydat ~ ((lo*tdat*Pgmax) / (lo*tdat+Pgmax)) - RD, 
                   f=rhs1, qy = "(lo*(Pgmax^2)) / ((lo*tdat+Pgmax)^2)")
  
  elist[['2']] <- c(eq = ydat ~ ((tdat * Pgmax)/(tdat + lo)) - RD, 
                   f=rhs2, qy = "(lo*Pgmax) / ((tdat+lo)^2)")
  elist[['10']] <- c(eq = ydat ~ P25*exp((c-ha)/(8.314*tdat)), 
                   f=rhsfit1, qy = "")
  return(elist)
}

# Calculate growth parameters
nlsLMFunction <- function(x,y, tname, spc, tempr, repl, fn)
{
  b <- c()
  tdat <- x
  ydat <- y
  elist <- makeFunctions(ydat, tdat)
  weeddata <- data.frame(y=ydat, x=tdat)
  
  for(n in fn)
  {
    optx <- optima(st(tdat, ydat, n), weeddata, elist[[n]]$f)
    
    o <- nlsLM(elist[[n]]$eq, 
                start=list(lo = optx$lo[1], Pgmax = optx$Pgmax[1], 
                           RD=optx$RD[1]), trace = F)
    lo <- optx$lo[1]
    Pgmax <- optx$Pgmax[1]
    RD<-optx$RD[1]
    b <- rbind(b, c(eq = n, rep = repl, lo = lo, Pgmax=Pgmax, RD=RD, 
                   specie = spc, temperature = tempr))
    plotCurve(tdat, ydat, fitted(o), n, 'NET PHOTOSYNTHESIS RATE',spc, tempr, 
              repl)
    plotCurve(tdat, eval(parse(text=elist[[n]]$qy)), 
              eval(parse(text=elist[[n]]$qy)), n, 
              'QUANTUM YIELD', spc, tempr, repl)
  
  }
  return(b)
}

# calculate temperature
getPar <- function(data, varn)
{
  data <- data.frame(data)
  data <- data.frame(data, tk = sapply(data[[varn]], function(x) x+273.15))
  return(data)
}

convert2factor <- function(data, variable_list)
{
  
  copydata <- data
  for(v in variable_list)
  {
    copydata[[v]]  <- factor(copydata[[v]])
  }
  return(copydata)
}

#Equations to fit curve to temp
eqFitc <- function(type, vl)
{
  v1 <- 8.314
  v2 <- 298.15
  if(type == 2)
  {
    return(log(1 + exp((vl[['S']]*v2-vl[['hd']])/(v1*v2))) +
           (vl[['ha']]/(v1*v2)))
  }
  if(type == 1)
  {
    return(vl[['ha']]/(v1*v2))
  }
  else
    print('unknown parameter')
}

eqEstimate <- function(t, vl, tp)
{
  print
  v1 <- 8.314
  v2 <- 273.15
  if(tp == 'Pgmax' | tp == 'RD')
  {
    return(exp((vl[['c']] - (vl[['ha']]) / (v1*(t + v2))) / 
           (1 + exp((vl[['S']]*(t + v2) - vl[['hd']]) / (v1*(t+v2)))))) 
  }
  if(tp == 'lo')
  {
    return(exp((vl[['c']]-(vl[['ha']])/(v1*(t+v2)))))
  }
  else
    print('unknown variable')
  
}

#' ArrheniusFit Calculations
#' @param sdata curve parameters to calculate ArrheniusFit
ArrheniusFit <- function(sdata)
{
  # Pgmax
  x <- sdata[['tk']]
  y <- sdata[['Pgmax']]
  m <- eqFit(x,y, 'Pgmax')
  m1 <- data.frame(Estimate = coef(m$oo))
  print(paste('Finished Pgmax'))
  rm(m)
  
  # lo
  x <- sdata[['tk']]
  y <- sdata$lo
  m <- eqFit(x[-5],y[-5], 'lo')
  m2 <- data.frame(Estimate = coef(m$oo))
  for(p in c('P25', 'ha', 'S', 'hd')){
    if(is.na(match(p, rownames(m2)))) {  
      m2 <- rbind(m2, K=NA)
      rownames(m2)[which(rownames(m2) == 'K')] = p
    }
  }
  print(paste('Finished lo'))
  rm(m)
  
  # RD
  x <- sdata[['tk']]
  y <- sdata$RD
  m <- eqFit(x,y, 'RD')
  m3 <- data.frame(Estimate = coef(m$oo))
  for(p in c('P25', 'ha', 'S', 'hd')){
    if(is.na(match(p, rownames(m3)))) {  
      m3 <- rbind(m3, K=NA)
      rownames(m3)[which(rownames(m3) == 'K')] = p
    }
  }
  print(paste('Finished RD'))
  
  
  mf <- cbind(Pgmax=m1$Estimate, lo=m2$Estimate, RD=m3$Estimate)
  rownames(mf) <- rownames(m1)
  mf <- rbind(mf, c=c(eqFitc(2, mf[,'Pgmax']),eqFitc(1, mf[,'lo']), 
                     eqFitc(2, mf[,'RD'])))
  return(mf)
}

# Process data from curve fitting
processCurveParam <- function(b)
{
  
  b <- convert2numeric(data.frame(b), c("eq","rep","lo","Pgmax","RD",
                                        "temperature"))
  b <- getPar(b, 'temperature')
  b <- data.frame(b)
  b <- b[order(b$tk),]
  b$RD <- -1*(b$RD)
  return(b)
}


calculateEstimates <- function(m, b)
{
  b <- data.frame(b, PgmaxAdj = sapply(b$temperature, 
                                       function(x) eqEstimate(x, m[,'Pgmax'], 
                                                              'Pgmax')))
  b <- data.frame(b, loAdj = sapply(b$temperature, 
                                    function(x) eqEstimate(x, m[,'lo'], 'lo')))
  b <- data.frame(b, RDAdj = sapply(b$temperature, 
                                    function(x) eqEstimate(x, m[,'RD'], 'RD')))
  return(b)
  
}

calculateEstimatesTemp <- function(b)
{
  b <- data.frame(b, PgmaxTemp = sapply(as.numeric(rownames(b)), function(x)  
    b[['Pgmax']][x]/b[['PgmaxAdj']][x]))
  b <- data.frame(b, loTemp = sapply(as.numeric(rownames(b)), function(x)  
    b[['lo']][x]/b[['loAdj']][x]))
  b <- data.frame(b, RDTemp = sapply(as.numeric(rownames(b)), function(x)  
    b[['RD']][x]/b[['RDAdj']][x]))
  b <- data.frame(b, Icom = sapply(as.numeric(rownames(b)), function(x)  
    (b[['RDTemp']][x] * b[['PgmaxTemp']][x]) / 
      (b[['loTemp']][x] * (b[['PgmaxTemp']][x] - b[['RDTemp']][x]))))
  

  return(b)
  
}

#' Calculate light saturation point
#' @param b matrix with curve parameters
#' @param I intensity value
calculateEstimatesTempIsat <- function(b, I)
{
  b <- data.frame(b, Isat = sapply(as.numeric(rownames(b)), function(x)
    (b[['RDTemp']][x] * b[['PgmaxTemp']][x] * (I-1)-I* (b[['PgmaxTemp']][x]^2)) /
      (b[['loTemp']][x] * (b[['PgmaxTemp']][x]*(I-1) + 
                             b[['RDTemp']][x]*(I-0.5)))))
  i <- which(colnames(b) == 'Isat')
  colnames(b)[i] <- paste('Isat', (I*100),sep='')
  return(b)

}

#' Quantum yield at specific I
#' @param b matrix with curve parameters
calculateEstimatesTempQuantYield <- function(b)
{
  b <- data.frame(b, oIcomp = sapply(as.numeric(rownames(b)), function(x)
    (b[['loTemp']][x]*(b[['PgmaxTemp']][x]^2)) / 
      ((b[['loTemp']][x] * b[['Icom']][x] + b[['PgmaxTemp']][x])^2)))
    return(b)
}

#'Plot NET PHOTOSYNTHESIS RATE
#' @param v1 curve parameters
#' @param sub subset
#' @param t temperature
plotSelNetTemp <- function(vl, sub, t)
{
  k <-((vl[['loTemp']] *  sub$I * vl[['PgmaxTemp']]) / (vl[['loTemp']] * 
       sub$I + vl[['PgmaxTemp']])) - vl[['RDTemp']]
  plot(sub$I, k, xlab='PHOTOSYNTHETIC PHOTON FLUX DENSITY [mmol (photons) m-2 s-1]
', ylab='NET PHOTOSYNTHESIS RATE', type='l',col='red', 
       main = paste('temperature = ', t, sep =''),
       cex.axis=0.8, cex.lab=0.6)
  
}

#' Plot QUANTUM YIELD
#' @param v1 curve parameters
#' @param sub subset
#' @param t temperature
plotSelQuaTemp <- function(vl, sub, t)
{
  k <-(vl[['loTemp']] *  (vl[['PgmaxTemp']] ^2)) / 
      ((vl[['loTemp']] *  sub$I * vl[['PgmaxTemp']])^2)
  plot(sub$I, k, xlab='PHOTOSYNTHETIC PHOTON FLUX DENSITY [mmol (photons) m-2 s-1]
', ylab='QUANTUM YIELD', type='l', col='red', main = paste('temperature = ', 
                                                           t, sep =''), 
       cex.axis=0.8, cex.lab=0.6)
  
}

r <- read.csv('merge_data.csv', header=TRUE, sep=',')
r$I <- as.numeric(as.character(r$I))
r <- r[complete.cases(r),]
row.names(r) <- 1:nrow(r)
b <- c()

# 
pdf('curves.pdf')
par(mfrow=c(2,2))
for(sp in unique(r$specie))
{
  sub <- r[which(r$specie == sp),]
  for(rep in unique(sub$r))
  {
    sub1 <- sub[which(sub$r == rep),]
    for(temp in unique(sub1$t))
    {
      sub2 <- sub1[which(sub1$t == temp),]
      tname<- paste(unique(as.character(sub2$specie)), 'Temp:', unique(sub2$t),
                   'rep: ', unique(sub2$r), sep=' ')
      b <- rbind(b, nlsLMFunction(sub2$I, sub2$Pn, tname,
                                  unique(as.character(sub2$specie)),
                                 unique(sub2$t), unique(sub2$r), c('1') ))
    }
  }
}



# Transform '-'
b <- processCurveParam(b)

#save(b, file = 'b.dat')
#load('b.dat')


#' Calculate calculate parameters for RD, lo, and Pgmax from b for each species 
#' and rep 
#'TODO look for confidence intervals for Pgmax, lo and RD's parameters (P25
#' and c)
#' 
#' 
# Run loop for both species and one rep
m_all <- data.frame()

par(mfrow=c(2,2))
for(sp in unique(b$specie)){
  for(rept in b$rep[1]){
  print(paste(sp, rept))
   sub <- subset(b, specie == sp & rep == rept)
   
   m <- ArrheniusFit(sub)
   m <- data.frame(rep = rept, species = sp, m)
   m_all <- rbind(m_all, m)
  }
}


break

b <- calculateEstimates(m, b)
b <- calculateEstimatesTemp(b)

# Calculate light saturation point
b <- calculateEstimatesTempIsat(b, 0.50)
b <- calculateEstimatesTempIsat(b, 0.85)
b <- calculateEstimatesTempIsat(b, 0.90)
b <- calculateEstimatesTempIsat(b, 0.95)

#'Quantum yield at specific I
b <- calculateEstimatesTempQuantYield(b)

# 
# Plot NET PHOTOSYNTHESIS RATE
plotSelNetTemp(b[1,], sub2, 25)
plotSelQuaTemp(b[1,], sub2, 25)
dev.off()
# break()
# 
# write.table(b, file='results_curvefitting.csv', sep=',')

       