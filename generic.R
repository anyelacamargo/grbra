library('ggplot2')
library('data.table')
library('grid')

resize_image <- function(output, scale){

image_names = list.files(path = '.', pattern = '*.JPG')

for(image_name in image_names){
  
  image <- load.image(image_name)
  resiz = imresize(image, 0.125)
  #map_il(1,~ resiz) %>% imappend("x") %>% plot
  save.image(resiz, paste('scaled_images/scaled_',
                          strsplit(image_name, '\\.')[[1]][1], '.jpg', sep=''), 
             quality = 0.7)
}

}

plotQTLParents <- function(u, fdr.level = 0.05, 
                           n, chrname, flab, tarchro = flab){
  u <- read.table('dummysample.csv', header=TRUE, sep=',')
  usub <- filter(u, chrom == 'chr2D' | chrom ==  'chr4D')
  input <- usub[, c('marker','chrom', 'pos','A', 'genomewide.permute.pval')]
  input <- mutate(input, pos1 = log(pos))
  
    #pv <- input[,5]
    input <- input[order(input[, 'chrom'], input[, "pos"]), ]
    chroms <- unique(input[, 'chrom'])
    n.chrom <- length(chroms)
    chrom.start <- rep(0, n.chrom)
    chrom.mid <- rep(0, n.chrom)
    if (n.chrom > 1) {
      for (i in 1:(n.chrom - 1)) {
        chrom.start[i + 1] <- chrom.start[i] + max(input[which(input[,
                                                                     'chrom'] == chroms[i]),"pos1"]) + 1
      }
    }
    
    x.max <- chrom.start[n.chrom] + max(input[which(input[,
                                                          'chrom'] == chroms[n.chrom]), "pos1"])
    plot(0, 0, type = "n", xlim = c(0, x.max), ylim = c(0,
                                                        max(input[, 4]) + 1), ylab = "name here", xlab = 'chromosome',
         xaxt = "n", cex.axis=0.8)
    #legend(2,(max(input[, 4]+0.8)), phenotypename)
     for (i in seq(1, n.chrom, by = 2)) {
       ix <- which(input[, 'chrom'] == chroms[i])
       chrom.mid[i] <- median(chrom.start[i] + input[ix,
                                                     "pos1"])
       points(chrom.start[i] + input[ix, "pos1"], input[ix, 'A'],
             col = "dark blue", pch = 16, cex=0.6)
       print(chrom.start[i] + input[ix, 3])
       #if(readline(i) == 'q') { break()}
       abline(v = max(chrom.start[i] + input[ix, "pos1"]), col = 'green')
     }
     if (n.chrom > 1) {
      for (i in seq(2, n.chrom, by = 2)) {
        ix <- which(input[, 'chrom'] == chroms[i])
        chrom.mid[i] <- median(chrom.start[i] + input[ix,
                                                      "pos1"])
        points(chrom.start[i] + input[ix, "pos1"], input[ix, 'A'],
               col = "dark blue", pch = 16, cex=0.6)
     
        abline(v = max(chrom.start[i] + input[ix, "pos1"]), col = 'red')
      }
    }
    
     pv <- which(input[, 5] <= 0.05)
    # for(ch in unique(input[pv, 2])){
    #   
    #   tcpv = which(chroms == ch)
    #   pv1 <- intersect(intersect(which(input[, 5] <= 0.05), which(input[, 2] == ch)),
    #                    which(input[, 4] > 4))
    #   points( chrom.start[tcpv] + input[pv1, 3], input[pv1,4], col = "purple", pch = 16, cex=0.8)
    #   
    # }
    
    # q.ans <- qvalue(10^-input[, 4])
    # 
    # temp <- cbind(q.ans, input[, 4])
    # temp <- temp[order(temp[, 1]), ]
    # if (temp[1, 1] < fdr.level) {
    #   temp2 <- tapply(temp[, 2], temp[, 1], mean)
    #   qvals <- as.numeric(rownames(temp2))
    #   x <- which.min(abs(qvals - fdr.level))
    #   
    #   first <- max(1, x - 2)
    #   last <- min(x + 2, length(qvals))
    #   if ((last - first) < 4) {
    #     last <- first + 3
    #   }
    #   splin <- smooth.spline(x = qvals[first:last], y = temp2[first:last],
    #                          df = 3)
    #   
    #   #lines(x = c(0, x.max), y = rep(predict(splin, x = fdr.level)$y,
    #   #                               2), lty = 2)
    #   #print(rep(predict(splin, x = fdr.level)$y, 2))
    #   #print(max(input[, 4]))
    #   lines(x = c(0, x.max), y = c(4,4), lty = 2)
    #   
    # }
    axis(side = 1, at = chrom.mid, labels = chroms, cex.axis=0.8, las=2)
    d = list()
    d$temp = temp
    #d$fdr = rep(predict(splin, x = fdr.level)$y,2)
    return(temp)
  #}
}
  
  
# @p
reprojection <- function(geodata){


  #data(world.cities)
  #uk = world.cities[world.cities$country.etc == 'UK' & world.cities$pop > 100000,]
  #head(uk)
  # coerce to a spatial object
  #coordinates(uk) = c('long','lat')
  #plot(uk)

  # key proj4 transformations for exercise
  wgs84 = '+proj=longlat +datum=WGS84' #World Geodetic System
  bng = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=
  airy +datum=OSGB36 +units=m +no_defs' #Ordnance Survey Great Britain datum
  mrc = '+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=
  m +nadgrids=@null +wktext +no_defs'

  # useful table of proj4 transformations
  #epsg = make_EPSG()
  #View(epsg[grep("OSGB", epsg$note),])

  #  BNG to WGS84 conversion
  uk_wgs84 = spTransform(geodata, CRS(wgs84))
  head(coordinates(uk_wgs84))

  return(uk_wgs84)
}

combMatrix <- function(N, vec){

  #N   <- 3
  #vec <- c(2, 0)
  lst <- lapply(numeric(N), function(x) vec)
  return(as.matrix(expand.grid(lst)))

}


#' format allele
#' @param x, snp vector
convert.snp2 <- function(x) {
  #print(x)
  #convert to {-1,0,1,NA}
  alleles <- c('AA', 'BB', 'AB'); # 0=AA, 2=BB
  y <- rep(NA,length(x))
  #print(alleles);
  y[which(x==alleles[1])] <- 0
  y[which(x==alleles[2])] <- 2
  y[which(x==alleles[3])] <- 1
  #break();
  return(y)
}


#' process genodata
#' @param dataril, genetic map
#' @param datamap, map
process_rils = function(dataril, datamap)
{
  f = list()
  rils = dataril

  rils_raw_clean = merge(datamap[,1:2], rils, by.x='Marker', by.y='original.order')
  clr = ncol(rils_raw_clean)
  rr = nrow(rils_raw_clean)# Cols dataframe
  f$R <- apply(rils_raw_clean[,4:clr],2, convert.snp2)
  #v = as.vector(rils_raw_clean[2, 3:(clr-1)]);
  markername = as.vector(rils_raw_clean[1:rr, 'Marker'])
  machinename = as.vector(rils[1, 3:(clr-1)])
  machinename = unlist(t(machinename))
  genotypename = as.vector(rils[2, 3:(clr-1)])
  genotypename = unlist(t(genotypename))

  f$rildataV = data.frame(markername, f$R)
  colnames(f$rildataV) = c('Marker', genotypename)
  rownames(f$rildataV) = 1:nrow(f$rildataV)
  f$rildataH = data.frame(genotypename, t(f$R))
  colnames(f$rildataH) = c('genotype', markername)
  rownames(f$rildataH) = 1:nrow(f$rildataH)
  return(f$rildataH)

}

#' getView for multigrid plot
#' @param xloc location of plot on x
#' @param yloc location of plot on y
#' @param h plot height
#' @param y location
#' @param x location
#' @param w width
#'
#' getView()
getView <- function(p, xloc, yloc, h, y, x, w){

  l <- viewport(height = unit(h, "npc"), width = unit(w, "npc"),
                just = c(xloc, yloc),
                y = y, x = x)
  print(p, vp = l)
}

#' createRasterBrick
#' @param pathname
#' @param kw file extension or keyword
#' createRasterBrick()
createRasterBrick = function(pathname, kw){


  filenames <- list.files(path = pathname, pattern = paste('\\.',kw, sep=''))

  RAST = vector("list", 0)
  for (p in filenames) {
    pname <- paste(pathname, p, sep='')
    RAST[[p]] = raster(pname)
  }
  RS = stack(RAST)
  RB = brick(RS, values = TRUE)
  return(RB)
}

#' Test plots
#' @param sdat dataset
#' @cnames colnames
#' plot_test
plot_test = function(sdat, cnames){
  dev.off()
  par(mfrow = c(length(cnames)/2, 2))
  for(cn in cnames){

    plot(sdat[[cn]], main=cn, ylab=cn)
  }

}


# Export traits to for QTL analysis
#' @param data data list
#' @param grdata contains the traits
#' @param tname output file
exportTraitsforQTLRange = function(data, grdata, tname)
{

  wholedata = merge(data$riltable, grdata, by.x='genotype', 'genotype')
  colnames(wholedata)[2] = c('SUBJECT.NAME')
  write.table(wholedata[, -1], file=paste(tname, '.phenotype', sep=''),
              sep='\t', row.names = F, quote = F)

}



#' replaceSmallValue
#' @param sdat dataset
#' replaceSmallValue()
replaceSmallValue <- function(sdat){

  zerovalues <- which(sdat <= 0.01)
  for(z in zerovalues){
    if(z > 3){
      sdat[z] <- sdat[z-1] - 0.01
    }

  }
  return(sdat)
}



goodnessoffit = function(dobs, dpred){

  gof <- 1 - (sum(( dobs - dpred )^2)/
                    sum((dobs - mean(dobs))^2))
  gof1 <-sqrt(mean((dpred - dobs)^2))

  return(gof)

}



optimiseFP <- function(p, X, Y){

  #
  A <- exp(p[1])
  B <- exp(p[2])
  xmid <- exp(p[3])
  scal <- exp(p[4])

  y.pred <- (A-B) /(1 + exp((X-xmid)/scal))
  RSS    <- sum((Y-y.pred)*(Y-y.pred))/length(X)
  return(RSS)
}




#' Fit logist 3 parm model
#' @param sdat data
#' @param trait to be tested
fitLogis = function(sdat, trait){

  y = as.formula(paste(trait, '~', 'SSlogis(DAS, Asym, xmid, scal)', sep=''))
  tmp.logis  <- getInitial(y, data = sdat)
  fit.logis  <- nlsLM(y, trace = F, control = list(maxiter=500), data = sdat)
  return(fit.logis)

}


#' Fit logist 4 parm model
#' @param sdat data
#' @param trait to be tested
fitFpl = function(sdat, trait){

  y = as.formula(paste(trait, '~', 'SSfpl(DAS, A, B, xmid, scal)', sep=''))
  tmp.fpl  <- getInitial(y, data = sdat)
  fit.fpl <- nlsLM(y, data = sdat)
  return(fit.fpl)

}

Init.exp <- function(mCall,LHS,data) {
  xy    <- sortedXyData(mCall[["DAS"]],LHS,data)
  r     <- (coef(lm(log(y) ~ x, xy))[2]) # Use the slope from a linear fit to the data as an initial guess for r
  M0    <- min(xy$y)		             # Use the minimum y value as an initial guess for A
  value <- c(M0, r)
  names(value) <- mCall[c("M0", "r")]
  return(value)
}


#' Fit Gompertz
#' @param sdat data
#' @param trait to be tested
fitGomp = function(sdat, trait){

  y = as.formula(paste(trait, '~', 'SSgompertz(DAS, Asym, b2, b3)', sep=''))
  tmp.gomp <- getInitial(y, data = sdat)
  fit.gomp <- nlsLM(y, data = sdat)
  return(fit.gomp)

}

#' Fit exponential
#' @param sdat data
#' @param trait to be tested
fitExp = function(sdat, trait){

  t = as.formula(paste(trait, '~', 'DAS', sep=''))
  mod1 <- lm(t, data = sdat)
  y = as.formula(paste(trait, '~', 'M0*exp(r*DAS)', sep=''))
  fit.exp <- nlsLM(y, data = sdat, start = list(M0 = exp(coef(mod1)[1]), r = coef(mod1)[2]))
  return(fit.exp)
}

#' fitBimodal
#' @param sdat data set
#' @title plot
#' @param mute produce a plot (TRUE) or not FALSE
fitBimodal = function(sdat, title=NULL, mute=FALSE){

  library(mixtools)
  set.seed(40)
  response = sdat
  mixmdl = normalmixEM(response, maxrestarts=100)
  if(mute == TRUE){
      plot(mixmdl, which=2, main2=title)
      lines(density(response), lty=2, lwd=2)
  }
  return(mixmdl)

}


testModels = function(sdat, trait)
{

  m = list()
  m[['3PLog']] = fitLogis(sdat, trait)
  m[['4PLog']] = fitFpl(sdat, trait)
  m[['Gompz']] = fitGomp(sdat, trait)
  plot(sdat[['DAS']], sdat[[trait]], pch=19, col="grey" , cex.main=1, cex.lab=1.2,
       main=trait, xlab="DAS", ylab="Observations")

  le <- c()
  for(i in names(m)){
    lines(sdat[['DAS']], predict(m[[i]]), col=(which(names(m) == i)+2), lwd=2)
    cod = goodnessoffit(sdat[[trait]], predict(m[[i]]))
    le <- c(le, sprintf("%s-: CoD=%s", i, round(cod, 3)))
  }
  legend("bottomright", legend=le, lwd=2, col=3:6, bty="n", cex=0.8)
}


output.logis.nls <- function(fit, times, CI = F, LOG = F, alpha = 0.05){
  coef <- coef(fit)
  params <- transform_param.logis(coef)
  K = params[1]
  r = params[2]
  M0 = params[3]
  fitted <- fit$m$fitted()
  resid  <- fit$m$resid()
  data <- data.frame(fitted = fitted, resid = resid)
  eq   <- bquote(paste(.(round(M0*r, 4)) /(.(round(M0, 3)) + .(round(K-M0, 2)) * e^{.(round(-r, 3))*t})))
  mss <- sum((fitted - mean(fitted))^2)
  rss <- sum(resid^2)
  R2  <- mss/(mss + rss)
  rmse <- sqrt(rss)
  AIC <- AIC(fit)
  summary <- c(R2 = R2, AIC = AIC, RMSE = rmse)
  rates = data.frame(
    times = times,
    M    = (M0*K)/(M0+(K-M0)*exp(-r*times)),
    AGR  = (r*M0*K*(K-M0)*exp(-r*times))/(M0+(K-M0)*exp(-r*times))^2)
  rates$RGRt <- rates$AGR/rates$M
  rates$RGRm <- r*(1 - rates$M/K)
  if(LOG ==T){
    rates$RGRt <- rates$AGR
    rates$RGRm <- r*rates$M*(1-rates$M/K)
    rates$AGR  <- rates$AGR*exp(rates$M)
  }
  if(CI == T){
    cov   <- summary(fit)$cov
    x <- y <- data.frame(rmvnorm(n=1000, mean=coef, sigma=cov))
    x$K  <- y$Asym
    x$r  <- 1/y$xmid
    x$M0 <- y$Asym/(1 + exp(y$xmid/y$scal)) #untransform best-fit parameters to K, r and M0
    M <- AGR <- RGRt <- RGRm <- matrix(NA, ncol = length(times), nrow = nrow(x))
    for(i in 1:nrow(x)){
      x.i  <- x[i,]
      M[i,]     <- (x.i$M0*x.i$K)/(x.i$M0+(x.i$K-x.i$M0)*exp(-x.i$r*times))
      AGR[i,]   <- (x.i$r*x.i$M0*x.i$K*(x.i$K-x.i$M0)*exp(-x.i$r*times))/(x.i$M0+(x.i$K-x.i$M0)*exp(-x.i$r*times))^2
      RGRt[i,] <- AGR[i,]/M[i,]
      RGRm[i,]  <-  x.i$r*(1 - M[i,]/x.i$K)
      if(LOG ==T){
        RGRt[i,] <- AGR[i,]
        RGRm[i,] <- x.i$r*M[i,]*(1 - M[i,]/x.i$K)
        AGR[i,]  <- AGR[i,]*exp(M[i,])
      }
    }
    CIs <- summarizer(list(M = M, AGR = AGR, RGRt = RGRt, RGRm = RGRm), alpha)
    out <- list(params = params, summary = summary, equation = eq, data = data, rates = cbind(rates, CIs))
  } else {
    out <- list(params = params, summary = summary, equation = eq, data = data, rates = rates)
  }
  return(out)
}



#' get averages
#' @param d the complete dataset
#' @param trait to plot
#' getAVG()
getAVG = function(d, trait){

  dat <- d$traitdata
  dat <- avgtrait(dat, trait, 'DAS')
  colnames(dat) = c('DAS', trait)

  return(dat)
}


#' findCurveType
#' @param sdata dataset
#' @param trait name
#' findCurveType()
findCurveType = function(sdat, trait){

  sdat[[trait]] = log(sdat[[trait]])
  sdat[[trait]] = convertToProp(sdat[[trait]])

  plot(sdat[['DAS']], sdat[[trait]], pch=19, col="grey" , cex.main=1, cex.lab=1.2,
       main=trait, xlab="DAS", ylab="Observations")
  le <- c()
  for(i in 3:4){
    test <- nplr(sdat[['DAS']], sdat[[trait]], npars = i, useLog = FALSE)
    lines(getXcurve(test), getYcurve(test), lwd = 2, col = i)
    goodness <- getGoodness(test)
    gof <- goodness$gof
    le <- c(le, sprintf("%s-P: GoF=%s", i, round(gof, 3)))
  }
  legend("bottomright", legend=le, lwd=2, col=2:5, bty="n", cex=0.8)
}



trcFunc <- function(DAS, L, k, x0, L1, k1, x01){
  (L / (1+exp(-k*(DAS - x0)))) + (L1 / (1+exp(-k1*(DAS - x01))))

}


trcFunc2 <- function(x, A,B,C,D,E){
  A*exp(-((x-B)/C)^2)+D*x+E

}


multimodal = function(){

  library(mixtools)
  wu = sub$water_amount
  mixmdl = normalmixEM(wu)
  plot(mixmdl,which=2)
  lines(density(wu), lty=2, lwd=2)

  library(mixtools)
  wait = faithful$waiting
  mixmdl = normalmixEM(wait)
  plot(mixmdl,density=TRUE)
  lines(density(wait), lty=2, lwd=2)



}

ba <- function(){


  trcFunc <- function(DAS, L, k, x0, L1, k1, x01){
    (L / (1+exp(-k*(DAS - x0)))) + (L1 / (1+exp(-k1*(DAS - x01))))

}




  as.formula(c((L / (1+exp(-k*(DAS - x0)))) + (L1 / (1+exp(-k1*(DAS - x01))))))


  y = as.formula(paste('water_amount', '~', 'trcFunc(DAS, L, k, x0, L1, k1, x01)', sep=''))
  weights = 1/exp(water_amount)

  getInitial(water_amount ~ trcFunc(DAS, L, k, x0, L1, k1, x01), data = sub)

  u = nlsLM(water_amount ~ trcFunc(DAS, L, k, x0, L1, k1, x01),
            data = sub ,start = list(L=0, k=0, x0=-3, L1=0, k1=0, x01=-1), trace = T,
            control = list(maxiter = 500, maxfev=2000))

  u1 = nlsLM(water_amount ~ trcFunc2(DAS, A, B,C,D,E),
             data = sub ,start = list(A=0.35, B=130, C=20, D=0, E=0), trace = T,
             control = list(maxiter = 500, maxfev=2000))


  plot(sub$water_amount, col='red')
  lines(predict(u1), col='blue')
  lines(predict(u), col='green')
}


# Smooth curve
#' @param data is the trait data
#' @param trait to be plotted
#' @param predictor is the predictor
SmoothCurve = function(dat, trait, predictor) {

  f = as.formula(paste(trait, '~', predictor))
  loessFit <- loess(f, data = dat, span = 0.35)
  loessFit <- data.frame(DAS=loessFit$x, trait=loessFit$fitted)
  loessFit <- loessFit[order(loessFit$DAS),]
  colnames(loessFit) = c(predictor, trait)
  #approxFit <- approx(u,n = 15)
  #lowessFit <-data.frame(lowess(u,f = .6,iter=1))
  #colnames(lowessFit) = colnames(loessFit)

  #plot(u,col='gray')
  #curve(sigmoid,0,200,add=TRUE,col='blue',)
  #lines(lowessFit,col='red')
  #lines(loessFit,col='green')
  #lines(approxFit,col='purple')
  #plot(u, las=2, cex.lab = 1.5, cex.axis=1.5,
  #            col=4, ylab =NA, xaxt="n", xlab = group, type='l', lwd=2)
  #      axis(1, at=x,labels=x, las=2, cex.axis=1.5)

  #legend(150,.6,
  #      legend=c("Sigmoid","Loess","Lowess",'Approx'),
  #       lty=c(1,1),
  #      lwd=c(2.5,2.5),col=c("blue","green","red","purple"))
  return(loessFit)
}


#' plotSmoothCurve function
#' @param dat, the data
#' printSmoothCurve()
plotSmoothCurve = function(dat, trait, group){

  x = unique(dat[['DAS']])[seq(1, length(unique(dat[['DAS']])), by=4)]
  par(oma=c(2,1,0,0))
  plot(dat, las=2, cex.lab = 1.5, cex.axis=1.5,
       col=4, ylab =NA, xlab = group, xaxt="n", type='l', lwd=2)
  axis(1, at=x, labels = x, las=2, cex.axis=1.5)

}

#' Get avg trait curve
#' @param data = traitdata
#' @param trait = trait
#' @param group = grouping factor
avgtrait = function(dat, trait, group)
{
  i = which(dat[[trait]] <= 0.01)
  if(length(i) > 0){
    dat = dat[-i, ]
  }

  u = aggregate(dat[[trait]], by=list(dat[[group]]), FUN=mean)
  return(u)
}




# Function that returns Root Mean Squared Error
rmse <- function(error)
{
  sqrt(mean(error^2))
}

# Function that returns Mean Absolute Error
mae <- function(error)
{
  mean(abs(error))
}


# Read excel files
readExcel = function(filename, indexname=NULL, nv=NULL)
{
  library('xlsx')
  library('rJava')
  library('readxl')
  #library('openxlsx')

  wdata = list()
  for(i in indexname)
  {

    xlsx_in = read_excel(filename, sheet = i[['index']], range = cell_cols(i[['range']]),
                         col_types = NULL, na=nv)
    wdata = rbind(wdata, xlsx_in)
  }
  return(wdata)
}

# Convert to DOY


convertDOY = function(data, fname, fdate, cname)
{

  data[[fname]] = as.Date(data[[fname]], format = fdate)
  data = data.frame(data, tname = sapply(data[[fname]],
                                         function(x) as.numeric(strftime(x, format = "%j"))))
  i = which(colnames(data) == 'tname')
  colnames(data)[i] = cname
  return(data)

}


converIdtag = function(x)
{
  if(x < 10000)
  {
    p = paste('W8-','00',x, sep='')
  }

  else if(x >= 10000 & x < 100000)
  {
    p = paste('W8-','0',x, sep='')
  }

  else
  {
    p = paste('W8-',x, sep='')
  }

  return(p)
}

#' stackFiles
#' @param patternword keyword
#' use ^ to restrict the keyword
stackFiles = function(kw, pathname)
{
  files <- list.files(path = pathname, pattern = paste('\\.',kw, sep=''))

  return(files)
}


stackFiles1 = function(kw, pathname)
{
  files <- list.files(path = pathname, pattern = kw)

  return(files)
}

# stackDataSets = function(filelist, sk=NULL, h=TRUE)
# {
#
#   myfiles = do.call(rbind, lapply(filelist,
#                                   function(x) cbind(x, read.table(x, stringsAsFactors = TRUE,
#                                                                   header=h))))
# }


#' stackDataSets
#' @param filelist of files to stack
#'
stackDataSets = function(filelist, sk=0, h=TRUE)
{

  myfiles = do.call(rbind, lapply(filelist,
                                  function(x) cbind(x, read.table(x, stringsAsFactors = TRUE,
                                                                  header=h, skip=sk))))
}

stackDataSetsExcel = function(filelist, wslist, navalue = '---')
{

  myfiles = do.call(rbind, lapply(filelist,
                                  function(x) cbind(x, readExcel(x, wslist, navalue))))
}


loadPheno = function(file_list, h, sepc, cl=NULL, b)
{

  f = read.table(file_list[1], sep=sepc, header=h)
  for(fname in file_list[2:length(file_list)])
  {

    f1 =  read.table(fname, sep=sepc, header=F)
    f = merge(f, f1, by=b, by.y=b)

  }
  return(f)
}


loadMarkerFile = function(file_list, h, sepc, traitname)
{

  f = read.table(file_list[1], sep=sepc, header=h)
  f = f[, c(1,8,6)]


  for(fname in file_list)
  {
    f1 =  read.table(fname, sep=sepc, header=h)
    colnames(f1)[9] = strsplit(fname, '\\.')[[1]][1]
    f = merge(f, f1[,c(1,8,9,11)], by.x=c('from', 'chromosome'),
              by.y=c('from', 'chromosome'))
  }

  colnames(f)[1:4] = c('marker', 'chrom',	'pos', traitname)
  return(f)
}


# w = what (data.frame)
# b = by: Should be a list
averageValues = function(w, b, funtype='mean')
{
  copydata = aggregate(w, by = b,
                       FUN= "mean", na.rm = TRUE)

  return(copydata)
}


# Merge files
datatamerge_files <- function(p = 'txt')
{
  filenames <- list.files(pattern = p)
  c <- do.call("rbind", lapply(filenames, read.csv, header = F))
  c <- read.table('output.csv', sep=',', header=F)
}


# Function to read data from external files
read_data = function(datafilename, h, sepc)
{
  idata <- read.table(datafilename, sep=sepc, header=h, na.strings = "NA", stringsAsFactors=TRUE)

  return(idata)

}


# Change colname
change_colname = function(data, oldp, newp)
{
  copydata  = data
  setnames(copydata, oldp, newp)
  return(copydata)

}

# Change row value
change_rowvalue = function(data, fieldname, oldp, newp )
{
  copydata  = data
  copydata[[fieldname]][copydata[[fieldname]] == oldp ] = newp
  return(copydata)
}


## Create a plot, passes character strings. Fit a line

create_figure= function(data, attrib_list, data2)
{
  copydata = data
  p = ggplot(na.omit(copydata), aes_string(x=attrib_list[['aes.x']], y=attrib_list[['aes.y']],
                                           group = 1)) +
    geom_point(shape=1) +    # Use hollow circles
    geom_smooth(method=lm) +
    geom_point(data=data2, aes_string(y=attrib_list[['aes.y']], x=attrib_list[['aes.x']]),
               color='red', size = 3) +
    ggtitle(attrib_list[['x.title']]) +
    xlab('DAS') + ylab(attrib_list[['y.title']]) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))


  return(p)
}

calibratedata = function(data)
{
  copydata = data
  for(tdate in unique(copydata$time))
  {
    t = which(copydata$time == tdate)
    for(i in 3:(ncol(copydata)-2))
    {
      if(as.Date(tdate) <= as.Date("2015-04-01"))
      {
        copydata[t,i] = copydata[t,i]*0.27
      }
      else
      {
        copydata[t,i] = copydata[t,i]*0.67 # 0.56
      }
    }
  }
  return(copydata)
}


convert2numeric = function(data, col_list)
{
  copydata=data
  for(cname in col_list)
  {
    copydata[[cname]] = as.numeric(as.character(copydata[[cname]]))
  }
  return(copydata)
}

# Convert to number
convert2number = function(data, col_list)
{
  copydata=data
  for(cname in col_list)
  {
    copydata[[cname]] = as.numeric(copydata[[cname]])
  }
  return(copydata)
}

# subset of data
average_similar_measures = function(data, featurename)
{
  copydata = data
  copydata[[featurename]] = factor(copydata[[featurename]])
  #print(plantid)
  # Arrange new dataset
  l <- split(copydata, copydata[[featurename]])
  h <- do.call(rbind, lapply(1:length(l), function(x) takeMeans(l[x])))
  h = h[,c(-1,-2)]
  h = h[,c(dim(h)[2], 1:(dim(h)[2]-1))]
  h <- data.frame(idtag=unique(copydata$idtag), h)
  colnames(h)[3:(dim(h)[2])] <- colnames(scsub)[3:(dim(scsub)[2])]
  return(h)
}

# Create plot but adding a vertical line to indicate something
## TODO: MAKE It GENERIC
print_plot = function(data, ldate, ytrait, xtrait)
{

  #days before phenotyping

  dbp = seq(as.Date("2014-10-25"), as.Date("2015-04-27"), by="days")
  i = which(as.Date(unique(data$time)) == as.Date(dbp[ldate]))
  print(i)
  if(length(i) == 0) { i = which(as.Date(unique(data$time)) == as.Date(dbp[ldate+1]))}
  copydata = data
  p = ggplot(copydata,aes(x = time, y = value, fill=intensity, group=intensity)) +
    geom_area() + geom_line(aes(ymax=value), position="stack") +
    geom_vline(xintercept = i, colour='red', size=2) +
    scale_fill_manual(values = as.character(unique(copydata$intensity))) +
    labs(y = ytrait, x = xtrait ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          text = element_text(size=8)) +
    ggtitle(idtagname) +
    #theme(axis.text.x = element_blank() ) + theme(legend.position="none")
    theme(axis.text.x = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  return(p)
}


sum_intensities <- function(data)
{

  copydata = data
  copydata = data.frame(copydata, sumint = sapply(rownames(copydata),
                                                  function(x) sum(copydata[x,3:131])))
  return(copydata)

}

#Select subset of data
selectSubset = function(data, fieldname, targetname)
{
  return(data[which(data[[fieldname]]  == targetname),])
}

convert2hex <- function(data1, data2)
{
  copydata1 = data1
  copydata2 = data2
  hexlist = c()

  for(i in 1:length(copydata1))
  {
    idx = which(copydata2$i == copydata1[i])
    hexlist = append(hexlist, as.character(copydata2[idx,5]))

  }
  return(hexlist)
}


get_subdata = function(data, conver_table, intens)
{


  copydata = data
  # I us want to use the same intensities # Vey crude way but it seems they're the most relevants.

  u = colnames(copydata) %in% intens # Search for intens
  scc1= data.frame(idtag = copydata$idtag,
                   time = copydata$time, copydata[,u]) # Create frame with sel intensi
  freq = colnames(scc1)[3:dim(scc1)[2]]
  freq = as.numeric(sub("X", "", freq))
  colnames(scc1)[3:dim(scc1)[2]] = freq # Delete X in columns
  #write.table(scc1, file='t.csv', sep=',', row.names=F)

  mm = as.matrix(scc1[, 3:dim(scc1)[2]]) # Deletes -4
  vector_data = unmatrix(mm,byrow=T) # Convert matrix into vector
  vector_data = as.numeric(vector_data) # Convert values to numeric
  idtag = as.character(unique(scc1$idtag)) # Unique idtags
  time = as.character(unique(scc1$time)) # unique time points
  hexlist = convert2hex(freq, conver_table) #
  conveg = expand.grid(idtag=idtag, intensity=hexlist,time=time)
  cdata = data.frame(conveg, value=vector_data)
  #cdatacal = calibrate_data(cdata)
  return(cdata)
}

# Take row means

takeMeans <- function(x)
{
  x1 = rapply(x, mean)
  x1['time'] = names(x)
  return(x1)
}

# subset of data
get_subset = function(data)
{
  copydata = data
  #print(plantid)
  # Arrange new dataset
  l <- split(copydata, copydata$timestamp)
  h <- do.call(rbind, lapply(1:length(l), function(x) takeMeans(l[x])))
  h = h[,c(-1,-2)]
  h = h[,c(dim(h)[2], 1:(dim(h)[2]-1))]
  h <- data.frame(idtag=unique(copydata$idtag), h)
  colnames(h)[3:(dim(h)[2])] <- colnames(scsub)[3:(dim(scsub)[2])]
  return(h)
}

### DB

connect2DB = function(query_list, db_list)
{
  library('RPostgreSQL')
  library('RODBC')
  #library('RMySQL')
  install.packages('XLConnect')



  drv <- dbDriver(db_list[['dbdriver']])
  dbhost <- db_list[['dbip']]
  dbport <- db_list[['dbport']]
  con <- dbConnect(drv, host=dbhost, port=dbport, dbname=db_list[['dbname']],
                   user=db_list[['dbuser']], password=db_list[['dbpass']])
  rs = dbSendQuery(con, query_list)
  data <- fetch(rs, n = -1)
  dbDisconnect(con)
  return(data)
}

# Create formula
creafeformula = function()
{

  formul = as.formula(paste('y ~', paste(rownames(data.frame(summary(qtl))), collapse='+')))

}


#Export data in weka format

create_query_temp = function(tdate)
{
  s = '%'
  st = sprintf("SELECT logdate, s3_1, s3_2,s3_3,s3_4,s3_5,s3_6,s3_7
               FROM phenomics.glasshouse_sensors where logdate > '%s' and compartent = '5' ",tdate)
  return(gsub("\r?\n|\r", " ", st))
}

convert2factor = function(data, variable_list)
{

  copydata = data
  for(v in variable_list)
  {
    copydata[[v]]  = factor(copydata[[v]])
  }
  return(copydata)
}

#Date format
changeDateFormat = function(data, fname)
{
  copydata = data
  copydata[[fname]] = as.Date(copydata[[fname]],format = "%Y-%m-%d")
  copydata[[fname]] = format(copydata[[fname]], "%d/%m/%Y")
  return(copydata)

}

# Convert from DOY to Date
convert2Date = function(x, or="2016-01-01", fmt='%Y-%m-%d')
{
  d = as.Date(x, format = fmt, origin = or)

  return(d)

}

createDAS = function(ds="2014-10-20", de="2015-04-30")
{
  dv = seq(as.Date(ds), as.Date(de), by="days")
  dv = format(dv, "%d/%m/%Y")
  dv = data.frame(date=dv, DAS = 1:length(dv))
  return(dv)

}

resetlevels = function(data, fname)
{
  copydata = data
  copydata[[fname]] = factor(copydata[[fname]], levels=unique(copydata[[fname]]))
  row.names(copydata) = 1:nrow(copydata)
  return(copydata)
}



# split timestamp into date and time_ Lemnatec
process_timestamp = function(data, timename, sp, timeaction)
{
  copydata = data
  d <- data.frame(copydata, date=sapply(copydata[[timename]],
                                        function(x) strsplit(as.character(x), sp)[[1]][1]))
  d <- data.frame(d, time=sapply(d[[timename]],
                                 function(x) strsplit(strsplit(as.character(x), ' ')[[1]][2],':')[[1]][1]))
  i = which(colnames(d) == 'time')
  colnames(d)[i] = timeaction

  return(d)
}

# Merge two datasets
merge_data <- function(d1, d2, genoname, gname)
{

  i = which(colnames(d2) == genoname)
  colnames(d2)[i] = 'genotypename'

  print(colnames(d2))
  d <- data.frame(d1, genotypename=sapply(d1$id_tag, function(x)
    d2[which(d2$barcode == as.character(x)), colnames(d2)[i]]))
  d <- data.frame(d, water.treatment=sapply(d$id_tag, function(x)
    d2[which(d2$barcode == as.character(x)), 'control']))
  d <- data.frame(d, rep=sapply(d$id_tag, function(x)
    d2[which(d2$barcode == as.character(x)), 'rep']))
  return(d)
}


#
# Calculate growth rate
# data
calculateGrowthRate <- function(data, genotype_list, DAS_list, trait)
{
  rownames(data) = 1:nrow(data)
  copydata = data
  copydata = data.frame(copydata, GR=rep(0, nrow(copydata)))

  for(gname in genotype_list)
  {
    for(dname in 2:length(DAS_list))
    {
      np <- intersect(which(copydata[['genotype']] == gname), which(copydata[['DAS']] == DAS_list[dname-1]))
      nn <- intersect(which(copydata[['genotype']] == gname), which(copydata[['DAS']] == DAS_list[dname]))
      if(!is.na(data[np, trait]) & !is.na(data[nn, trait]))
      {
        copydata[nn, 'GR'] <- log(data[nn, trait]) - log(data[np, trait]) / ( DAS_list[dname] - DAS_list[dname-1])
        print(copydata[nn, 'GR'])
      }
    }
  }



  return(copydata)
}



scaleRange = function(min)
{
  (x-min(x))/(max(x)-min(x))
}


growthModels = function()
{
  fm <- lme(Area.in.square.mm ~ DAS, ss, random = ~ DAS | genotype)
  fm1 <- lme(log(Area.in.square.mm) ~ DAS,  ss, random = ~ DAS | genotype)
  fm2 <- lme(log(Area.in.square.mm) ~ DAS + I(DAS^2), ss, random = ~ DAS | genotype)
  fm3 <- lme(log(Area.in.square.mm) ~ DAS + I(DAS^2) + I(DAS^3), ss, random = ~ DAS | genotype)

  plot(fm3, log(Area.in.square.mm) ~ fitted(.) | genotype, abline = c(0,1))

  u =predict(fm3, ss, level = 0:1)
  with(ss, plot(log(Area.in.square.mm) ~ DAS))
  points(ss$DAS, u$predict.genotype, col='red', lwd=5)

}


# b biomass w water n = index
wue <- function(data, genotype_list, das_list)
{
  copydata <- data
  #time_list <- c(163, 164,165,166)#unique(datacopy$DAS)
  #line_list <- as.character(unique(datacopy$genotype))
  #line_list1 = c('MEL 046-7')

  for(line in genotype_list)
  {
    for(time in 2:(length(das_list)-1))
    {
      nb <- intersect(which(datacopy[['genotype']] == line), which(datacopy[['DAS']] == das_list[time-1]))
      na <- intersect(which(datacopy[['genotype']] == line), which(datacopy[['DAS']] == das_list[time]))
      bb = datacopy[['Area.in.square.mm']][nb]
      ba = datacopy[['Area.in.square.mm']][na]
      wb = datacopy[['water_amount']][nb]
      wa = datacopy[['water_amount']][na]

      if (ba < bb){
        datacopy[['Area.in.square.mm']][na] = bb
        ba = datacopy[['Area.in.square.mm']][na]
      }
      dtb = ba/1000 - bb/1000
      dtw = wa - wb
      if(dtw <= 0)
        dtw = 1

      wue = dtb / dtw

      datacopy[na, 'wue'] <- round(wue, 2)
      #if(wue == Inf | is.na(wue))
      #{
      #print(paste(line, time_list[time], nb, bb, ba, wb, wa, wue, sep=': '))
      #}
    }
  }
  datacopy = datacopy[-(which(datacopy[['DAS']] == time_list[1])),]
  return(datacopy)
}


# fname = 'genotype'
# genotype_list = unique(o$genotype)
#das_list = unique(o$DAS)

WUEDaily = function(data)
{
  wuetable = data
  wuetable$water_amount[wuetable$water_amount == 0.01] = NA
  wuetable = wuetable[!is.na(wuetable$water_amount),]
  wuetable = data.frame(wuetable,
                        WUE = sapply(row.names(wuetable),
                                     function(x) wuetable[x,'Area.in.square.mm']/wuetable[x,'water_amount'])) # area/water_amount
  m = which(wuetable$DAS <= 101)
  wuetable = wuetable[-m,] # delete rows whose time is <= 2015-05-11
  m = which(wuetable$DAS <= 112)
  wuetable = wuetable[-m,] #
  return(wuetable)

}



deleteUncomplete = function(genotype_list, das_list, data, fname)
{
  glist = c()
  copydata = o
  for(dname in das_list)
  {
    subdata = copydata[which(copydata[[fname[1]]] == dname),]
    if(length(unique(subdata[[fname[2]]])) !=  length(genotype_list))
    {
      glist = append(glist, dname)

    }

  }
  i = match(copydata$DAS, glist)
  copydata = copydata[is.na(i),]

  return(copydata)
}



# # ------------- Read metadata
#
#
# tidy_cols = function(data, clist)
# {
#   copydata = data
#   ttc = clist['X1']
#   i = match(ttc, colnames(copydata))
#   colnames(copydata)[i] =  clist['X2']
#   return(copydata)
#
# }
#
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}
#
#
# # Plot water use efficiency
plot_wue = function(data)
{
  copydata = data
  for(ename in unique(copydata$genotypename))
  {
    fname = paste("Intplot_AreaPerGenotype_WUE", ename, ".png", sep="")

    i = which(copydata$genotype == ename)
    sdata = copydata[i, ]
    pd <- position_dodge(0.1) # move them .05 to the left and right
    tgc <- summarySE(sdata, measurevar="wue",
                     groupvars=c('water.treatment', 'genotypename', 'date', 'temp'))
    png(fname)
    mn = min(sdata$temp)
    mx = max(sdata$temp)
    p = ggplot(data = tgc,
               aes(x = date, y = wue, colour=temp, shape=water.treatment,
                   group=interaction(water.treatment, genotypename))) +
      ggtitle(ename) +
      geom_path(alpha = 0.5) +
      geom_errorbar(aes(ymin=wue-se, ymax=wue+se), width=.1, position=pd) +
      stat_summary(fun.y=mean, geom="point", size = 3)+
      stat_summary(fun.y=mean, geom="line") +
      theme(axis.text.x=element_text(angle=90),
            axis.text.x = element_text(size=20)) +
      scale_colour_gradient(limits=c(mn, mx), low="yellow", high="red")
    print(p)
    #if(readline(fname) == 'q') { break()}
    dev.off()
  }
}
#
# Plot temperature
plot_temp = function(data)
{
  copydata = data
  for(ename in unique(copydata$genotypename))
  {
    fname = paste("Intplot_ArePerGenotype_temp", ename, ".png", sep="")

    i = which(copydata$genotype == ename)
    sdata = copydata[i, ]
    pd <- position_dodge(0.1) # move them .05 to the left and right
    tgc <- summarySE(sdata, measurevar="area", groupvars=c('water.treatment',
                                                           'genotypename', 'date', 'temp'))
    png(fname)
    mn = min(sdata$temp)
    mx = max(sdata$temp)
    p = ggplot(data = tgc, aes(x=date, y=area, colour=temp, shape=water.treatment, group=
                                 interaction(water.treatment, temp))) +
      ggtitle(ename) + geom_path(alpha = 0.5) +
      geom_errorbar(aes(ymin=area-se, ymax=area+se), width=.1, position=pd) +
      stat_summary(fun.y=mean, geom="point", size=4) +
      stat_summary(fun.y=mean, geom="line") +
      theme(axis.text.x=element_text(angle=90),
            axis.text.x = element_text(size=20)) +
      scale_colour_gradient(limits=c(mn, mx), low="yellow", high="red")
    print(p)
    #if(readline(fname) == 'q') { break()}
    dev.off()
  }
}

#
#average data
avgdata <- function(data, fields, variables)
{
  print(fields)
  print(variables)
  copydata <- data
  copydata = aggregate(cbind(fields) ~ variables, data = copydata, FUN= "mean" )
  return(copydata)
}

#
# # Split temp data by date by hour
# splittime = function(data)
# {
#
#   copydata = data
#   d = strsplit(copydata, ' ')[[1]][1]
#   h = strsplit(strsplit(copydata, ' ')[[1]][2], ':')[[1]][1]
#   return(paste(d,h))
# }
#

PlotCorrelation = function(data)
{
  dta = data
  # Plot correlations
  dta.r <- cor(dta, use = 'complete.obs') # get correlations
  dta.col <- dmat.color(dta.r) # get colors
  dta.o <- order.single(dta.r)
  cpairs(dta, dta.o, panel.colors=dta.col, gap=.5, cex.labels=1.2)

}

pcacomplete <- function(data)
{
  copydata <- data
  res <- PCA(copydata, quali.sup=c(1:5), graph=F, scale.unit=TRUE)
  #png(fname, width = 2000, height = 880, pointsize = 22)
  plot(res, choix = "var", cex=0.9, label='var',  col.var = rainbow(10), lwd=2,
       font.lab=2, font=2, ylab='hh')
  #dev.off()
  plotloadings(res$var$contrib, time)
  #lo <- sweep(res$var$coord,2,sqrt(res$eig[1:ncol(res$var$coord),1]),FUN="/")
  #plotloadings(lo, time)
  res <- PCA(copydata[6:25], graph=F)
  u = dimdesc(res, axes=c(1,2))
  return(res)

}

#' Check file list
#' @param filelist: List of files
#' @param pathname= where files are
checkFileSize = function(filelist, pathname)
{
  m = c()
  for(i in 1:length(filelist))
  {

    fname = paste(pathname, filelist[i], sep='')
    if(file.info(fname)$size <= 2)
    {
      m = append(m, i)
    }
    filelist[i] = fname
  }
  if(!is.null(m)){
    filelist = filelist[-m]
  }
  return(filelist)
}


# order dataset
orderData = function(data, fieldnames)
{


  data = data[order(data[, 'Field1'], data[, 'Field2']), ]


}



manhattan <- function(input, fdr.level = 0.05, phenotypename, chromsl=NULL, flab=1, tarchro=NULL) {
  #pv <- input[,5]
  input <- input[order(input[, 2], input[, 3]), ]
  chroms <- unique(input[, 2])
  n.chrom <- length(chroms)
  chrom.start <- rep(0, n.chrom)
  chrom.mid <- rep(0, n.chrom)
  if (n.chrom > 1) {
    for (i in 1:(n.chrom - 1)) {
      chrom.start[i + 1] <- chrom.start[i] + max(input[which(input[,
                                                                   2] == chroms[i]), 3]) + 1
    }
  }
  x.max <- chrom.start[n.chrom] + max(input[which(input[,
                                                        2] == chroms[n.chrom]), 3])
  plot(0, 0, type = "n", xlim = c(0, x.max), ylim = c(0,
                                                      max(input[, 4]) + 1), ylab = "-log(p)", xlab = 'chromosome',
       xaxt = "n", cex.axis=0.8)
  legend(2,(max(input[, 4]+0.8)), phenotypename)
  for (i in seq(1, n.chrom, by = 2)) {
    ix <- which(input[, 2] == chroms[i])
    chrom.mid[i] <- median(chrom.start[i] + input[ix,
                                                  3])
    points(chrom.start[i] + input[ix, 3], input[ix, 4],
           col = "dark blue", pch = 16, cex=0.6)
  }
  if (n.chrom > 1) {
    for (i in seq(2, n.chrom, by = 2)) {
      ix <- which(input[, 2] == chroms[i])
      chrom.mid[i] <- median(chrom.start[i] + input[ix,
                                                    3])
      points(chrom.start[i] + input[ix, 3], input[ix,
                                                  4], col = "red", pch = 16, cex=0.6)
    }
  }


  if(!is.null(tarchro))
  {
    # New line to add color
    ix <- which(input[, 2] == tarchro)

    tc = which(chroms == tarchro)

    chrom.mid[tc] <- median(chrom.start[tc] + input[ix,3])
    points(chrom.start[tc] + input[ix, 3], input[ix,4], col = "green", pch = 16, cex=0.6)

  }


  pv <- which(input[, 5] <= 0.05)
  for(ch in unique(input[pv, 2])){

    tcpv = which(chroms == ch)
    pv1 <- intersect(intersect(which(input[, 5] <= 0.05), which(input[, 2] == ch)),
                     which(input[, 4] > 4))
    points( chrom.start[tcpv] + input[pv1, 3], input[pv1,4], col = "purple", pch = 16, cex=0.8)

  }

  q.ans <- qvalue(10^-input[, 4])

  temp <- cbind(q.ans, input[, 4])
  temp <- temp[order(temp[, 1]), ]
  if (temp[1, 1] < fdr.level) {
    temp2 <- tapply(temp[, 2], temp[, 1], mean)
    qvals <- as.numeric(rownames(temp2))
    x <- which.min(abs(qvals - fdr.level))

    first <- max(1, x - 2)
    last <- min(x + 2, length(qvals))
    if ((last - first) < 4) {
      last <- first + 3
    }
    splin <- smooth.spline(x = qvals[first:last], y = temp2[first:last],
                           df = 3)

    #lines(x = c(0, x.max), y = rep(predict(splin, x = fdr.level)$y,
    #                               2), lty = 2)
    #print(rep(predict(splin, x = fdr.level)$y, 2))
    #print(max(input[, 4]))
    lines(x = c(0, x.max), y = c(4,4), lty = 2)

  }
  axis(side = 1, at = chrom.mid, labels = chromsl, cex.axis=flab, las=2)
  d = list()
  d$temp = temp
  #d$fdr = rep(predict(splin, x = fdr.level)$y,2)
  return(temp)
}


qvalue <- function(p)
{
  smooth.df = 3
  if (min(p) < 0 || max(p) > 1) {
    print("ERROR: p-values not in valid range.")
    return(0)
  }
  lambda = seq(0, 0.9, 0.05)
  m <- length(p)
  pi0 <- rep(0, length(lambda))
  for (i in 1:length(lambda)) {
    pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
  }
  spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
  pi0 <- predict(spi0, x = max(lambda))$y
  pi0 <- min(pi0, 1)
  if (pi0 <= 0) {
    print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values.")
    return(0)
  }
  u <- order(p)
  qvalue.rank <- function(x) {
    idx <- sort.list(x)
    fc <- factor(x)
    nl <- length(levels(fc))
    bin <- as.integer(fc)
    tbl <- tabulate(bin)
    cs <- cumsum(tbl)
    tbl <- rep(cs, tbl)
    tbl[idx] <- tbl
    return(tbl)
  }
  v <- qvalue.rank(p)
  qvalue <- pi0 * m * p/v
  qvalue[u[m]] <- min(qvalue[u[m]], 1)
  for (i in (m - 1):1) {
    qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]],
                        1)
  }
  return(qvalue)
}
