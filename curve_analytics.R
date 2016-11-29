library('minpack.lm');
library("nlmrt")
library("optimx")
library(ggplot2);

source('../senescence_disease/generic.R');

# [64-bit] C:\Program Files\R\R-3.3.1

rhs1 <- function(b, mydata) {
  sum((mydata$y - ((b[1]*mydata$x*b[2]) / (b[1]*mydata$x+b[2])) - b[3])^2)
}


rhs2 <- function(b, mydata) {
  #eunsc = ydat ~ ((tdat * Pgmax)/(tdat + lo)) - RD;
  sum((mydata$y - (((mydata$x*b[2]) / (mydata$x+b[1])) - b[3]))^2)
}

rhsfit1 <- function(b, mydata) {
  #ydat ~ P25*exp((ha/(8.314*298.15)) - ha/(8.314*tdat))
  sum((mydata$y - (b[1]*exp((b[2]/(8.314*298.15)) - b[2] /(8.314*mydata$x)) ) )^2)
}

rhsfit2 <- function(b, mydata) {
  #(P25*exp( (log(1 + exp((S*298.15-hd)/(8.314*298.15))) + (ha/(8.314*298.15))) 
   #         - (ha/(8.314*tdat)))) / (1 + exp((S*tdat-hd)/(8.312*tdat)));
  sum((mydata$y - (  (b[1]*exp(( log(1+ exp((b[3]*298.15-b[4])/(8.314*298.15))) + (b[2]/(8.314*298.15)) 
                               - (b[2]/(8.314*mydata$x)))) / (1 + exp((b[3]*mydata$x-b[4])/(8.314*mydata$x)))  )) )^2)
}

optima = function(start1, w, fc)
{
  ## Optim
  w.optx <- optimx(par=start1, fn=fc, mydata=w, 
                   control=list(all.methods=TRUE, save.failures=TRUE, maxit=10000))
  return(w.optx)
}
  
plotCurve = function(x,y, fm, eq, yl, spc, tempr, repl)
{
  plot(x, y, main = paste('eq: ', eq, sep = ''), ylab = yl,
       xlab = 'PHOTOSYNTHETIC PHOTON FLUX DENSITY', cex.axis=0.7);
  lines(x, fm, col = 2, lwd = 2);
  text(max(x)/2, max(y)/2, paste('specie: ', spc,
    '\ntemp:', tempr, '\n',  'rep:', repl, sep=''), cex = 0.9);
}


searchOutlier = function(data)
{
  library('outliers');
  data = data[order(x),];
  row.names(data) = 1:nrow(data);
  u = chisq.out.test(data$y, variance=var(data$y), opposite = FALSE);
  i = which(data$y == as.numeric(strsplit(u$alternative, ' ')[[1]][3]));
  data = data[-i,];
  row.names(data) = 1:nrow(data);
  s = c();
  mnn = round(data$y[1],3);
  mno=mnn;
  for(i in 2:nrow(data))
  {
    if(mno > mnn)
    {
      mno = mno;
    }
    else
    {
      mno = mnn;
    }

    if(mno > round(data$y[i],3) )
    {
  
      s = append(s, i);
      mnn = round(data$y[i-1],3);
    }
    else
    {
      mnn = round(data$y[i],3);
    }
    
  }
  
  data = data[-s,];
  row.names(data) = 1:nrow(data);
  return(data);
}


eqFitcl1 = function(tdat, ydat, weeddata)
{
  o = list();
  weeddata = data.frame(y=ydat, x=tdat);
  o$start= c(P25 = 0, ha=38000);
  o$eunsc = ydat ~ P25*exp((ha/(8.314*298.15)) - ha/(8.314*tdat));
  o$optx1 = optima(o$start, weeddata, rhsfit1);
  o1= nlsLM(o$eunsc, 
              start=list(P25 = o$optx1$P25[1], ha = o$optx1$ha[1]), trace = TRUE);
  return(o1);
}

eqFitcl2 = function(tdat, ydat, weeddata)
{
  o = list();
  o$start= c(P25 = 10, ha=20000, S=800, hd=258000);
  o$eunsc = ydat ~ (P25*exp( (log(1 + exp((S*298.15-hd)/(8.314*298.15))) + (ha/(8.314*298.15))) 
                             - (ha/(8.314*tdat)))) / (1 + exp((S*tdat-hd)/(8.312*tdat)));
  o$optx1 = optima(o$start, weeddata, rhsfit2);
  o2 = nlsLM(o$eunsc, 
               start=list(P25 = o$optx$P25[1], ha = o$optx$ha[1], S = o$optx$S[1], hd = o$optx$hd[1]), 
               trace = TRUE);
  return(o2);
}


getlogistic = function(x, y)
{
  library('nlstools');
  library('nlme');
  weeddata = data.frame(y=y, x=x);
  f = as.formula(paste('y', '~', 'SSlogis(x, phi1, phi2, phi3)', sep=''));
  gr = nls(f, data = weeddata);
  r = nlsResiduals(gr);
  i = union(which(r$resi2[1:nrow(weeddata),2] > 1), which(r$resi2[1:nrow(weeddata),2] < -1));
  # plot(weeddata[['x']], weeddata[['y']], xlab = "DAS", ylab = 'ta');
  # lines(weeddata[['x']], as.vector(predict(gr)), col="red",lty=2,lwd=3);
  # text(170, min(log(sub[[trait]]))+1.5, paste(ename, '\n', 'phi1: ', round(alpha[1],2), 
  #                                             '\n  phi2:', round(alpha[2],2), '\n',  ' phi3:', round(alpha[3],2), sep=''),cex = 0.8);
  # 
  return(i);
}
# Equation 1
eqFit = function(xp,yp)
{
  r = list()
  i = getlogistic(xp, yp);
  tdat = xp[-i];
  ydat = yp[-i];
  weeddata = data.frame(y=ydat, x=tdat);
  r$o1 = eqFitcl1(tdat, ydat, weeddata);
  r$o2 = eqFitcl2(tdat, ydat, weeddata);
  
  
  k = min(ydat) + (min(ydat)/3);
  plot(tdat, ydat);
  lines(tdat, fitted(r$o1), col = 'red', lwd = 2);
  lines(tdat, fitted(r$o2), col = 'blue', lwd = 2);
  legend(290,k, legend = c('eq1', 'eq2'), col=c('red', 'blue'), lty=1);
  r$ydat = ydat;
  return(r)
}

st = function(x,y, opt)
{
  if(opt == '1')
  {
    mx = max(y); mn=min(y); mm=0;
    start1= c(lo = mm, Pgmax = mx, RD=mn);
  }
  if(opt == '2')
  {  
    mx = max(y); mn=min(y); mm=5;
    start1= c(lo = mm, Pgmax = mx, RD=mn);
  }
  return(start1)
}

makeFunctions = function(ydat, tdat)
{
  
  elist = list();
  elist[['1']] = c(eq = ydat ~ ((lo*tdat*Pgmax) / (lo*tdat+Pgmax)) - RD, 
                   f=rhs1, qy = "(lo*(Pgmax^2)) / ((lo*tdat+Pgmax)^2)");
  
  elist[['2']] = c(eq = ydat ~ ((tdat * Pgmax)/(tdat + lo)) - RD, 
                   f=rhs2, qy = "(lo*Pgmax) / ((tdat+lo)^2)");
  elist[['10']] = c(eq = ydat ~ P25*exp((c-ha)/(8.314*tdat)), 
                   f=rhsfit1, qy = "");
  return(elist);
}

# Calculate growth parameters
nlsLMFunction = function(x,y, tname, spc, tempr, repl, fn)
{
  b = c();
  tdat = x;
  ydat = y;
  elist = makeFunctions(ydat, tdat);
  weeddata = data.frame(y=ydat, x=tdat);
  
  for(n in fn)
  {
    optx = optima(st(tdat, ydat, n), weeddata, elist[[n]]$f);
    
    o = nlsLM(elist[[n]]$eq, 
                start=list(lo = optx$lo[1], Pgmax = optx$Pgmax[1], RD=optx$RD[1]), trace = F);
    lo = optx$lo[1]; Pgmax = optx$Pgmax[1]; RD=optx$RD[1];
    b = rbind(b, c(eq = n, rep = repl, lo = lo, Pgmax=Pgmax, RD=RD, 
                   specie = spc, temperature = tempr));
    plotCurve(tdat, ydat, fitted(o), n, 'NET PHOTOSYNTHESIS RATE',spc, tempr, repl);
    plotCurve(tdat, eval(parse(text=elist[[n]]$qy)), eval(parse(text=elist[[n]]$qy)), n, 
              'QUANTUM YIELD', spc, tempr, repl);
  
  }
  return(b)
}

# calculate temperature
getPar = function(data, varn)
{
  data = data.frame(data)
  data = data.frame(data, tk = sapply(data[[varn]], function(x) x+273.15));
  return(data);
}

convert2factor = function(data, variable_list)
{
  
  copydata = data;
  for(v in variable_list)
  {
    copydata[[v]]  = factor(copydata[[v]])
  }
  return(copydata)
}

r = read.table('merge_data.csv', header=TRUE, sep=',');
r$I = as.numeric(as.character(r$I));
r = r[complete.cases(r),];
row.names(r) = 1:nrow(r);
b = c();



#pdf('curves.pdf');
#par(mfrow=c(2,2));
for(sp in unique(r$specie)[1])
{
  sub = r[which(r$specie == sp),];
  for(rep in unique(sub$r))
  {
    sub1 = sub[which(sub$r == rep),];
    for(temp in unique(sub1$t))
    {
      sub2 = sub1[which(sub1$t == temp),];
      tname= paste(unique(as.character(sub2$specie)), 'Temp:', unique(sub2$t), 'rep: ', unique(sub2$r), sep=' ');
      b = rbind(b, nlsLMFunction(sub2$I, sub2$Pn, tname, unique(as.character(sub2$specie)), 
                                 unique(sub2$t), unique(sub2$r), c('1') ));
    }
  }
}
break();
#dev.off()

b = convert2numeric(data.frame(b), c("eq","rep","lo","Pgmax","RD","temperature"));
b = getPar(b, 'temperature');
b = data.frame(b);
b = b[order(b$tk),];
b$RD = -1*(b$RD);
# Pgmax
x = b$tk; y = b$Pgmax;
m = eqFit(x,y);
# lo
x = b$tk; y = b$lo;
m = eqFit(x,y);
# RD
x = b$tk; y = -1*(b$RD);
m = eqFit(x,y);

fname = paste('gr_param', '.pdf', sep='');
pdf(fname);
for(tname in c('lo', 'Pgmax', 'RD'))
{
  print(xyplot(round(as.numeric(b[[tname]])) ~ 
                 b[['temperature']] | b[['specie']], data=b, ylab=tname, xlab = 'temperature'));
}
dev.off()

ggplot(b, aes(x=temperature, y = round(RD) )) + 
       geom_point() + facet_wrap("area") 
       