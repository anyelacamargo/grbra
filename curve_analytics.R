library('minpack.lm');
library("nlmrt")
library("optimx")
library(ggplot2)

# [64-bit] C:\Program Files\R\R-3.3.1

rhs1 <- function(b, mydata) {
  sum((mydata$y - ((b[1]*mydata$x*b[2]) / (b[1]*mydata$x+b[2])) - b[3])^2)
}


rhs2 <- function(b, mydata) {
  #eunsc = ydat ~ ((tdat * Pgmax)/(tdat + lo)) - RD;
  sum((mydata$y - (((mydata$x*b[2]) / (mydata$x+b[1])) - b[3]))^2)
}

rhsfit <- function(b, mydata) {
  #ydat ~ P25*exp((c-ha)/(8.314*tdat)
  sum((mydata$y - ((b[1]*exp( ((b[2]/(8.314*298.15)) - b[2]) /(8.314*mydata$x)) ))) ^2)
}


optima = function(start1, w, fc)
{
  ## Optim
  w.optx <- optimx(par=start1, fn=fc, mydata=w, 
                   control=list(all.methods=TRUE, save.failures=TRUE, maxit=2500))
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
  data = weeddata;
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
}




# Equation 1
eqFit = function(x,y)
{
  #b1=lo, b2=pgmax, b3=RD
  r = list()
  #mx = max(y); mn=min(y); mm=0;
  tdat = x;
  ydat = y;
  start2= c(P25 = 0.0, ha=20000);
  weeddata = data.frame(y=ydat, x=tdat);
  searchOutlier(weeddata);
  tdat = weeddata$x ;
  ydat = weeddata$y;
  optx = optima(start2, weeddata, rhsfit);
  optx
  eunsc10 = ydat ~ P25*exp((  (ha/(8.314*298.15)) - ha)/(8.314*tdat));
  o = nlsLM(eunsc10, 
            start=list(P25 = optx$P25[1], ha = optx$ha[1]), trace = TRUE);
  o
  plot(tdat, ydat);
  lines(tdat, fitted(o), col = 2, lwd = 2);
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
                   f=rhsfit, qy = "");
  return(elist);
}

# Calculation
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
break();


pdf('curves.pdf');

par(mfrow=c(2,2));
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

b = convert2numeric(data.frame(b), c("eq","rep","lo","Pgmax","RD","temperature"))
b1 = getPar(b, 'temperature');
dev.off()


b= data.frame(b)
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
       