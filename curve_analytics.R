library('minpack.lm');
library("nlmrt")
library("optimx")
library(ggplot2)

rhs <- function(b, mydata) {
  sum((mydata$y - ((b[1]*mydata$x*b[2]) / (b[1]*mydata$x+b[2])) - b[3])^2)
  
}

optima = function(start1, w, coeflist)
{
  ## Optim
  w.par <- coeflist$coefficients;
  w.optx <- optimx(par=start1, fn=rhs, mydata=w, 
                   control=list(all.methods=TRUE, save.failures=TRUE, maxit=2500))
  return(w.optx)
}
  
plotCurve = function(x,y, fm, tname)
{
  plot(x, y, main = tname)
  lines(x, fitted(fm), col = 2, lwd = 2);
}


nlsLMFunction = function(x,y, tname)
{
  mx = max(y); mn=min(y); mm=0;
  start1= c(b1 = mm, b2 = mx, b3=mn);
  tdat = x;
  ydat = y;
  weeddata = data.frame(y=ydat, tt=tdat)
  eunsc = y ~ ((b1*tt*b2) / (b1*tt+b2)) - b3;
  anlxb1g =try(nlxb(eunsc, start=start1, trace=FALSE, data=data.frame(y=ydat, tt=tdat)))
  anlxb1g$coefficients
  
  weeds <- data.frame(y=ydat, x=weeddata$tt)
  weed.par <- anlxb1g$coefficients
  w.optx = optima(start1, weeds,anlxb1g)

  o = nlsLM(ydat ~ ((lo*tdat*Pgmax) / (lo*tdat+Pgmax)) - RD, 
  start=list(lo = w.optx$b1[1], Pgmax = w.optx$b2[1], RD=w.optx$b3[1]), trace = TRUE);
  
  ## plot fitted values
  plotCurve(tdat, ydat, o, tname);
  return(w.optx)
  #print(cor(ydat,fitted(o)));
  
  
}



r = read.table('merge_data.csv', header=TRUE, sep=',');
#break();
b = c();
pdf('curves.pdf');
par(mfrow=c(4,2));
for(cl in unique(r$m)[c(1:3, 6:9)])
{
  for(region in unique(r$location))
  {
    for(temp in c(15,25, 30, 35))
    {
      for(rept in c(1:3))
      {
        i = intersect(intersect(intersect(which(r$r == rept), which(r$t == temp)), 
        which(r$location == region)), 
        which(r$m == cl));
        tna= paste(region, 'Temp: ', temp,'rep: ',rept, 'cl:', cl, sep=' ');
        print(tna)
        o = nlsLMFunction(r$n[i], r$v[i], tna);
        b = rbind(b, c(lo = o$b1[1], Pgmax=o$b2[1], RD=o$b3[1], 
                       area = region, temperature = temp))
        
      }
    }
  }
}
dev.off()

b= data.frame(b)
fname = paste('gr_param', '.pdf', sep='');
pdf(fname);
for(tname in c('lo', 'Pgmax', 'RD'))
{
  print(xyplot(round(as.numeric(b[[tname]])) ~ 
                 b[['temperature']] | b[['area']], data=b, ylab=tname, xlab = 'temperature'));
}
dev.off()

ggplot(b, aes(x=temperature, y = round(RD) )) + 
       geom_point() + facet_wrap("area") 
       