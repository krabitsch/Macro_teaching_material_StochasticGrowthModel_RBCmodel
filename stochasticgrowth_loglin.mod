// StochasticGrowth_loglinear.mod
//
// Solving the stochastic growth model with Dynare, model specified in log-linear terms
//
// Foundations of Macroeconomics, WU Vienna, Nov. 2023


var c k a y inve;
varexo e;
parameters alfa betta delta gam rhoA Ass css yss kss sigma;

 alfa   = 0.33;
 betta  = 0.99;
 delta  = 0.025;
 gam    = 2;
 rhoA   = 0.95;
 Ass    = 1;
 kss    = ((alfa*Ass)/(((1-betta)/betta)+delta))^(1/(1-alfa));
 css    = Ass*kss^alfa - delta*kss;
 iss    = delta*kss;
 yss    = Ass*kss^alfa;
 sigma  = 0.007;

model(linear);
  -gam * c = -gam * c(+1) + (1-betta*(1-delta)) * ((alfa-1) * k + a(+1));
  (css/yss) * c + (kss/yss) * k = alfa * k(-1) + a + (1-delta) * (kss/yss) * k(-1);
  y    = a + alfa * k(-1);
  inve = (1/delta) * k - ((1-delta)/delta) * k(-1); 
  a    = rhoA * a(-1) + e;
end;


initval;
  k    = 0;
  c    = 0;
  a    = 0;
  y    = 0;
  inve = 0;
end;

steady;

check;

shocks;
var e = sigma^2;
end;

stoch_simul(periods=2100);
//stoch_simul(periods=2100,irf=30,order=1) y c inve;



