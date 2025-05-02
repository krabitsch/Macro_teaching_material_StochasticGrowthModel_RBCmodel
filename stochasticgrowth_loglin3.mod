// StochasticGrowth_loglinear2.mod
//
// Solving the stochastic growth model with Dynare, does Taylor expansion in logs
//
// Foundations of Macroeconomics, WU Vienna, Nov. 2023


var c k a y inve;
varexo e;
parameters alfa betta delta gam rhoA Ass css yss kss iss sigma;

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

model;
  css^(-gam)*exp(c)^(-gam)   = betta * css^(-gam)*exp(c(+1))^(-gam) * (1-delta + alfa*Ass*exp(a(+1))*kss^(alfa-1)*exp(k)^(alfa-1));
  css*exp(c) + kss*exp(k) = exp(a)*kss^alfa*exp(k(-1))^alfa + (1-delta) * kss*exp(k(-1));
  yss*exp(y)          = exp(a)*kss^alfa*exp(k(-1))^alfa;
  kss*exp(k)          = iss*exp(inve) + (1-delta)*kss*exp(k(-1));
  (a)             = rhoA * (a(-1)) + e;
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

stoch_simul(periods=2100,order=1);


