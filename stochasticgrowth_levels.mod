// StochasticGrowth_levels.mod
//
// Solving the stochastic growth model with Dynare, does Taylor expansion in levels
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

model;
  c^(-gam) = betta * c(+1)^(-gam) * (1-delta + alfa*a(+1)*k^(alfa-1));
  c + k    = a*k(-1)^alfa + (1-delta) * k(-1);
      log(a)    = rhoA * log(a(-1))  + e;
      y    = a*k(-1)^alfa;
      k    = inve + (1-delta)*k(-1);
end;


initval;
  k    = kss;
  c    = css;
  a    = Ass;
  y    = yss;
  inve = iss;
end;

steady;

check;

shocks;
var e = sigma^2;
end;

stoch_simul(periods=2100,irf =30,order=1);


