// standard RBC model with log preferences, social planner problem, gives solution in terms of percentage deviation from stst variables
// used in some of my courses at WU (e.g., Foundations of Macro, Macro Models and Methods, Advanced Macro II)
// Katrin Rabitsch
// last update fall term 2024

var C K N A;

varexo eA;

parameters theta betta delta alfa rhoA sigA A_ss KN_ss CN_ss N_ss K_ss C_ss YN_ss Y_ss;

theta = 3.48;
betta = 0.99;
alfa  = 1/3;
delta = 0.025;
rhoA  = 0.97;
sigA  = 0.007;

A_ss  = 1;
KN_ss = ( (1/betta-(1-delta))/(alfa*A_ss) )^(1/(alfa-1));
YN_ss = A_ss * KN_ss^alfa;
CN_ss = KN_ss^(alfa) - delta*KN_ss;
N_ss  = ( ( 1 / ((1-alfa)*KN_ss^alfa) ) * theta * CN_ss + 1 )^(-1);
K_ss  = KN_ss * N_ss;
C_ss  = CN_ss * N_ss;
Y_ss  = YN_ss * N_ss;

model;

theta/(1-exp(N)) * exp(C) = (1-alfa)*exp(A)*(exp(K(-1))/exp(N))^alfa;

1/exp(C) = betta * 1/exp(C(+1)) * ( alfa*exp(A(+1))*exp(K)^(alfa-1)*exp(N(+1))^(1-alfa) + 1 - delta );

exp(C) + exp(K) - (1-delta)*exp(K(-1)) = exp(A)*exp(K(-1))^alfa*exp(N)^(1-alfa);

(A) = rhoA*(A(-1)) + eA;

end;

initval;
C  = log(C_ss);
K  = log(K_ss);
A  = log(A_ss);
N  = log(N_ss);
end;

shocks;
var eA = sigA^2;
end;

stoch_simul(order=1, irf=20);

