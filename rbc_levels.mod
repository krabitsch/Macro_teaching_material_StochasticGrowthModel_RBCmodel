// standard RBC model with log preferences, social planner problem, gives solution in terms of level variables (absolute deviations from stst)
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

theta/(1-N) * C = (1-alfa)*A*(K(-1)/N)^alfa;

C^(-1) = betta * C(+1)^(-1) * ( alfa*A(+1)*K^(alfa-1)*N(+1)^(1-alfa) + 1 - delta );

C + K - (1-delta)*K(-1) = A*K(-1)^alfa*N^(1-alfa);

log(A) = rhoA*log(A(-1)) + eA;

end;

initval;
C  = C_ss;
K  = K_ss;
A  = A_ss;
N  = N_ss;
end;

shocks;
var eA = sigA^2;
end;

stoch_simul(order=1, irf=20);

