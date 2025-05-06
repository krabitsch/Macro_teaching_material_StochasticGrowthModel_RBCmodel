// standard RBC model with log preferences, decentralized economy, gives solution in terms of level variables (absolute deviations from stst)
// used in some of my courses at WU (e.g., Foundations of Macro, Macro Models and Methods, Advanced Macro II)
// Katrin Rabitsch
// last update fall term 2024

var C K N A W RK Y INVE;

varexo eA;

parameters theta betta delta alfa rhoA sigA A_ss KN_ss CN_ss N_ss K_ss C_ss W_ss RK_ss Y_ss INVE_ss;

theta = 3.48;
betta = 0.99;
alfa  = 1/3;
delta = 0.025;
rhoA  = 0.97;
sigA  = 0.007;

A_ss  = 1;
KN_ss = ( (1/betta-(1-delta))/(alfa*A_ss) )^(1/(alfa-1));
CN_ss = KN_ss^(alfa) - delta*KN_ss;
N_ss  = ( ( 1 / ((1-alfa)*KN_ss^alfa) ) * theta * CN_ss + 1 )^(-1);
K_ss  = KN_ss * N_ss;
C_ss  = CN_ss * N_ss;
RK_ss  = alfa*A_ss*K_ss^(alfa-1)*N_ss^(1-alfa);
W_ss = (1-alfa)*A_ss*K_ss^(alfa)*N_ss^(-alfa);
Y_ss    = A_ss * K_ss^alfa * N_ss^(1-alfa);
INVE_ss = K_ss - (1-delta)*K_ss;

model;

// labor-leisure choice, labor supply
theta/(1-N) * C = W;

// Euler equation, capital supply
1/C = betta * (1/C(+1)) * ( RK(+1) + 1 - delta );

// resource constraint
C + INVE = Y;

// firms problem: labor demand
RK = alfa*A*K(-1)^(alfa-1)*N^(1-alfa);

// firms problem: capital demand
W  = (1-alfa)*A*K(-1)^(alfa)*N^(-alfa);

// exogenous AR(1) for TFP:
log(A) = rhoA*log(A(-1)) + eA;

// definitions, output:
Y = A * K(-1)^alfa * N^(1-alfa);

// definitions, investment:
INVE = K - (1-delta)*K(-1);

end;


initval;
C  = C_ss;
K  = K_ss;
A  = A_ss;
N  = N_ss;
W  = W_ss;
RK = RK_ss;
Y  = Y_ss;
INVE = INVE_ss;
end;

resid;

steady;


shocks;
var eA = sigA^2;
//var eA = 1;
end;


stoch_simul(order=1,irf=20);
