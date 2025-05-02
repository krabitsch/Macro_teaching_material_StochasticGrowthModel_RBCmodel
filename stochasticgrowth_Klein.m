% Solving the log-linearized stochastic growth model with the Klein algorithm:
%
% uses solab.m 
%
% Foundations of Macroeconomics, WU Vienna, Nov. 2023
%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 format compact;
 clc; clear all; close all;
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Growth model without population growth or technological progress
 % 
 % model in log-linear form:
 %
 % -gam * c(t) = -gam * c(t+1) + [1-betta*(1-delta)] * [(alfa-1) * k(t+1) + a(t+1)]
 % (css/yss) * c(t) + (kss/yss) * k(t+1) = alfa * k(t) + a(t) + (1-delta) * (kss/yss) * k(t)
 % a(t+1) = rhoA * a(t)
 
 % Parameter values
 alfa   = 0.33
 betta  = 0.99
 delta  = 0.025
 gam    = 2
 rhoA   = 0.95
 Ass    = 1
 
 % Steady state capital, consumption and investment:
  disp('steady state values:')
 kss    = ((alfa*Ass)/(((1-betta)/betta)+delta))^(1/(1-alfa))
 css    = Ass*kss^alfa - delta*kss
 iss    = delta*kss
 yss    = Ass*kss^alfa 
 
 % Write log-linearized model as:
 %
 % A * E(t)z(t+1) = B * z(t)
 %
 % construct matrices A and B:
 A  = zeros(3,3); 
 B  = zeros(3,3);
 
 ik = 1; % column index for capital      in matrix A and B
 ia = 2; % column index for productivity in matrix A and B
 ic = 3; % column index for consumption  in matrix A and B
 
 % equation 1: -gam * c(t) = -gam * c(t+1) + [1-betta*(1-delta)] * [(alfa-1) * k(t+1) + a(t+1)]
 A(1,ik) = (1-betta*(1-delta))*(alfa-1);
 A(1,ia) = (1-betta*(1-delta));
 A(1,ic) = -gam;
 B(1,ic) = -gam;
 
 % equation 2: (css/yss) * c(t) + (kss/yss) * k(t+1) = alfa * k(t) + a(t) + (1-delta) * (kss/yss) * k(t)
 A(2,ik) =  kss/yss;
 B(2,ik) =  alfa+(1-delta)*(kss/yss);
 B(2,ia) =  1;
 B(2,ic) = -(css/yss);

 % equation 3: a(t+1) = rhoA * a(t)
 A(3,ia) =  1;
 B(3,ia) = rhoA;
 

% call solab.m (by Paul Klein) to do generalized Schur decomposition on model  A * E(t)[z(t+1)] = B * z(t)
%
% solab.m need to be stored in same directory as "growth.m" or be set as path: in command bar: -> 'File' -> 'Set Path'
 
[G,H]=solab(A,B,2)



%%
% adding definitions: -> vector of controls now contains y(t) = [c(t); y(t); i(t)] 
% -> matrix A is no longer invertible

A2 = zeros(5,5);
B2 = zeros(5,5);

A2(1:3,1:3) = A;
B2(1:3,1:3) = B;

iy = 4;  % column index for output     in matrix A and B
ii = 5;  % column index for investment in matrix A and B

% equation 4: y(t) = a(t) + alfa * k(t)
B2(4,iy) = -1;
B2(4,ia) =  1;
B2(4,ik) = alfa;

% equation 5: i(t) = (1/delta) * k(t+1) - ((1-delta)/delta) * k(t)
B2(5,ii) =  1;
A2(5,ik) = (1/delta);
B2(5,ik) = ((1-delta)/delta);

% call solab.m (by Paul Klein) to do generalized Schur decomposition on model  A * E(t)[z(t+1)] = B * z(t)
%
% solab.m need to be stored in same directory as "growth.m" or be set as path: in command bar: -> 'File' -> 'Set Path'
 
[G,H]=solab(A2,B2,2)

