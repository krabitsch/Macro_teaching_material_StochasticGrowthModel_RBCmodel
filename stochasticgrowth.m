% stochasticgrowth.m
%
% Solving the log-linearized stochastic growth model with:
%
% a) method of undetermined coefficients
% b) matrix algebra by solving a simple eigenvalue problem (based on diagonalization)
% c) matrix algebra by solving a generalized eigenvalue problem (based on triagonalization, generalized Schur/ QZ decomposition)
%
% and using model solution to study impulse responses and simulation of the stochastic growth model
%
% Foundations of Macroeconomics, WU Vienna, Nov. 2023
%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 format compact;
 clc; clear all; close all;
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Growth model without population growth or technological progress
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 %      max E(0) SUM(t=0 to inf) {betta^(t)*u(C(t))}
 %
 % s.t. K(t) + K(t+1) = A(t)*K(t)^alfa + (1-delta) * K(t), K(0) given
 %
 % Utility function:    u(C)   = C^(1-gam)/(1-gam)
 % Production function: Y(K,A) = A*K^alfa

 % model is nonlinear system of difference equation given by:
 %
 % C(t)^(-gam)   = betta * C(t+1)^(-gam) * [1-delta + alfa*A(t+1)*K(t+1)^(alfa-1)]
 % C(t) + K(t+1) = A(t)*K(t)^alfa + (1-delta) * K(t)
 % ln[A(t+1)]    = rhoA * ln[A(t)] + eps(t+1)
 %
 % in log-linear form:
 %
 % -gam * c(t) = -gam * c(t+1) + [1-betta*(1-delta)] * [(alfa-1) * k(t+1) + a(t+1)]
 % (css/yss) * c(t) + (kss/yss) * k(t+1) = alfa * k(t) + a(t) + (1-delta) * (kss/yss) * k(t)
 % a(t+1) = rhoA * a(t)
 
 % Parameter values
 alfa   = 0.33;
 betta  = 0.99;
 delta  = 0.025;
 gam    = 2;
 rhoA   = 0.95;
 Ass    = 1; 
 
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

 disp('Model can be written in matrix notation as: A *E(t) z(t+1) = B * z(t), where')
 
 A  = [(1-betta*(1-delta))*(alfa-1)       (1-betta*(1-delta))        -gam  ;
            0                              1                          0    ;
        kss/yss                            0                          0    ]
    
 B  = [       0                       0              -gam       ;
              0                       rhoA            0         ;  
        alfa+(1-delta)*(kss/yss)      1              -(css/yss) ]
 
 W  = inv(A)*B

 
% *************************************************************************
 disp('Computing solution by method of undetermined coefficients')
% *************************************************************************
 % method of undetermined coefficients   
 %
 % denote W = [w11 w12] 
 %            [w21 w22]
 %
 % [k(t+1)]   [w11  w12  w13]   [k(t)]
 % [a(t+1)] = [w21  w22  w23] * [a(t)]
 % [c(t)  ] = [w31  w32  w33] * [c(t)]
 %
 % guess that solution is of form:
 % k(t+1) = phiKK *k(t) + phiKA *a(t)
 % c(t)   = phiCK *k(t) + phiCA *a(t)
 %
 % plug in guesses to find system of 4 equations in unknown coefficients phiK and phiC:
 %
 % phiKK                       = w11 + w13 * phiCK
 % phiKA                       = w12 + w13 * phiCA
 % phiCK * phiKK               = w31 + w33 * phiCK
 % phiCK * phiKA + phiCK * rho = w21 + w22 * phiC

 phiCK = (-(W(1,1)-W(3,3)) - sqrt((W(1,1)-W(3,3))^2+4*W(1,3)*W(3,1)))/(2*W(1,3))  % select stable root
 phiKK =  W(1,1) + W(1,3)*phiCK
 phiCA = (W(3,2) - phiCK*W(1,2))/(phiCK*W(1,3)+rhoA-W(3,3))
 phiKA =  W(1,2) + W(1,3)*phiCA

  

 % use model solution to generate impulse responses or simulate model
 % Let x = [k a] and y = [c y i]
 
 H = [phiKK phiKA;
      0     rhoA]
  
 G = [phiCK phiCA];
 
 G= [G;                                       % policy for c
     alfa                      1;             % policy for y
     (phiKK-(1-delta))/delta  phiKA/delta]    % policy for i
 
%
% *************************************************************************
% disp('Computing solution by simple matrix algebra (based on diagonalizing W)')
% ************************************************************************* 
 % diagonalization of matrix W to find eigenvalues and eigenvectors
 % W = Q * D * inv(Q)

%   % diagonalization of matrix W to find eigenvalues and eigenvectors
%  [Q,D]      = eig(W)
 
 % ... follow similar steps as in 'growth.m' ... 
 
%% Impulse responses
 x(:,1) = [0;.01] % capital (first  variable in state vector x starts at 0 (stst), 
                  % product.(second variable in state vector x starts at 1% above stst,  

 for i=1:40
%  i 
 x(:,i+1) = H*x(:,i);
 end
 %%
 y = (G*x)';
 x = x';
 
 % Plot Impulse Responses
 figure
 plot(x(:,1)*100,'k-s','MarkerSize',3)
 hold on 
 plot(y(:,1)*100,'r-+','MarkerSize',3)
 hold on 
 plot(y(:,2)*100,'b-d','MarkerSize',3)
 hold on 
 plot(y(:,3)*100,'m-^','MarkerSize',3)
 hold on 
 plot(x(:,2)*100,'c-v','MarkerSize',3)
 hold on 
 ylabel('Percent Deviations')
 xlabel('Quarters')
 hold off
 legend('k','c','y','i','a')
 title('Impulse responses to 1% productivity increase, percentage deviations from stst') 
 filename=['IR_stochgrowth.eps'];
 print ('-depsc2','-loose',filename);

 
%% Simulated Time Series
 clear x y
 sigmae = 0.007
 e      = sigmae*randn(100,1);
 x(:,1) = [0 e(1)];
 for i=1:length(e)
 x(:,i+1) = H*x(:,i)+[0;1]*e(i);
 end
 
 y = (G*x)';
 x = x';
 
 %% Plot simulated series
 figure
 plot((y(1:100,1))*100,'r-','MarkerSize',3)
 hold on 
 plot((y(1:100,2))*100,'b-','MarkerSize',3)
 hold on 
 plot((y(1:100,3))*100,'m-','MarkerSize',3)
 ylabel('percent deviations')
 xlabel('Quarters')
 hold off
 legend('c','y','i')
 title('Simulated series of y, c and i, percentage deviations from stst') 
 filename=['SIM_stochgrowth.eps'];
 print ('-depsc2','-loose',filename);



% *************************************************************************
 disp('Solving with Klein algorithm (generalized Schur decomposition)')
% ************************************************************************* 
 % call solab.m (by Paul Klein) to do generalized Schur decomposition on model  A * E(t)[z(t+1)] = B * z(t)
 %
 % solab.m need to be stored in same directory as "growth.m" or be set as path: in command bar: -> 'File' -> 'Set Path'
 
[G,H]=solab(A,B,2)
