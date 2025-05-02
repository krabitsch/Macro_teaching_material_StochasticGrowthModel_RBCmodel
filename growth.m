% growth.m
%
% Solving the log-linearized growth model with:
%
% a) method of undetermined coefficients
% b) matrix algebra by solving a simple eigenvalue problem (based on diagonalization)
%
% and using model solution to study transitional dynamics of the simple growth model
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
 %      max SUM(t=0 to inf) {betta^(t)*u(C(t))}
 %
 % s.t. K(t) + K(t+1) = K(t)^alfa + (1-delta) * K(t), K(0) given
 %
 % Utility function:    u(C) = C^(1-gam)/(1-gam)
 % Production function: Y(K) = K^alfa

 % model is nonlinear system of difference equation given by:
 %
 % C(t)^(-gam)   = betta * C(t+1)^(-gam) * [1-delta + alfa*K(t+1)^(alfa-1)]
 % C(t) + K(t+1) = K(t)^alfa + (1-delta) * K(t)
 %
 % in log-linear form:
 %
 % -gam * c(t) = -gam * c(t+1) + [1-betta*(1-delta)] * (alfa-1) * k(t+1)
 % (css/yss) * c(t) + (kss/yss) * k(t+1) = alfa * k(t) + (1-delta) * (kss/yss) * k(t)
 
 % Parameter values
 alfa   = 0.33;
 betta  = 0.99;
 delta  = .025;
 gam    = 2; 
 

 
 % Steady state capital, consumption and investment:
 disp('steady state values:')
 kss    = (alfa/(((1-betta)/betta)+delta))^(1/(1-alfa))
 css    = kss^alfa - delta*kss
 iss    = delta*kss
 yss    = kss^alfa 
 
%%
 % Write log-linearized model as:
 %
 % A * z(t+1) = B * z(t)
 %
 % construct matrices A and B:
 disp('Model can be written in matrix notation as: A * z(t+1) = B * z(t), where')
 
 A  = [(1-betta*(1-delta))*(alfa-1)       -gam  ;
        kss/yss                            0    ]
    
 B  = [       0                           -gam  ;
        alfa+(1-delta)*(kss/yss)      -(css/yss)]
 
% obtain model solution:    
%
% matrix A here is invertible: ->

 W  = inv(A)*B

 % z(t+1) = W * z(t)
 %%
% *************************************************************************
 disp('Computing solution by method of undetermined coefficients')
% *************************************************************************
 % method of undetermined coefficients   
 %
 % denote W = [w11 w12] 
 %            [w21 w22]
 %
 % [k(t+1)]   [w11  w12]   [k(t)]
 % [c(t)  ] = [w21  w22] * [c(t)]
 %
 % guess that solution is of form:
 % k(t+1) = phiK *k(t)
 % c(t)   = phiC *k(t)
 %
 % plug in guesses to find system of 2 equations in unknown coefficients phiK and phiC:
 %
 % [phiK       ]   [    ]   [w11 + w12 * phiC] * [    ] 
 % [phiC * phiK] * [k(t)] = [w21 + w22 * phiC] * [k(t)] 
 %
 % ->
 %
 % phiK        = w11 + w12 * phiC
 % phiC * phiK = w21 + w22 * phiC
 %
 % reduce system to a quadratic equation in phiK:
 %
 % phiK^2 - (w11+w22) * phiK + (w11*w22-w12*w21)
 % 
 % solve according to formula: 
 % phiK = [-b +- sqrt( b^2 - 4*a*c)] / 2a 
 %
 
 a =   1;
 b = -(W(1,1)+W(2,2));
 c =   W(1,1)*W(2,2)-W(1,2)*W(2,1);
 
 phiK = (-b - sqrt( b^2 - 4*a*c )) / (2*a)  % select stable root
 phiC = -(W(1,1)/W(1,2)) + phiK/W(1,2)

 
 %%
% *************************************************************************
 disp('Computing solution by simple matrix algebra (based on diagonalizing W)')
% ************************************************************************* 
 % diagonalization of matrix W to find eigenvalues and eigenvectors
 % W = Q * D * inv(Q)
 %
 % can write solution as: 
 %   
 % [k(t) c(t)]' = Q* D^(t) *inv(Q) * [k(0) c(0)]'
 % 
 % denote elements of Q and D by: Q = [q11 q12; q21 q22], D=[d1 0; 0 d2]
 %
 % -> multiplying out gives:
 % 
 % k(t) =   [q11* d1^(t) *(q22*k(0)-q12*c(0)) + q12* d2^(t) *(-q21*k(0)+q11*c(0))] * (1/det(Q)) 
 % c(t) =   [q21* d1^(t) *(q22*k(0)-q12*c(0)) + q22* d2^(t) *(-q21*k(0)+q11*c(0))] * (1/det(Q)) 
 % 
 % with d1 being the eigenvalue the eigenvalue <1 (stable) and d2 being the
 % eigenvalue >1 (unstable) need to set jump variable c(0) such that
 % unstable dynamics are wiped out. That is:
 %
 % c(0) = (q21/q11) * k(0)
 %
 % then solution can be written as:
 % 
 % k(t+1) = coeffK * k(t) = [(1/det(Q)) * q11 *(q22-q12*(q21/q11)) * d1] * k(t)
 % c(t+1) = coeffC * k(t) = [(1/det(Q)) * q21 *(q22-q12*(q21/q11)) * d1] * k(t)
 % 
 % then:
 %
 % k(t+1) = phiK * k(t) 
 % c(t)   = phiC * k(t) 
 %
 % where phiK = coeffK, and phiC = (coeffC/coeffK)
 
  [Q,D]      = eig(W)
 
 if D(1,1)<D(2,2)
  
     coeffK = (1/det(Q))* (Q(1,1) * (Q(2,2)-(Q(1,2)*Q(2,1))/Q(1,1))) * D(1,1);
     coeffC = (1/det(Q))* (Q(2,1) * (Q(2,2)-(Q(1,2)*Q(2,1))/Q(1,1))) * D(1,1);

     phiK = coeffK
     phiC = coeffC/coeffK

 elseif D(2,2)<D(1,1)
     
     coeffK = (1/det(Q))* (Q(1,2) * (-Q(2,1)+(Q(1,1)*Q(2,2))/Q(1,2))) * D(2,2);
     coeffC = (1/det(Q))* (Q(2,2) * (-Q(2,1)+(Q(1,1)*Q(2,2))/Q(1,2))) * D(2,2);

     phiK = coeffK
     phiC = coeffC/coeffK

 end

 
 
%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Compute the transition path when economy starts with a capital stock 10% below its steady state value
 k=zeros(201,1);
 k(1,1) = -.1
 tic
 for i=1:201
 k(i+1,1) = phiK*k(i,1);
 end
 toc
 c = phiC.*k;
 y = alfa.*k;
 i = (1/delta)*phiK*k(2:end,1)-((1-delta)/delta)*k(1:end-1,1);
%%
 % Plot Transition Paths, percent deviations from steady state of variables
 figure
 plot(k*100,'k-s','MarkerSize',3)
 hold on 
 plot(c*100,'r-+','MarkerSize',3)
 hold on 
 plot(y*100,'b-d','MarkerSize',3)
 hold on 
 plot(i*100,'m-^','MarkerSize',3)
 ylabel('Percent Deviations')
 xlabel('Quarters')
 hold off
 legend('k','c','y','i')
 title('Transition dynamics, percentage deviations from stst (no long run growth)') 

 % Plot Transition Paths, levels of variables
 figure
 plot(kss*exp(k),'k-s','MarkerSize',3)
 hold on 
 plot(css*exp(c),'r-+','MarkerSize',3)
 hold on 
 plot(yss*exp(y),'b-d','MarkerSize',3)
 hold on 
 plot(iss*exp(i),'m-^','MarkerSize',3)
 ylabel('Levels')
 xlabel('Quarters')
 hold off
 legend('k','c','y','i')
%  legend('c','y','i')
 title('Transition dynamics, levels (no long run growth)')
 
 