%
% Function: solab
%
% Purpose: Solves for the recursive representation of the stable solution to a system
% of linear difference equations.
%
% Inputs: Two square matrices a and b and a natural number nk
%
% a and b are the coefficient matrices of the difference equation
%
% a*x(t+1) = b*x(t)
% 
% where x(t) is arranged so that the state variables come first, and
%
% nk is the number of state variables.
%
% Outputs: the decision rule f and the law of motion p. If we write
%
% x(t) = [k(t);u(t)] where k(t) contains precisely the state variables, then
% 
% u(t)   = f*k(t) and
%
% k(t+1) = p*k(t).
%
% Calls: qzdiv
%

function [f,p] = solab(a,b,nk);

[s,t,q,z] = qz(a,b);            % upper triangular factorization of the matrix pencil b-za
[s,t,q,z] = qzdiv(1,s,t,q,z);   % reordering of generalized eigenvalues in ascending order

z21 = z(nk+1:end,1:nk);
z11 = z(1:nk,1:nk);
rango=rank(z11);
if rank(z11)<nk;
	error('Invertibility condition violated')
end

z11i = z11\eye(nk);
s11 = s(1:nk,1:nk);
t11 = t(1:nk,1:nk);
malo1=t(nk,nk);
malo2=(s(nk,nk));
malisimo1=t(nk+1,nk+1);
malisimo2=(s(nk+1,nk+1));
if abs(t(nk,nk))>abs(s(nk,nk)) | abs(t(nk+1,nk+1))<abs(s(nk+1,nk+1));
   warning('Wrong number of stable eigenvalues.');
end

dyn = s11\t11;

f = real(z21*z11i);
p = real(z11*dyn*z11i);