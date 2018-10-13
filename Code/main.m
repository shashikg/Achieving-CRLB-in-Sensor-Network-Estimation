% System Formulation, Data Generation --------------------
clear
clc
M = 100;
N = 100;
n = -50:1:50;
nsq = n.*n;
alpha = 1.0;
rm = -1 + 2.*rand(100,1);
B = -1 + 2.*rand(100,1);
A = 0.1.*(10).^(rm);
V = 0.1.*randn(100, 101);
one = ones(1,101);
X = alpha.*(A*nsq) + (A.*B)*n + A*one + V;

% Constant Calc -----------------------------------
N2 = sum(nsq);
N4 = sum(nsq.*nsq);
c1 = 20*(N+1)*(N+2);
c2 = (3*N*N + 6*N -4)*(N+1)*(N+2);

%Estimation Start ----------------------------------
Ahat = ones(100,1);
Bhat = ones(100,1);
for iter=1:3
  %alpha estimate ---------------------------------
  Xn = sum(Ahat.*X,1);
  s1 = sum(Xn);
  s2 = sum(Xn.*nsq);
  alpha_hat = (c1.*s1 - 240.*s2)./(c1.*s2 - c2.*s1)
  
  %Other Est
  X0 = sum(X,2);
  X1 = sum(X.*n,2);
  X2 = sum(X.*nsq,2);
  Bhat = (X1./N2).*((alpha_hat^2.*N4 + 2.*alpha_hat.*N2 + N)./(alpha_hat.*X2 + X0));
  Ahat = (alpha_hat.*X2 + X0)./(alpha_hat^2.*N4 + 2.*alpha_hat.*N2 + N);
end



