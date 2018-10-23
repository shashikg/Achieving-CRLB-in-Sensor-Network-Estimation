% System Formulation, Data Generation --------------------
clear
clc

Niter = 1000;
At = zeros(100, Niter);
Bt = zeros(100, Niter);
A_CRLB = zeros(100, Niter);
A_MSEi = zeros(100, Niter);
A_MSEr = zeros(100, Niter);
B_CRLB = zeros(100, Niter);
B_MSEi = zeros(100, Niter);
B_MSEr = zeros(100, Niter);
alpha_CRLB = zeros(1, Niter);
alpha_MSEi = zeros(1, Niter);
alpha_MSEr = zeros(1, Niter);
M = 100;
N = 100;
n = -50:1:50;
nsq = n.*n;
alpha = 1.0;
rm = -1 + 2.*rand(100,1);
B = -1 + 2.*rand(100,1);
A = 0.1.*(10).^(rm);
k = randn(100,1);

for iter = 1:Niter

  V = 0.1.*randn(100, 101);
  one = ones(1,101);
  X = alpha.*(A*nsq) + (A.*B)*n + A*one + V;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% Initial Estimates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Constant Calc -----------------------------------
  N2 = sum(nsq);
  N4 = sum(nsq.*nsq);
  c1 = 20*(N+1)*(N+2);
  c2 = (3*N*N + 6*N -4)*(N+1)*(N+2);
  A_CRLB(:,iter) = 9.*(A.^2).*0.01./(4.*N.*sum(A.^2));
  alpha_CRLB(:,iter)   =   (9*(alpha.^2)*0.01)/(4*N*sum(A.^2));
  B_CRLB(:,iter)      =   ((9*(B.^2)*0.01)./(4*N*sum(A.^2))) + ((12*0.01)./((A.^2)*N.^3));

  %Estimation Start ----------------------------------
  Xn = sum(X,1);
  s1 = sum(Xn);
  s2 = sum(Xn.*nsq);
  alpha_hat = (c1.*s1 - 240.*s2)./(c1.*s2 - c2.*s1)

  %Other Est
  X0 = sum(X,2);
  X1 = sum(X.*n,2);
  X2 = sum(X.*nsq,2);
  Bhat = (X1./N2).*((alpha_hat^2.*N4 + 2.*alpha_hat.*N2 + N)./(alpha_hat.*X2 + X0));
  Ahat = (alpha_hat.*X2 + X0)./(alpha_hat^2.*N4 + 2.*alpha_hat.*N2 + N);

  %MSE Calculation-------------------------------
  Ae = A - Ahat;
  A_MSEi(:,iter) = 10.*log10(Ae.*Ae);
  Be = B - Bhat;
  B_MSEi(:,iter) = 10.*log10(Be.*Be);
  alphae = alpha - alpha_hat;
  alpha_MSEi(:,iter) = 10.*log10(alphae.*alphae);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for i=1:30
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Re Estimates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %alpha estimate ---------------------------------
    Xn = sum(Ahat.*X,1);
    s1 = sum(Xn);
    s2 = sum(Xn.*nsq);
    alpha_hat = (c1.*s1 - 240.*s2)./(c1.*s2 - c2.*s1);

    %Other Est
    X0 = sum(X,2);
    X1 = sum(X.*n,2);
    X2 = sum(X.*nsq,2);
    Bhat = (X1./N2).*((alpha_hat^2.*N4 + 2.*alpha_hat.*N2 + N)./(alpha_hat.*X2 + X0));
    Ahat = (alpha_hat.*X2 + X0)./(alpha_hat^2.*N4 + 2.*alpha_hat.*N2 + N);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end 
  %MSE Calculation-------------------------------
  Ae = A - Ahat;
  A_MSEr(:,iter) = 10.*log10(Ae.*Ae);
  Be = B - Bhat;
  B_MSEr(:,iter) = 10.*log10(Be.*Be);
  alphae = alpha - alpha_hat;
  alpha_MSEr(:,iter) = 10.*log10(alphae.*alphae);
  
  At(:,iter) = A;
  Bt(:,iter) = B;
end
%alpha_hat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  PLOTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
plot(abs(mean(At, axis=2)), 10.*log10(mean(alpha_CRLB, axis=2)));
hold on;
plot(abs(mean(At, axis=2)), mean(alpha_MSEi, axis=2));
hold on;
plot(abs(mean(At, axis=2)), mean(alpha_MSEr, axis=2));
hold on;

figure(2);
plot(abs(mean(At, axis=2)), 10.*log10(mean(A_CRLB, axis=2)), '*');
hold on;
plot(abs(mean(At, axis=2)), mean(A_MSEi, axis=2), 'o');
hold on;
plot(abs(mean(At, axis=2)), mean(A_MSEr, axis=2), '+');
hold on;

figure(3);
plot(abs(mean(Bt, axis=2)), 10.*log10(mean(B_CRLB, axis=2)), '*');
hold on;
plot(abs(mean(Bt, axis=2)), mean(B_MSEi, axis=2), 'o');
hold on;
plot(abs(mean(Bt, axis=2)), mean(B_MSEr, axis=2), '+');
hold on;



