function Estimate_alp_Linear_Model_full_network_complete_V1
close all;
clear all;
clc;

format long;

n           =   -50:1:50;
N           =   length(n);

n_sq        =   (n.^2);
n_pow4      =   (n.^4);
% sum_n2      =   (1/12)*N*(N+1)*(N+2);%sum(n_sq);
%
% sum_n4      =   (1/240)*N*(N+1)*(N+2)*(3*N^2 +6*N - 4);%sum(n_pow4);
sum_n2      =   sum(n_sq);
sum_n4      =   sum(n_pow4);

Noise_Pow   =   0.01;
M           =   100;
alp         =   1;
NIter       =   50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Alp         = 1;
rm          = -1:0.00001:1;
A_m_Random  = 0.1*10.^(rm);
B_m_Random  = -1:0.00001:1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Case1 : A is fixed and varying 0.1:0.01:1, and B are random
mat_A = [(n.^2).'   n.'  ones(N,1)] ;
c1 = 20*(N+1)*(N+2);
c2 = (3*N^2 +6*N - 4)*(N+1)*(N+2);
if  1
    h = waitbar(0,'Please wait...');
    for iter = 1:NIter
        waitbar(iter / NIter);
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         aa = 0.1;
%         bb = 1;
%         A_m              =  aa + (bb-aa).*rand(1,M);
%         [am1,I]            =  sort(A_m);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         aa1 = -1;
%         bb1 = 1;
%         B_m              =  aa1 + (bb1-aa1).*rand(1,M);
%         [bm1,I1]         =  sort(B_m);
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ran_index        =  ceil(length(A_m_Random)*rand(1,M));
                A_m              =  A_m_Random(ran_index);
                [am1,I]          =  sort(A_m);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ran_index1       =  ceil(length(B_m_Random)*rand(1,M));
                B_m              =  B_m_Random(ran_index1);
                [bm1,I1]         =  sort(B_m);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Noise               =   sqrt(Noise_Pow)*((randn(M,N) + 1j*randn(M,N))/sqrt(2));
        
        xx_vec              =   [am1.*alp;am1.*bm1;am1];
        x                   =   mat_A*xx_vec + Noise.' ;
        x_sum               =   sum(x.');
        XX_alp              =   ((pinv(mat_A'*mat_A))*mat_A')*x_sum.';
        XX                  =   ((pinv(mat_A'*mat_A))*mat_A')*x;
        
        
        ALP_ntw(iter)       =   XX_alp(1,:)./XX_alp(3,:);
        
        ALP_ntw_sin(iter,:) =    XX(1,:)./XX(3,:);
%         a_est_d(iter,:)     =    real(XX(3,:));
%         Beta_est_d(iter,:)  =    ((XX(2,:))./(XX(3,:)));
        Bmm(iter,:)         =    bm1;
        Amm(iter,:)         =    am1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%  Method in paper
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Computing A_Est for ntwk
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sum_xn              = sum(x);
        sum_xn_n2           = (n_sq*x);
        alp_est             = (ALP_ntw(iter));
        Num_aest            = alp_est*sum_xn_n2 + sum_xn;
        den_aest            = alp_est.^2*sum_n4 + 2*alp_est*sum_n2 + N;
        Am_est_pap(iter,:)  = Num_aest./den_aest;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Computing B_Est
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sum_xn_n              = (n*x);
        Num_best              = alp_est.^2*sum_n4 + 2*alp_est*sum_n2 + N;
        den_best              = alp_est*sum_xn_n2 + sum_xn;
        bm_est_pap(iter,:)    = (Num_best./den_best).*(sum_xn_n./sum_n2);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Computing A_Est for single node
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sum_xn              = sum(x);
        sum_xn_n2           = (n_sq*x);
        
        alp_est1             = real(ALP_ntw_sin(iter,:));
        Num_aest            = alp_est1.*sum_xn_n2 + sum_xn;
        den_aest            = alp_est1.^2.*sum_n4 + 2*alp_est1.*sum_n2 + N;
        Am_est_pap_sin_node(iter,:)  = Num_aest./den_aest;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Computing B_Est
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sum_xn_n              = (n*x);
        Num_best              = alp_est1.^2.*sum_n4 + 2*alp_est1.*sum_n2 + N;
        den_best              = alp_est1.*sum_xn_n2 + sum_xn;
        bm_est_pap_sin_node(iter,:)    = (Num_best./den_best).*(sum_xn_n./sum_n2);
        
        
        
        % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % %             %%%%  Method in paper
        % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sum_xn              =    ((ones(1,N)*(x_sum)'));
        sum_n_xn            =    (n*(x_sum'));
        n_2_xn              =    ((n_sq*(x_sum)'));
        alp_est_d(iter)     =    (c1*sum_xn - 240*n_2_xn) / (c1*n_2_xn - c2*sum_xn);
        % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%  CRLB Computation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        CRLB_ALP(iter,:)      =   (9*(Alp.^2)*Noise_Pow)/(4*N*sum(am1.^2));
        CRLB_Am1(iter,:)      =   (9*(am1.^2)*Noise_Pow)./(4*N*sum(am1.^2));
        CRLB_bm1(iter,:)      =   ((9*(bm1.^2)*Noise_Pow)./(4*N*sum(am1.^2))) + ((12*Noise_Pow)./((am1.^2)*N.^3));
        
    end
    alp_mse                     = mean((real(ALP_ntw)- alp.*ones(1,NIter)).^2);
    alp_sin_mse                 = mean((real(ALP_ntw_sin)- alp.*ones(NIter,M)).^2);
    size(ALP_ntw_sin)
    
    
    %     Amp_MSE  = mean((real(a_est_d)- Amm).^2);
    %     Beta_MSE = mean((real(Beta_est_d)- Bmm).^2);
    
    Amp_MSE_paper  = mean((real(Am_est_pap)- Amm).^2);
    Beta_MSE_paper = mean((real(bm_est_pap)- Bmm).^2);
    
    
    %     Amp_MSE_paper_sin_node  = mean((real(a_est_d)- Amm).^2);
    %     Beta_MSE_paper_sin_node = mean((real(Beta_est_d)- Bmm).^2);
    Amp_MSE_paper_sin_node  = mean((real(Am_est_pap_sin_node)- Amm).^2);
    Beta_MSE_paper_sin_node = mean((real(bm_est_pap_sin_node)- Bmm).^2);
    
    
    % Beta_MSE = mean((real(Beta_est_d)- Bmm).^2);
    % Amp_MSE   = mean((real(a_est_d)- Bmm).^2);
    %mean(Amm)-am1    
    
    
    
    
    close(h);
    
    figure(100);
    plot(am1,((10*log10(alp_sin_mse ))),'b*-');
    hold on;
    plot(mean(Amm),ones(1,M).*10*log10(alp_mse),'r*-');
    hold on;
    plot(am1,ones(1,M).*10*log10(mean(CRLB_ALP)),'c*-');
    xlabel('Amplitude');
    ylabel('MSE(dB)');
    legend('Alpha-Estimation(Single sensor)',...
        'Alpha-Estimation(network sensor)','Alpha-Estimation(CRLB)');
    xlim([0.1 1]);
    
    
    
    figure(200);
    plot(am1,((10*log10(Amp_MSE_paper ))),'b*-');
    hold on;
    plot(am1,10*log10((Amp_MSE_paper_sin_node)),'r*-');
    hold on;
    plot(am1,10*log10(mean(CRLB_Am1)),'c*-');
    xlabel('Amplitude');
    ylabel('MSE(dB)');
    legend('Amplitude-Estimation(network sensor)',...
        'Amplitude-Estimation(Single sensor)','Amplitude-Estimation(CRLB)');
        xlim([0.1 1]);

    
    figure(300);
    plot(abs(bm1),((10*log10(Beta_MSE_paper ))),'b*-');
    hold on;
    plot(abs(bm1),10*log10((Beta_MSE_paper_sin_node)),'r*-');
    hold on;
    plot(abs(bm1),10*log10(mean(CRLB_bm1)),'c*-');
    xlabel('Amplitude');
    ylabel('MSE(dB)');
    legend('Bm-Estimation(network sensor)',...
        'Bm-Estimation(Single sensor)','Bm-Estimation(CRLB)');
        xlim([0.1 1]);

    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
