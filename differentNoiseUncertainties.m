L=1500;

snr_dB = -15;

snr = 10.^(snr_dB./10);

Pf = 0.01:0.01:0.9; % Pf = Probability of False Alarm

 thresh1 = (qfuncinv(Pf)./sqrt(L))+ 1;
 Pd_the1 = qfunc(((thresh1 - (snr + 1)).*sqrt(L))./(sqrt(2).*(snr + 1)));
 thresh2 = (qfuncinv(Pf).*1.01./sqrt(L))+ 1.01;
 Pd_the2 = qfunc(((thresh2 - (snr + (1/1.01))).*sqrt(L))./(sqrt(2).*(snr + (1/1.01))));
 thresh3 = (qfuncinv(Pf).*1.03./sqrt(L))+ 1.03;
 Pd_the3 = qfunc(((thresh3 - (snr + (1/1.03))).*sqrt(L))./(sqrt(2).*(snr + (1/1.03))));
 thresh4 = (qfuncinv(Pf).*1.05./sqrt(L))+ 1.05;
 Pd_the4 = qfunc(((thresh4 - (snr + (1/1.05))).*sqrt(L))./(sqrt(2).*(snr + (1/1.05))));

figure 
plot(Pf,Pd_the1,Pf,Pd_the2,Pf,Pd_the3,Pf,Pd_the4)

hold on


title('Fig.3. ROC curves for different noise uncertainties at SNR=-15dB and N=1500')
ylabel('P_D')
xlabel('P_{FA}')
legend('\rho=1','\rho=1.01','\rho=1.03','\rho=1.05','Location','southeast')
hold on

