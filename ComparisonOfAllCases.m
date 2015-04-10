L=1500;
snr_dB = -15;

snr = 10.^(snr_dB./10);

Pf = 0.01:0.01:0.9; % Pf = Probability of False Alarm
thresh = (qfuncinv(Pf)./sqrt(L))+ 1;
Pd_the = qfunc(((thresh - (snr + 1)).*sqrt(L))./(sqrt(2).*(snr + 1)));
 
 Pd_the2 = qfunc(((thresh - (snr + (1/1.05))).*sqrt(L))./(sqrt(2).*(snr + (1/1.05))));

 Pd_the3 = qfunc(((thresh/1.03 - (snr + (1/1.05))).*sqrt(L))./(sqrt(2).*(snr + (1/1.05))));
 Pd_the4 = qfunc(((thresh/1.04 - (snr + (1/1.05))).*sqrt(L))./(sqrt(2).*(snr + (1/1.05))));
 figure 
plot(Pf,Pd_the,Pf,Pd_the2,Pf,Pd_the3,Pf,Pd_the4)

hold on


title('Fig.5 ROC curves with different noise uncertainties and different variable thresholds')
ylabel('P_D')
xlabel('P_{FA}')
legend('\rho=1.00,\rho\prime=1.00','\rho=1.05,\rho\prime=1.00','\rho=1.05,\rho\prime=1.03','\rho=1.05,\rho\prime=1.04','Location','southeast')
hold on