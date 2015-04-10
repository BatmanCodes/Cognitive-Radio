L=500;
d=round(rand(1,L)); % Data sequence
b=2*d-1; % Convert unipolar to bipolar
T=1; % Bit duration
Eb=T/2; % This will result in unit amplitude waveforms
fc=3/T; % Carrier frequency
t=linspace(0,5,500); % discrete time sequence between 0 and 5*T (500 samples)
N=length(t); % Number of samples
Nsb=N/length(d); % Number of samples per bit
dd=repmat(d',1,Nsb); % replicate each bit Nsb times
bb=repmat(b',1,Nsb); dw=dd'; % Transpose the rows and columns
dw=dw(:)'; 
% Convert dw to a column vector (colum by column) and convert to a row vector
bw=bb';
bw=bw(:)'; % Data sequence samples
w=sqrt(2*Eb/T)*cos(2*pi*fc*t); % carrier waveform
bpsk_w=bw.*w; % modulated waveform



snr_dB = -10;
snr_dB2 = -15;
snr_dB3 = -20;% SNR in decibels
snr = 10.^(snr_dB./10);
snr2 = 10.^(snr_dB2./10);
snr3 = 10.^(snr_dB3./10);% Linear Value of SNR
Pf = 0.01:0.01:1; % Pf = Probability of False Alarm
%% Simulation to plot Probability of Detection (Pd) vs. Probability of False Alarm (Pf) 
for m = 1:length(Pf)
    
    i = 0;i2=0;i3=0;
for kk=1:10000 % Number of Monte Carlo Simulations
 n = randn(1,L); %AWGN noise with mean 0 and variance 1
 %s = sqrt(snr).*randn(1,L); % Real valued Gaussina Primary User Signal 
 s = sqrt(snr)*bpsk_w; s2 = sqrt(snr2)*bpsk_w;s3 = sqrt(snr3)*bpsk_w;
 y = s + n;y2 = s2 + n;y3=s3+n; % Received signal at SU
 energy = abs(y).^2;  energy2 = abs(y2).^2;energy3 = abs(y3).^2;% Energy of received signal over N samples
 energy_fin =(1/L).*sum(energy);energy_fin2 =(1/L).*sum(energy2);energy_fin3 =(1/L).*sum(energy3); % Test Statistic for the energy detection
 thresh(m) = (qfuncinv(Pf(m))./sqrt(L/2))+ 1; % Theoretical value of Threshold, refer, Sensing Throughput Tradeoff in Cognitive Radio, Y. C. Liang
 if(energy_fin >= thresh(m))  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) counter by 1
     i = i+1;
 end
  if(energy_fin2 >= thresh(m))  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) counter by 1
     i2 = i2+1;
  end
  if(energy_fin3 >= thresh(m))  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) counter by 1
     i3 = i3+1;
 end
end
Pd(m) = i/kk; 
Pd2(m) = i2/kk;
Pd3(m) = i3/kk;
end

figure 
plot(Pf, Pd,Pf,Pd2,Pf,Pd3)

hold on

thresh = (qfuncinv(Pf)./sqrt(L))+ 1;
Pd_the = qfunc(((thresh - (snr + 1)).*sqrt(L))./(sqrt(2).*(snr + 1)));
plot(Pf, Pd_the,'g')
title('Fig.1 ROC curves for different SNR at N=500')
ylabel('P_D')
xlabel('P_{FA}')
legend('SNR=-10dB','SNR=-15dB','SNR=-20dB','theoretical plot for SNR=-10dB','Location','southeast')
hold on

