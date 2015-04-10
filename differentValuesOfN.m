L=500;L2=1000;L3=1500;
d=round(rand(1,L)); 
d2=round(rand(1,L2));
d3=round(rand(1,L3));% Data sequence
b=2*d-1;b2=2*d2-1;b3=2*d3-1;% Convert unipolar to bipolar
T=1; % Bit duration
Eb=T/2; % This will result in unit amplitude waveforms
fc=3/T; % Carrier frequency
t=linspace(0,5,500);t2=linspace(0,5,1000);t3=linspace(0,5,1500); % discrete time sequence between 0 and 5*T
N=length(t);N2=length(t2);N3=length(t3); % Number of samples
Nsb=N/length(d);Nsb2=N2/length(d2); Nsb3=N3/length(d3);% Number of samples per bit
dd=repmat(d',1,Nsb);dd2=repmat(d2',1,Nsb2);dd3=repmat(d3',1,Nsb3); % replicate each bit Nsb times
bb=repmat(b',1,Nsb);bb2=repmat(b2',1,Nsb2); bb3=repmat(b3',1,Nsb3);
dw=dd';dw2=dd2'; dw3=dd3';% Transpose the rows and columns
dw=dw(:)'; dw2=dw2(:)';dw3=dw3(:)';
% Convert dw to a column vector (colum by column) and convert to a row vector
bw=bb';bw2=bb2';bw3=bb3';
bw=bw(:)'; bw2=bw2(:)';bw3=bw3(:)';% Data sequence samples
w=sqrt(2*Eb/T)*cos(2*pi*fc*t);w2=sqrt(2*Eb/T)*cos(2*pi*fc*t2);w3=sqrt(2*Eb/T)*cos(2*pi*fc*t3); % carrier waveform
bpsk_w=bw.*w; bpsk_w2=bw2.*w2;bpsk_w3=bw3.*w3;% modulated waveform


snr_dB = -15;
% SNR in decibels
snr = 10.^(snr_dB./10);

Pf = 0.01:0.01:0.9; % Pf = Probability of False Alarm
%% Simulation to plot Probability of Detection (Pd) vs. Probability of False Alarm (Pf) 
for m = 1:length(Pf)
    
    i = 0;i2=0;i3=0;
for kk=1:10000 % Number of Monte Carlo Simulations
 n = randn(1,L);n2 = randn(1,L2);n3 = randn(1,L3); %AWGN noise with mean 0 and variance 1
 %s = sqrt(snr).*randn(1,L); % Real valued Gaussina Primary User Signal 
 s = sqrt(snr)*bpsk_w; s2 = sqrt(snr)*bpsk_w2;s3 = sqrt(snr)*bpsk_w3;
 y = s + n;y2 = s2 + n2;y3=s3+n3; % Received signal at SU
 energy = abs(y).^2;  energy2 = abs(y2).^2;energy3 = abs(y3).^2;% Energy of received signal over N samples
 energy_fin =(1/L).*sum(energy);energy_fin2 =(1/L2).*sum(energy2);energy_fin3 =(1/L3).*sum(energy3); % Test Statistic for the energy detection
 thresh(m) = (qfuncinv(Pf(m))./sqrt(L))+ 1; thresh2(m) = (qfuncinv(Pf(m))./sqrt(L2))+ 1; thresh3(m) = (qfuncinv(Pf(m))./sqrt(L3))+ 1; % Theoretical value of Threshold, refer, Sensing Throughput Tradeoff in Cognitive Radio, Y. C. Liang
 if(energy_fin >= thresh(m))  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) counter by 1
     i = i+1;
 end
  if(energy_fin2 >= thresh2(m))  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) counter by 1
     i2 = i2+1;
  end
  if(energy_fin3 >= thresh3(m))  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) counter by 1
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
title('Fig.2. ROC curves for different N at SNR= -15dB')
ylabel('P_D')
xlabel('P_{FA}')
legend('N=500','N=1000','N=1500','Location','southeast')
hold on

