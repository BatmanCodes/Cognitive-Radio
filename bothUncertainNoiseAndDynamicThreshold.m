L=1500;
d=round(rand(1,L)); % Data sequence
b=2*d-1; % Convert unipolar to bipolar
T=1; % Bit duration
Eb=T/2; % This will result in unit amplitude waveforms
fc=3/T; % Carrier frequency
t=linspace(0,5,1500); % discrete time sequence between 0 and 5*T (15000 samples)
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

% plotting commands follow


snr_dB = -15;

snr = 10.^(snr_dB./10);

Pf = 0.01:0.01:0.9; % Pf = Probability of False Alarm
%% Simulation to plot Probability of Detection (Pd) vs. Probability of False Alarm (Pf) 
for m = 1:length(Pf)
    
    i = 0;i2=0;i3=0;i4=0;
for kk=1:10000 % Number of Monte Carlo Simulations
 n = randn(1,L); n2 = sqrt(1.01).*randn(1,L);%AWGN noise with mean 0 and variance 
 %s = sqrt(snr).*randn(1,L); % Real valued Gaussina Primary User Signal 
 s = sqrt(snr).*bpsk_w;
 %s = 1/sqrt(2)*sqrt(snr).*randn(1,L)+randn(1,L);
 y = s + n;y2 = s + n2; % Received signal at SU
 energy = abs(y).^2;  energy2 = abs(y2).^2;% Energy of received signal over N samples
 energy_fin =(1/L).*sum(energy);energy_fin2 =(1/L).*sum(energy2); % Test Statistic for the energy detection
 thresh(m) = (qfuncinv(Pf(m))./sqrt(L))+ 1;  thresh2(m) = (qfuncinv(Pf(m)).*1.02./sqrt(L))+ 1.02;thresh3(m) = ((qfuncinv(Pf(m))./sqrt(L))+ 1)./1.002;thresh4(m) = ((qfuncinv(Pf(m)).*1.02./sqrt(L))+ 1.02)./1.001;% Theoretical value of Threshold, refer, Sensing Throughput Tradeoff in Cognitive Radio, Y. C. Liang
 if(energy_fin >= thresh(m))  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) counter by 1
     i = i+1;
 end
  if(energy_fin2 >= thresh2(m))  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) counter by 1
     i2 = i2+1;
  end
  if(energy_fin >= thresh3(m))  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) counter by 1
     i3 = i3+1;
  end
  if(energy_fin2 >= thresh4(m))  % Check whether the received energy is greater than threshold, if so, increment Pd (Probability of detection) counter by 1
     i4 = i4+1;
  end
  
end
Pd(m) = i/kk; 
Pd2(m) = i2/kk;
Pd3(m) = i3/kk;
Pd4(m) = i4/kk;

end


figure 
plot(Pf,Pd,Pf,Pd2,Pf,Pd3,Pf,Pd4)

hold on


title('Figure 4. ROC curves of energy detection scheme with no noise uncertainty,with noise uncertainty, and with dynamic threshold')
ylabel('P_D')
xlabel('P_{FA}')
legend('no noise uncertainty','with noise uncertainty','dynamic threshold','both noise uncertainty and dynamic threshold','Location','southeast')
hold on

