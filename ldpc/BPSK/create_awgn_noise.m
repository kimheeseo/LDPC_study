%URL
%https://kr.mathworks.com/matlabcentral/fileexchange/40450-create-awgn-noise?s_tid=srchtitle_awgn_4
clc; clear all; close all;
a=randi(1,300);
m=4;                    % m= 4 means QPSK
sx=pskmod(a,m);
for snr= -15:15;
rxsx=awgn(sx,snr);  % output at channel
rx=pskdemod(rxsx,m); % at receiver
[ber,err]=biterr(a,rx);
ber1=berawgn(snr,'psk',m,'diff');
subplot(311)
plot(snr,ber1,'*','MarkerSize',10);
title('Uncoded');
xlabel('SNR');
ylabel('BER');
hold on
subplot(312)
plot(snr,ber,'*','MarkerSize',10);
title('Coded   QPSK');
xlabel('SNR');
ylabel('BER');
hold on 
end
m=2;            % m= 2 means BPSK
sx=pskmod(a,m);
for snr= -15:15;
rxsx=awgn(sx,snr);  % output at channel
rx=pskdemod(rxsx,m); % at receiver
[ber,err]=biterr(a,rx);
ber1=berawgn(snr,'psk',m,'diff');
subplot(313)
plot(snr,ber,'*','MarkerSize',10);
title('Coded    BPSK');
xlabel('SNR');
ylabel('BER');
hold on 
end