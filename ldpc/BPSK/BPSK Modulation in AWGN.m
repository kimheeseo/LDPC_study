%URL :
% https://kr.mathworks.com/matlabcentral/fileexchange/33284-bpsk-modulation-in-awgn?s_tid=srchtitle_awgn_7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is used to calculate the bit error rate (BER) of BPSK 
% modulation scheme at the given energy per bit to noise power spectral
% density ratio (EbNo). 
%
% By Yasir Ahmed
% www.raymaps.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;                                  %Clear all variables
close all;                                  %Close all figures
l=1e6;                                      %Number of bits or symbols
EbNodB=0:2:10;                              %Range of EbNo in dB
for n=1:length(EbNodB); 
    s=2*(round(rand(1,l))-0.5);             %Random symbol generation
    w=(1/sqrt(2*10^(EbNodB(n)/10)))*randn(1,l);  %Random noise generation
    r=s+w;                                  %Received signal
    s_est=sign(r);                          %Demodulation
    BER(n)=(l-sum(s==s_est))/l;             %BER calculation
end
semilogy(EbNodB, BER,'o-');                 %Plot 
xlabel('EbNo(dB)')                          %Label for x-axis    
ylabel('BER')                               %Label for y-axis
grid on