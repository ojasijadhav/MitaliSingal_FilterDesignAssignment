clear all;
m = 46;
q_m = floor(0.1*m-0.0001);
r_m = m - 10*q_m;
BL = 10+5*q_m + 13*r_m;
BH = BL + 45;
trans_bw = 3*10^3;

% Band Edge specifications
fs1 = BL*10^3-trans_bw;
fp1 = BL*10^3;
fp2 = BH*10^3;
fs2 = BH*10^3+trans_bw;
f_samp = 540e3;
delta = 0.15;
Wc1 = fp1*2*pi/f_samp;
Wc2  = fp2*2*pi/f_samp;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

N_min = ceil((A-7.95) / (2.285*(6/540)*pi));           %empirical formula for N_min

%Window length for Kaiser Window
n=N_min;

%Ideal bandpass impulse response of length "n"
bp_ideal = ideal_lp(((fp2+fs2)/f_samp)*pi,n) - ideal_lp(((fs1+fp2)/f_samp)*pi,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandPass = bp_ideal .* kaiser_win;

fvtool(FIR_BandPass);         %frequency response

%magnitude response
[H,f] = freqz(FIR_BandPass,1,1024, f_samp);
plot(f,abs(H))
grid