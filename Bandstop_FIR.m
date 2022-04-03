clear all;
m = 46;
q_m = floor(0.1*m);
r_m = m - 10*q_m;
BL = 5+3*q_m + 11*r_m;
BH = BL + 25;
trans_bw = 3*10^3;

% Band Edge specifications
fp1 = BL*10^3-trans_bw;
fs1 = BL*10^3;
fs2 = BH*10^3;
fp2 = BH*10^3+trans_bw;
f_samp = 400e3;
delta = 0.15;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

Wn = [(fs1+fp1)/2 (fs2+fp2)/2]*2/f_samp;        %average value of the two paramters
N_min = ceil((A-7.95) / (2.285*(6/400)*pi));       %empirical formula for N_min

%Window length for Kaiser Window
n=N_min+15;

%Ideal bandstop impulse response of length "n"

bs_ideal =  ideal_lp(pi,n) -ideal_lp(((fp2+fs2)/f_samp)*pi,n) + ideal_lp(((fs1+fp2)/f_samp)*pi,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandStop = bs_ideal .* kaiser_win;
fvtool(FIR_BandStop);         %frequency response

%magnitude response
[H,f] = freqz(FIR_BandStop,1,1024, f_samp);
plot(f,abs(H))
grid