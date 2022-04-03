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
fns1 = fs1/f_samp*2;          
fnp1 = fp1/f_samp*2;
fnp2 = fp2/f_samp*2;
fns2 = fs2/f_samp*2;
twn = trans_bw/f_samp*2;

%Transformed specs using Bilinear Transformation
ws1 = tan(fs1/f_samp*pi);          
wp1 = tan(fp1/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);
ws2 = tan(fs2/f_samp*pi);

%Parameters for Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;
wl_s1 = ((ws1)^2 - (W0)^2)/(B*ws1);
wl_s2 = ((ws2)^2 - (W0)^2)/(B*ws2);
wl_s = min(abs(wl_s1),abs(wl_s2));


%Butterworth LPF parameters
D1 = 1/((1-delta)^2)-1;
D2 = 1/(delta)^2 - 1;

N = log(sqrt(D2/D1))/log(wl_s); %Butterworth approximation
N = ceil(N);         %order
wcl = (1/D1^(1/(2*N)));
wcu =(wl_s/(D2^(1/(2*N))));
wc = ((1/D1^(1/(2*N))) + (wl_s/(D2^(1/(2*N)))))/2;

syms x;
lpf2 = (1+(x/(1i*wc))^(2*N));
roots_lpf2 = double(solve(lpf2));
%disp(roots_lpf2)

real_lpf2 = real(roots_lpf2);
img_lpf2  = imag(roots_lpf2);
scatter(real_lpf2,img_lpf2);
%disp(length(roots_lpf2))
roots_lpf = zeros([N 1]);
j=1;
for i=1:2*N
    if real_lpf2(i)<0 
        roots_lpf(j) = roots_lpf2(i);
        j=j+1;
    end
    
end
%disp(roots_lpf)
s = tf('s');
z = tf('z');
den_lpf = 1;
den_bpf = 1;
den_dis = 1;
for i=1:N
    den_lpf = den_lpf*(s-roots_lpf(i));
    den_bpf = den_bpf*(((s^2+W0^2)/(B*s))-roots_lpf(i));
    den_dis = den_dis*((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-roots_lpf(i));
end
lpf = (wc^N)/den_lpf;
bpf = (wc^N)/den_bpf;
dis_bpf = (wc^N)/den_dis;
disp(dis_bpf);s
%[num,den]  = tfdata(lpf,'v');
%freqs(num,den);
%[num,den]  = tfdata(bpf,'v');
%freqs(num,den);
%
%freqs(num,den);
%discrete
[num,den]  = tfdata(dis_bpf,'v');
fvtool(num,den);
