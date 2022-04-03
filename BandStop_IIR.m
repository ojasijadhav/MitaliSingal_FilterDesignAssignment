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
fns1 = fs1/f_samp*2;          
fnp1 = fp1/f_samp*2;
fnp2 = fp2/f_samp*2;
fns2 = fs2/f_samp*2;
%Transformed specs using Bilinear Transformation
ws1 = tan(fs1/f_samp*pi);          
wp1 = tan(fp1/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);
ws2 = tan(fs2/f_samp*pi);

%Parameters for Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;
wl_s1 = (B*ws1)/((W0)^2 - (ws1)^2);
wl_s2 = (B*ws2)/((W0)^2 - (ws2)^2);
wl_s = min(abs(wl_s1),abs(wl_s2));


%Butterworth LPF parameters
D1 = 1/((1-delta)^2)-1;
D2 = 1/(delta)^2 - 1;

N = acosh(sqrt(D2/D1))/acosh(wl_s); %Chebyschev approximation
N = ceil(N);
wc = ((1/D1^(1/(2*N))) + (wl_s/(D2^(1/(2*N)))))/2;

syms x;
lpf2 = (1+(D1*(chebyshevT(N,x/1i))^2));
roots_lpf2 = double(solve(lpf2));

%disp(roots_lpf2)

real_lpf2 = real(roots_lpf2);
img_lpf2  = imag(roots_lpf2);


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
den_bsf = 1;
den_des = 1;
norm = 1;
for i=1:N
    den_lpf = den_lpf*(s-roots_lpf(i));
    den_bsf = den_bsf*(((B*s)/(s^2+W0^2))-roots_lpf(i));
    den_des = den_des*(((B*((z-1)/(z+1)))/(((z-1)/(z+1))^2+W0^2))-roots_lpf(i));
    norm = norm*roots_lpf(i);
end
lpf = (norm)/(den_lpf);
bsf = (norm)/(den_bsf);
des_bsf = (norm)/(den_des);
%[num,den]  = tfdata(lpf,'v');
%freqs(num,den);
%[num,den]  = tfdata(bpf,'v');
%freqs(num,den);
%
%freqs(num,den);
%discrete
[num,den]  = tfdata(des_bsf,'v');
fvtool(num,den);