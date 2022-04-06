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


Gp = 0.85; Gs =0.15 ; % filter specifications
Wp = 1; Ws = wl_s;
fp = Wp/(2*pi);
fs = Ws/(2*pi);
ep = sqrt(1/Gp^2 - 1); es = sqrt(1/Gs^2 - 1); % ripples εp =

k = Wp/Ws; % k = 
k1 = ep/es; % k1 = 

[K,Kp] = ellipk(k); % elliptic integrals K = , K' = 
[K1,K1p] = ellipk(k1); % elliptic integrals K1 = , K'1 = 

Nexact = (K1p/K1)/(Kp/K); N = ceil(Nexact); % Nexact = , N = 

k = ellipdeg(N,k1); % recalculated k = 

fs_new = fp/k; % new stopband fs = 

L = floor(N/2); r = mod(N,2); q = (1:L); % L = , r = , i = [1; ]
u = (2*q-1)/N; zeta_i = cde(u,k); % ui = [; ], ζi = [;]

za = Wp * 1j./(k*zeta_i); % filter zeros

v0 = -1j*asne(1j/ep, k1)/N; % v0 = 

pa = Wp * 1j*cde(u-1j*v0, k); % filter poles
pa0 = Wp * 1j*sne(1j*v0, k);

s= tf('s');
z=tf('z');
lpf = ((s-za(1))*(s-conj(za(1))))/(((s-pa0)^r)*(s-pa(1))*(s-conj(pa(1))));
norm = (za(1)*conj(za(1)))/((conj(pa(1)))*(pa(1))*(pa0^r));
dc_gain  = Gp^(1-r);
lpf = dc_gain*lpf/norm;
[num,den]  = tfdata(lpf,'v');
freqs(num,den);
[h,w] =freqs(num,den);
plot(w,abs(h));
figure;
plot(w,angle(h));

bsf =  ((((B*s)/(s^2+W0^2))-za(1))*(((B*s)/(s^2+W0^2))-conj(za(1))))/((((B*s)/(s^2+W0^2))-pa0)*(((B*s)/(s^2+W0^2))-pa(1))*(((B*s)/(s^2+W0^2))-conj(pa(1))));
bsf = bsf*dc_gain/(norm);
 
[num1,den1] = tfdata(bsf,'v');
[h1,w1]= freqs(num1,den1);
figure;
plot(w1,abs(h1));
figure;
plot(w1,angle(h1));
 
 
 
 bsf_dis = ((((B*((z-1)/(z+1)))/(((z-1)/(z+1))^2+W0^2))-za(1))*(((B*((z-1)/(z+1)))/(((z-1)/(z+1))^2+W0^2))-conj(za(1))))/((((B*((z-1)/(z+1)))/(((z-1)/(z+1))^2+W0^2))-pa0)*(((B*((z-1)/(z+1)))/(((z-1)/(z+1))^2+W0^2))-pa(1))*(((B*((z-1)/(z+1)))/(((z-1)/(z+1))^2+W0^2))-conj(pa(1))));
 bsf_dis = bsf_dis*dc_gain/(norm);
 [num2,den2] = tfdata(bsf_dis,'v');
fvtool(num2,den2);