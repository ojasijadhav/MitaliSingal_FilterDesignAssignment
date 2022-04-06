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
lpf = ((s-za(1))*(s-za(2))*(s-conj(za(1)))*(s-conj(za(2))))/((s-pa(1))*(s-pa(2))*(s-conj(pa(1)))*(s-conj(pa(2))));
norm = (za(1)*za(2)*conj(za(1))*conj(za(2)))/((conj(pa(2))*(conj(pa(1)))*(pa(2))*(pa(1))));
dc_gain  = Gp^(1-r);
lpf = dc_gain*lpf/norm;
[num,den]  = tfdata(lpf,'v');
% freqs(num,den);
[h,w] =freqs(num,den);
plot(w,abs(h));
figure;
plot(w,angle(h));

bpf = ((((s^2+W0^2)/(B*s))-za(1))*(((s^2+W0^2)/(B*s))-za(2))*(((s^2+W0^2)/(B*s))-conj(za(1)))*(((s^2+W0^2)/(B*s))-conj(za(2))))/((((s^2+W0^2)/(B*s))-pa(1))*(((s^2+W0^2)/(B*s))-pa(2))*(((s^2+W0^2)/(B*s))-conj(pa(1)))*(((s^2+W0^2)/(B*s))-conj(pa(2))));
bpf = bpf*dc_gain/(norm);

[num1,den1] = tfdata(bpf,'v');
[h1,w1]= freqs(num1,den1);
figure;
plot(w1,abs(h1));
figure;
plot(w1,angle(h1));



bpf_dis = (((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-za(1))*((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-za(2))*((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-conj(za(1)))*((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-conj(za(2))))/(((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-pa(1))*((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-pa(2))*((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-conj(pa(1)))*((((z-1)^2 + (W0*(z+1))^2)/((z^2-1)*B))-conj(pa(2))));
bpf_dis = bpf_dis*dc_gain/(norm);

[num2,den2] = tfdata(bpf_dis,'v');
[h2,w2]= freqz(num2,den2);
figure;
plot(w1,abs(h1));
figure;
plot(w1,angle(h1));
fvtool(num2,den2);
