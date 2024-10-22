clear all
close all
derad = pi/180;        % deg -> rad
radeg = 180/pi;
twpi = 2*pi;
kelm = 8;               % 阵列数量
dd = 0.5;               % space 
d=0:dd:(kelm-1)*dd;     % 
iwave = 3;              % number of DOA
theta = [10 30 60];     % 角度
snr = 10;               % input SNR (dB)
n = 500;                 % 
A=exp(-1i*twpi*d.'*sin(theta*derad));%%%% direction matrix
S=randn(iwave,n);
X=A*S;
X1 = awgn(X,snr,"measured");
Rxx = X1*X1'/n;
InvS = inv(Rxx);%%%求逆
[EV,D] = eig(Rxx);%计算特征向量和特征值
EVA = diag(D)';%矩阵转为向量，特征值
[EVA,I] = sort(EVA);%特征值排序
EVA = fliplr(EVA);
EV = fliplr(EV(:,I));

%MUSIC

for iang = 1:361
    angle(iang) = (iang-181)/2;
    phim = derad *angle(iang);
    a = exp(-j*twpi*d*sin(phim)).';%转置需要加点？
    L = iwav;
    En = EV(:,L+1:kelm);
    SP(iang) = (a'*a)/(a'*En*En'*a);
end

% 
SP = abs(SP);
SPmax = max(SP);
SP = 10*log10(SP/SPmax);
h = plot(angle,Sp);
set(h,'Linewidth',2)
xlabel('angle (degree)')
ylabel('magnitude (dB)')
axis([-90 90 -60 0])
set(gca, 'XTick',[-90:30:90])
grid on  




