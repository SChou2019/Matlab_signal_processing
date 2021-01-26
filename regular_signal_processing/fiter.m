%A计权
%时域滤波器法
data = load("D:\software_learning\matlab_code\filter_test.mat");
Fs = 32768;
y = data.ExampleRotAnalysis_E_motor_(:,1);

%截断数据
y = y(1:5*Fs+1);
%
[B,A]= adsgn(Fs);
y_timefiter = filter(B,A,y);

%窗函数
window = hann(Fs);
fs = Fs;
noverlap = 0.5 * Fs;
nfft = Fs;
%[s,f,t] = spectrogram(y_timefiter,window,noverlap,f,fs);

[S,F,T,P]=spectrogram(y_timefiter,window,noverlap,nfft,fs);
%figure
%surf(T,F,abs(P));
axis tight;

view(0,90);

colormap = abs(P);
%线型平均频谱
Aver_lin = mean(colormap,2);
figure
plot(F,Aver_lin)
hold on

%能量平均
[m,n] = size(colormap);
Aver_Q = zeros(m,1);
for i = 1:m
    sum = 0;
    for j = 1:n
        sum = sum + power(colormap(i,j),2);
    end
    Aver_Q(i) = sqrt(sum / n);
end  

plot(F,Aver_Q)

%S = spectrogram(X,WINDOW,NOVERLAP,NFFT,Fs) 
function [B,A] = adsgn(Fs)
f1 = 20.598997;
f2 = 107.65265;
f3 = 737.86223;
f4 = 12194.217;
A1000 = 1.9997;
pi = 3.14159265358979;
NUMs = [(2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0];
DENs = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]); 
DENs = conv(conv(DENs,[1 2*pi*f3]),[1 2*pi*f2]);

[B,A] = bilinear(NUMs,DENs,Fs);
end

%C计权
function [B,A] = adsgn2(Fs)
f1 = 20.598997;
f4 = 12194.217;
C1000 = 0.062;
pi = 3.12159265358979;
NUMs = [ (2*pi*f4)^2*(10^(C1000/20)) 0 0 ];
DENs = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]); 

[B,A] = bilinear(NUMs,DENs,Fs); 
end
