Fs = 300;%采样频率
Fo = 60;%陷波频率
Q = 35;%品质因子
Wo = Fo/(Fs/2);
BW = Wo/Q;
[b,a] = iirnotch(Wo,BW);
figure(1)
freqz(b,a,1024);
[c,d] = iircomb(10,BW,"notch");
figure(2)
freqz(c,d,1024);