Fs = 300;%����Ƶ��
Fo = 60;%�ݲ�Ƶ��
Q = 35;%Ʒ������
Wo = Fo/(Fs/2);
BW = Wo/Q;
[b,a] = iirnotch(Wo,BW);
figure(1)
freqz(b,a,1024);
[c,d] = iircomb(10,BW,"notch");
figure(2)
freqz(c,d,1024);