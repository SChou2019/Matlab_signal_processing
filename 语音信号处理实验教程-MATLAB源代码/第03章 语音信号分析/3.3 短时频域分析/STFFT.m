%ʵ��Ҫ��һ����ʱ����Ҷ�任����
function d=STFFT(x,win,nfft,inc)
xn=enframe(x,win,inc)';
y=fft(xn,nfft);
d=y(1:(1+nfft/2),:);