%���׼��㺯��
function xhat = Nrceps(x)
fftxabs = abs(fft(x));
xhat = real(ifft(log(fftxabs)));
