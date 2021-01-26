% ifft≤‚ ‘
x = [1,2,3,4,5];
Y = fft(x);
Y;
Y1 = [1+3i, 2+2i,4,5+6i,7];
X_inverse = ifft(Y1);

X_inverse