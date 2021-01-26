%白噪声数据生成
%系统自带函数

wn = 10*rand(100000,1);%
sound(wn,44100)

%% mag = 100*ones(1000,1);
mag(1) = 0;
deg = rand(1000,1)*pi;
complex = mag + 1i*deg;

wn_sim = real(ifft(complex));
wn_avg = mean(wn_sim);
wn_s = std(wn_sim,1);
%figure()
%plot(1000*wn_sim)

