path = "C:\Users\Administrator\Desktop\export\CTC_BASE_CTC_97.5_1_112_50_30.mat";
tf = load(path);
data = tf.CTC_BASE_CTC_97;
%È¥³ý0Hz£¿
data_real = data(2:end);
%data_whole = [data_real; flipud(conj(data_real))];
%·ùÖµÐÞÕý£¿£¿
data_whole = [data_real; conj(flipud(data_real))];
coeff = real(ifft(data_whole));
figure(1)
%plot(abs(data_whole))
plot(coeff)