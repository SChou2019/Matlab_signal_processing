fs = 51200;
[num,txt,raw]=xlsread('D:\software_learning\matlab_code\Book2.xlsx');
x = num(:,1);
y = num(:,2);
%[a,b] = xcorr(x,'unbiased');
[a,b] = xcorr(x,y);
plot(b/fs,a)