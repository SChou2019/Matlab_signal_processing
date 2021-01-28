% 5.1 LPC_Levinson.m
% cacluate by autocorrelation 
% Levinson-Durbin
% split into frame, analysis order  by LPC
path = "D:\学习相关\语音信号处理，阵列\[数字语音处理及MATLAB仿真（第2版）\[数字语音处理及MATLAB仿真（第2版）][张雪英][程序源代码]\第5章\sx86.txt";
%fid = fopen("sx86.txt",'r');
fid = fopen(path,'r');
p1 = fscanf(fid,'%f');
fclose(fid);
p = filter([1 -0.68],1,p1);  %Pre emphasis filtering
x = fra(320,160,p); %split into frames, 320 points each frame,overlap points 160
x = x(60,:);        %take 60th frame, and x is row vector
s = x';             % transform
N = 16;             %order
p = N;
n = length(s);      %signal length
Rp = zeros(p,1);
for i = 1:p
    Rp(i,1) = sum(s(i+1,n).*s(1:n-i));  %autocorrealation
    %Rn(i) = sum(s(1:N-i).*s(1+i:N))
end
Rp = Rp(:);        %transform
Rp_0 = s'*s;       %equal to Rn(0)
Ep = zeros(p,1) ;       %最佳线性预测反滤波能量
%i = 1 specila case
Ep_0 = Rp_0;
k(1,1) = Rp(1,1)/Rp_0;
a(1,1) = k(1,1);
Ep(1,1) = (1 - k(1,1)^2)*Ep_0;
if p>1
    for i = 2:p
        k(i,1) = (Rp(i,1)-sum(a(1:i-1,i-1).*Rp(i-1:-1:1)))/Ep(i-1,1);
        a(i,i) = k(i,1);
        Ep(i,1) = (1 - k(i,1)^2)*Ep(i-1,1);
        for j = 1:i-1
            a(j,i) = a(j,i-1) - k(i,1) * a(i-j,i-1);
        end
    end
end
c = -a(:,p); 
a1(1,1) = 1.0;
for i = 2:p+1
    a1(1,i) = c(i-1,1);
end

function f = fra(len,inc,x)
fh = fix(((size(x,1)-len)/inc)+1);
f = zeros(fh,len);
i=1;n=1;
while i<fh
    j = 1;
    while j < len
        f(i,j) = x(n);
        j = j+1;
        n = n+1;
    end
    n = n-len+inc;
    i = i + 1;
end

end
