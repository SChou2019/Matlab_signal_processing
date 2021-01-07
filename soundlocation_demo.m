% sound loaction
clc
clear all;
close all
%��ʼ������
c = 334;%ʵ�����¶����
fs = 1000;%����Ƶ��
T = 0.1; 
t = 0:1/fs:T;  % ʱ�� [0, 0.1]
L = length(t); % ʱ�䳤��:101
f = 500; %����Ȥ��Ƶ��
w = 2*pi*f; %ת���ɽ�Ƶ��
k = w/c;%����k
% ��Ԫλ��
M = 17;%��Ԫ����
% Nmid = 12;      % �ο���
% d = 3;         % ��Ԫ���
% m = (0:1:M-1) 
yi = zeros(M,1);%��Ԫ��λ��yƽ����
zi = [ 0; 3; 6; 9;12;15;18;21;24;12;12;12;12;12;12;12;12];
xi = [12;12;12;12;12;12;12;12;12; 0; 3; 6; 9;15;18;21;24];
%xi = xi.'      % ������ m*d ��Ԫ��*��Ԫ���
figure(1)
plot(xi,zi,'r*');
title('ʮ������˷�����')

%��Դλ�ã�ʵ��Ӧ���в�֪����ɨ��ǶȺͼ��ȷ��

x1 = 12;
y1 = 10;
z1 = 12; %��Դλ�� ��12,10,12�� x,zΪˮƽ�棬�����ɸ�������ȷ��
 
x2 = 12;  % array center
y2 = 0;
z2 = 12;
Ric1 = sqrt((x1-xi).^2+(y1-yi).^2+(z1-zi).^2);%��Դ������Ԫ�ľ���
Ric2 = sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2); %sound to array center:10
Rn1 = Ric1 - Ric2; %��Դ������Ԫ��ο���Ԫ�����̲�ʸ��


% Ric1 = sqrt((x1-xi).^2+(y1-yi).^2+(z1-zi).^2); % ��Դ������Ԫ�ľ���
% Ric2 = sqrt((x1-x2).^2+(y1-y2).^2+(z1-z2).^2); %sound to array center:10
% Rn1 = Ric1 - Ric2; %��Դ������Ԫ��ο���Ԫ�����̲�ʸ��

s1 = cos(2*w*t);%�ο���Զ���յ���ʸ��
Am = 10^(-1);%���
n1 = Am * (randn(M, L) + 1i*randn(M, L)); % ����Ԫ��˹������
p1 = zeros(M,L);
% ����Ԫ����ʱ��ͣ�ͬʱ��ĸ���Ԫ���յ���ѹ�����źţ��Լ�Э�������
for k1 = 1:M
    p1(k1,:) =Ric2/Ric1(k1)*s1.*exp(-1i*w*Rn1(k1)/c);
    %p1(k1,:) = Ric2/Ric1(k1) * s1.*exp(-j*w*Rn1(k1)/c);
end


p = p1 + n1;%��������
R = p*p'/L;%�������ݵ���Э����


%ɨ�淶Χ
%��������0.1��ɨ�跶Χ20*20ƽ�棬˫��ѭ���õ�M*1ʸ���������õ������׾���cross spectrum matrix��
%��DSP���ۣ��õ�������
step_x = 0.1;
step_z = 0.1;
y = y1;
x = (9:step_x:15);
z = (9:step_z:15);
for k1 = 1:length(z)
    for k2 = 1:length(x)
        Ri = sqrt((x(k2)-xi).^2+(y-yi).^2+(z(k1)-zi).^2);  % ��ɨ��㵽����Ԫ�ľ۽�����ʸ��
        Ri2 = sqrt((x(k2)-x2).^2+(y-y2).^2+(z(k1)-z2).^2);  % 10.8628
        Rn = Ri - Ri2;%���̲�
        b = exp(-1i*w*Rn/c);%���۽�����ʸ��
        Pcbf(k1,k2) = abs(b'*R*b);% ???CSM,��ؼ�,(1,18)*(18,18)*(18,1)
    end
end



%��һ������
for k1 = 1:length(z)
    pp(k1) = max(Pcbf(k1,:));%Pcbfÿ�����ֵ
end
Pcbf = Pcbf/max(pp);%ȫ�ֹ�һ��
figure(2)
surf(x,z,Pcbf)
xlabel('x(m)');zlabel('z(m)')
title('��ά����Դͼ')
colorbar
 
figure(3)
pcolor(x,z,Pcbf);
shading interp;
xlabel('x(m)');
ylabel('z(m)');
title('����Դͼ')
colorbar


