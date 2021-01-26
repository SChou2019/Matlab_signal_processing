close  all
% �����źŵĲ��� 
t=0:199;
xs=5*sin(0.3*t);
figure;
subplot(2,1,1);
plot(t,xs);grid;
ylabel('��ֵ');
title('{�������Ҳ��ź�}');
 
% ��������źŵĲ���
randn('state',sum(100*clock));
xn=randn(1,200);
zn=randn(1,200);
xn=xn+zn;
subplot(2,1,2);
plot(t,xn);grid;
ylabel('��ֵ');
xlabel('ʱ��');
title('{������������ź�}');
 
% �ź��˲�
xn = xs+xn;
xn = xn.' ;   % �����ź�����
dn = xs.' ;   % Ԥ������������
M  = 23   ;   % �˲����Ľ���
rho_max = max(eig(xn*xn.'));   % �����ź���ؾ�����������ֵ
mu = rand()*(1/rho_max)   ;    % �������� 0 < mu < 1/rho_max
[yn,W,en] = LMS(xn,dn,M,mu);
 
% �����˲��������ź�
figure;
subplot(2,1,1);
plot(t,xn);grid;
ylabel('��ֵ');
xlabel('ʱ��');
title('{�˲��������ź�}');
% ��������Ӧ�˲�������ź�
subplot(2,1,2);
plot(t,yn);grid;
ylabel('��ֵ');
xlabel('ʱ��');
title('{����Ӧ�˲�������ź�}');
% ��������Ӧ�˲�������ź�,Ԥ������źź����ߵ����
figure 
plot(t,yn,'b',t,dn,'g',t,dn-yn,'r',t,xn,'m');grid;
legend('����Ӧ�˲������','Ԥ�����','���','����Ӧ�˲�������');
ylabel('��ֵ');
xlabel('ʱ��');
title('{����Ӧ�˲���}');
%��������Ȩֵ��
figure
mm=0:M-1;
plot(mm,W(:,end)','m*');grid;
title('{����Ȩֵ��}'); 


function [yn,W,en]=LMS(xn,dn,M,mu,itr)
% LMS(Least Mean Squre)�㷨
% �������:
%     xn   ������ź�����      (������)
%     dn   ����������Ӧ����    (������)
%     M    �˲����Ľ���        (����)    �˲����Ľ���������ָ����г���Ĵ���,�����Խ�ߣ��˲�Ч����Խ��
%     mu   ��������(����)      (����)    Ҫ�����0,С��xn����ؾ����������ֵ�ĵ���    
%     itr  ��������            (����)    Ĭ��Ϊxn�ĳ���,M<itr<length(xn)
% �������:
%     W    �˲�����Ȩֵ����     (����)
%          ��СΪM : itr,
%     en   �������(itr : 1)    (������)  
%     yn   ʵ���������         (������)

% ������������Ϊ4����5��
if nargin == 4                 % 4��ʱ�ݹ�����Ĵ���Ϊxn�ĳ��� 
    itr = length(xn);
elseif nargin == 5             % 5��ʱ����M<itr<length(xn)
    if itr>length(xn) || itr<M
        error('��������������С!');
    end
else
    error('������������ĸ���!');
end
% ��ʼ������
en = zeros(itr,1);             % �������,en(k)��ʾ��k�ε���ʱԤ�������ʵ����������
W  = zeros(M,itr);             % ÿһ�д���һ����Ȩ����,ÿһ�д���-�ε���,��ʼΪ0

% ��������
for k = M:itr                  % ��k�ε���
    x = xn(k:-1:k-M+1);        % �˲���M����ͷ������
    y = W(:,k-1).' * x;        % �˲��������
    en(k) = dn(k) - y ;        % ��k�ε�������� 
    % �˲���Ȩֵ����ĵ���ʽ
    W(:,k) = W(:,k-1) + 2*mu*en(k)*x;
end

% ������ʱ�˲������������
yn = inf * ones(size(xn));
for k = M:length(xn)
    x = xn(k:-1:k-M+1);
    yn(k) = W(:,end).'* x;
end 
end
