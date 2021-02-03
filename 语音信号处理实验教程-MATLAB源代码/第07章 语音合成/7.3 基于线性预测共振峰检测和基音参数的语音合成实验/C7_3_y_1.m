%ʵ��Ҫ�󣺻�������Ԥ�⹲�����ͻ��������������ϳ�
clear all; clc; close all;

[xx,fs]=wavread('C7_3_y.wav');                     % ��ȡ�ļ�
xx=xx-mean(xx);                           % ȥ��ֱ������
x1=xx/max(abs(xx));                       % ��һ��
x=filter([1 -.99],1,x1);                  % Ԥ����
N=length(x);                              % ���ݳ���
time=(0:N-1)/fs;                          % �źŵ�ʱ��̶�
wlen=240;                                 % ֡��
inc=80;                                   % ֡��
overlap=wlen-inc;                         % �ص�����
tempr1=(0:overlap-1)'/overlap;            % б���Ǵ�����w1
tempr2=(overlap-1:-1:0)'/overlap;         % б���Ǵ�����w2
n2=1:wlen/2+1;                            % ��Ƶ�ʵ��±�ֵ
wind=hamming(wlen);                       % ������
X=enframe(x,wlen,inc)';                   % ��֡
fn=size(X,2);                             % ֡��
Etemp=sum(X.*X);                          % ����ÿ֡������
Etemp=Etemp/max(Etemp);                   % ������һ��
T1=0.1; r2=0.5;                           % �˵������
miniL=10;                                 % �л������֡��
mnlong=5;                                 % Ԫ���������֡��
ThrC=[10 15];                             % ��ֵ
p=12;                                     % LPC�״�
frameTime=FrameTimeC(fn,wlen,inc,fs);     % ����ÿ֡��ʱ��̶�
Doption=0;                                

% ������-���취�������
[voiceseg,vosl,SF,Ef,period]=pitch_Ceps(x,wlen,inc,T1,fs); %���ڵ��׷��Ļ������ڼ��
Dpitch=pitfilterm1(period,voiceseg,vosl);       % ��T0����ƽ�����������������T0

%% �������ȡ
for i=1:length(SF)
    [Frmt(:,i),Bw(:,i),U(:,i)]=Formant_Root(X(:,i),p,fs,3);
end
%% �����ϳ�
zint=zeros(2,4);                          % ��ʼ��
tal=0;
for i=1 : fn
    yf=Frmt(:,i);                         % ȡ��i֡�����������Ƶ�ʺʹ���
    bw=Bw(:,i);
    [an,bn]=formant2filter4(yf,bw,fs);    % ת�����ĸ������˲���ϵ��
    synt_frame=zeros(wlen,1);
    
    if SF(i)==0                           % �޻�֡
        excitation=randn(wlen,1);         % ����������
        for k=1 : 4                       % ���ĸ��˲�����������
            An=an(:,k);
            Bn=bn(k);
            [out(:,k),zint(:,k)]=filter(Bn(1),An,excitation,zint(:,k));
            synt_frame=synt_frame+out(:,k); % �ĸ��˲������������һ��
        end
    else                                  % �л�֡
        PT=round(Dpitch(i));              % ȡ����ֵ
        exc_syn1 =zeros(wlen+tal,1);      % ��ʼ�����巢����
        exc_syn1(mod(1:tal+wlen,PT)==0)=1;% �ڻ������ڵ�λ�ò������壬��ֵΪ1
        exc_syn2=exc_syn1(tal+1:tal+inc); % ����֡��inc�����ڵ��������
        index=find(exc_syn2==1);
        excitation=exc_syn1(tal+1:tal+wlen);% ��һ֡�ļ�������Դ
        
        if isempty(index)                 % ֡��inc������û������
            tal=tal+inc;                  % ������һ֡��ǰ�����
        else                              % ֡��inc������������
            eal=length(index);            % �����м�������
            tal=inc-index(eal);           % ������һ֡��ǰ�����
        end
        for k=1 : 4                       % ���ĸ��˲�����������
            An=an(:,k);
            Bn=bn(k);
            [out(:,k),zint(:,k)]=filter(Bn(1),An,excitation,zint(:,k));
            synt_frame=synt_frame+out(:,k); % �ĸ��˲������������һ��
        end
    end
    Et=sum(synt_frame.*synt_frame);       % �����������ϳ�����
    rt=Etemp(i)/Et;
    synt_frame=sqrt(rt)*synt_frame;
        if i==1                           % ��Ϊ��1֡
            output=synt_frame;            % ����Ҫ�ص����,�����ϳ�����
        else
            M=length(output);             % �����Ա����ص���Ӵ���ϳ�����
            output=[output(1:M-overlap); output(M-overlap+1:M).*tempr2+...
                synt_frame(1:overlap).*tempr1; synt_frame(overlap+1:wlen)];
        end
end
ol=length(output);                        % �����output�ӳ����������ź�xx�ȳ�
if ol<N
    output=[output; zeros(N-ol,1)];
end
% �ϳ�����ͨ����ͨ�˲���
out1=output;
out2=filter(1,[1 -0.99],out1);
b=[0.964775   -3.858862   5.788174   -3.858862   0.964775];
a=[1.000000   -3.928040   5.786934   -3.789685   0.930791];
output=filter(b,a,out2);
output=output/max(abs(output));

subplot 211; plot(time,x1,'k'); title('ԭʼ��������');
axis([0 max(time) -1 1.1]); xlabel('ʱ��/s'); ylabel('��ֵ')
subplot 212; plot(time,output,'k');  title('�ϳ���������');
axis([0 max(time) -1 1.1]); xlabel('ʱ��/s'); ylabel('��ֵ')

