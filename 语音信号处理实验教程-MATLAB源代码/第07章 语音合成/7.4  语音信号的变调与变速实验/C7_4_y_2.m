% ʵ��Ҫ����������źű��ʵ��
clear all; clc; close all;

[xx,fs]=wavread('C7_4_y.wav');                     % ��ȡ�ļ�
xx=xx-mean(xx);                           % ȥ��ֱ������
x=xx/max(abs(xx));                        % ��ֵ��һ��
lx=length(x);                             % ���ݳ���
time=(0:lx-1)/fs;                         % �����Ӧ��ʱ������
wlen=240;                                 % �趨֡��
inc=80;                                   % �趨֡�Ƶĳ���  
overlap=wlen-inc;                         % �ص�����
tempr1=(0:overlap-1)'/overlap;            % б���Ǵ�����w1
tempr2=(overlap-1:-1:0)'/overlap;         % б���Ǵ�����w2
n2=1:wlen/2+1;                            % ��Ƶ�ʵ��±�ֵ
X=enframe(x,wlen,inc)';                   % ���ղ������з�֡
fn=size(X,2);                             % ��֡��
T1=0.1; r2=0.5;                           % �˵������
miniL=10;                                 % �л������֡��
mnlong=5;                                 % Ԫ���������֡��
ThrC=[10 15];                             % ��ֵ
p=12;                                     % LPC�״�
frameTime=FrameTimeC(fn,wlen,inc,fs);     % ����ÿ֡��ʱ��̶�
in=input('���������Ƶ�������ı���:','s');  % �������Ƶ����������
rate=str2num(in);

for i=1 : fn                              % ��ȡÿ֡��Ԥ��ϵ��������
    u=X(:,i);
    [ar,g]=lpc(u,p);
    AR_coeff(:,i)=ar;
    Gain(i)=g;
end

% �������
[voiceseg,vosl,SF,Ef,period]=pitch_Ceps(x,wlen,inc,T1,fs); %���ڵ��׷��Ļ������ڼ��
Dpitch=pitfilterm1(period,voiceseg,vosl);       % ��T0����ƽ�����������������T0

if rate>1, sign=-1; else sign=1; end
lmin=floor(fs/450);                       % �������ڵ���Сֵ
lmax=floor(fs/60);                        % �������ڵ����ֵ
deltaOMG = sign*100*2*pi/fs;              % ��ֵ˳ʱ�����ʱ����ת��d��
Dpitchm=Dpitch/rate;                      % ������Ļ�������
% Dfreqm=Dfreq*rate;                        % ������Ļ���Ƶ��

tal=0;                                    % ��ʼ��
zint=zeros(p,1); 
for i=1 : fn
    a=AR_coeff(:,i);                      % ȡ�ñ�֡��ARϵ��
    sigma=sqrt(Gain(i));                  % ȡ�ñ�֡������ϵ��

    if SF(i)==0                           % �޻�֡
        excitation=randn(wlen,1);         % ����������
        [synt_frame,zint]=filter(sigma,a,excitation,zint);
    else                                  % �л�֡
        PT=floor(Dpitchm(i));             % ������ֵ��Ϊ����
        if PT<lmin, PT=lmin; end          % �ж��޸ĺ������ֵ�з���
        if PT>lmax, PT=lmax; end
        ft=roots(a);                      % ��Ԥ��ϵ�����
        ft1=ft;
%���ӹ����Ƶ�ʣ�ʵ���Ϸ��ĸ�˳ʱ��ת���·��ĸ���ʱ��ת������µĸ�ֵ
        for k=1 : p
            if imag(ft(k))>0, 
                ft1(k) = ft(k)*exp(j*deltaOMG);
	        elseif imag(ft(k))<0 
                ft1(k) = ft(k)*exp(-j*deltaOMG);
	        end
        end
        ai=poly(ft1);                     % ���µĸ�ֵ�������Ԥ��ϵ��

        exc_syn1 =zeros(wlen+tal,1);      % ��ʼ�����巢����
        exc_syn1(mod(1:tal+wlen,PT)==0)=1;% �ڻ������ڵ�λ�ò������壬��ֵΪ1
        exc_syn2=exc_syn1(tal+1:tal+inc); % ����֡��inc�����ڵ��������
        index=find(exc_syn2==1);
        excitation=exc_syn1(tal+1:tal+wlen);% ��һ֡�ļ�������Դ
        
        if isempty(index)                 % ֡��inc������û������
            tal=tal+inc;                  % ������һ֡��ǰ�����
        else                              % ֡��inc������������
            eal=length(index);            % �����������
            tal=inc-index(eal);           % ������һ֡��ǰ�����
        end
        gain=sigma/sqrt(1/PT);            % ����
        [synt_frame,zint]=filter(gain,ai,excitation,zint);%�ü�������ϳ�����
    end
    
    if i==1                               % ��Ϊ��1֡
            output=synt_frame;            % ����Ҫ�ص����,�����ϳ�����
        else
            M=length(output);             % �ص����ֵĴ���
            output=[output(1:M-overlap); output(M-overlap+1:M).*tempr1+...
                synt_frame(1:overlap).*tempr2; synt_frame(overlap+1:wlen)];
        end
end
output(find(isnan(output)))=0;
ol=length(output);                        % �����output�ӳ����������ź�xx�ȳ�
if ol<lx
    output1=[output; zeros(lx-ol,1)];
else
    output1=output(1:lx);
end
bn=[0.964775   -3.858862   5.788174   -3.858862   0.964775]; % �˲���ϵ��
an=[1.000000   -3.928040   5.786934   -3.789685   0.930791];
output=filter(bn,an,output1);             % ��ͨ�˲�
output=output/max(abs(output));           % ��ֵ��һ��

subplot 211; plot(time,x,'k'); title('ԭʼ��������'); 
axis([0 max(time) -1 1]); xlabel('ʱ��/s'); ylabel('��ֵ')
subplot 212; plot(time,output,'k');  title('�ϳ���������');
xlim([0 max(time)]); xlabel('ʱ��/s'); ylabel('��ֵ')


