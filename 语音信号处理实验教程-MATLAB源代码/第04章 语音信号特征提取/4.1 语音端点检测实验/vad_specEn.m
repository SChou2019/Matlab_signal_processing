%���ط��˵��⺯��
function [voiceseg,vsl,SF,NF,Enm]=vad_specEn(x,wnd,inc,NIS,th1,th2,fs)

y=enframe(x,wnd,inc)';             % ��֡
fn=size(y,2);                           % ��֡��
if length(wnd)==1
    wlen=wnd;               % ���֡��
else
    wlen=length(wnd);
end
df=fs/wlen;                             % ���FFT��Ƶ�ʷֱ���
fx1=fix(250/df)+1; fx2=fix(3500/df)+1;  % �ҳ�250Hz��3500Hz��λ��
km=floor(wlen/8);                       % ������Ӵ�����
K=0.5;                                  % ����K
for i=1:fn
    A=abs(fft(y(:,i)));                 % ȡ��һ֡����FFT��ȡ��ֵ
    E=zeros(wlen/2+1,1);            
    E(fx1+1:fx2-1)=A(fx1+1:fx2-1);      % ֻȡ250��3500Hz֮��ķ���
    E=E.*E;                             % ��������
    P1=E/sum(E);                        % ��ֵ��һ��
    index=find(P1>=0.9);                % Ѱ���Ƿ��з����ĸ��ʴ���0.9
    if ~isempty(index), E(index)=0; end % ����,�÷�����0
    for m=1:km                          % �����Ӵ�����
        Eb(m)=sum(E(4*m-3:4*m));
    end
    prob=(Eb+K)/sum(Eb+K);              % �����Ӵ�����
    Hb(i) = -sum(prob.*log(prob+eps));  % �����Ӵ�����
end   
Enm=multimidfilter(Hb,10);              % ƽ������
Me=min(Enm);                            % ������ֵ
eth=mean(Enm(1:NIS));
Det=eth-Me;
T1=th1*Det+Me;
T2=th2*Det+Me;

[voiceseg,vsl,SF,NF]=vad_revr(Enm,T1,T2);