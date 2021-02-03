%�������˵��⺯��
function [voiceseg,vsl,SF,NF,Epara]=vad_pro(x,wnd,inc,NIS,th1,th2,mode)

y=enframe(x,wnd,inc)';             % ��֡
fn=size(y,2);                           % ��֡��
if length(wnd)==1
    wlen=wnd;               % ���֡��
else
    wlen=length(wnd);
end
if mode==1
    aparam=2; bparam=1;                     % ���ò���
    etemp=sum(y.^2);                        % ��������
    etemp1=log10(1+etemp/aparam);           % ���������Ķ���ֵ
    zcr=STZcr(x,wnd,inc);                          % ������ֵ
    Ecr=etemp1./(zcr+bparam);               % ���������
    Epara=multimidfilter(Ecr,10);             % ƽ������
    dth=mean(Epara(1:(NIS)));                % ��ֵ����
    T1=1.2*dth;
    T2=2*dth;
else
     aparam=2;                                    % ���ò���
    for i=1:fn
        Sp = abs(fft(y(:,i)));                   % FFT�任ȡ��ֵ
        Sp = Sp(1:wlen/2+1);	                 % ֻȡ��Ƶ�ʲ���
        Esum(i) = log10(1+sum(Sp.*Sp)/aparam);   % �����������ֵ
        prob = Sp/(sum(Sp));		             % �������
        H(i) = -sum(prob.*log(prob+eps));        % ������ֵ
        Ef(i) = sqrt(1 + abs(Esum(i)/H(i)));     % �������ر�
    end   
    Epara=multimidfilter(Ef,10);                   % ƽ���˲� 
%     Epara=Ef;
    Me=max(Epara);                                 % Enm���ֵ
    eth=mean(Epara(1:NIS));                        % ��ʼ��ֵeth
    Det=Me-eth;                                  % ���ֵ��������ֵ
    T1=th1*Det+eth;
    T2=th2*Det+eth;
end
[voiceseg,vsl,SF,NF]=vad_forw(Epara,T1,T2);% ����ȷ���˫���޶˵���