%���׷��������ڼ�⺯��
function [voiceseg,vsl,SF,Ef,period]=pitch_Ceps(x,wnd,inc,T1,fs,miniL)
if nargin<6, miniL=10; end
if length(wnd)==1
    wlen=wnd;               % ���֡��
else
    wlen=length(wnd);
end
y  = enframe(x,wnd,inc)';                  % ��֡
[voiceseg,vsl,SF,Ef]=pitch_vad(x,wnd,inc,T1,miniL);   % �����Ķ˵���
fn=length(SF);
lmin=fix(fs/500);                           % �������ڵ���Сֵ
lmax=fix(fs/60);                            % �������ڵ����ֵ
period=zeros(1,fn);                         % �������ڳ�ʼ��
for k=1:fn 
    if SF(k)==1                             % �Ƿ����л�֡��
        y1=y(:,k).*hamming(wlen);           % ȡ��һ֡���ݼӴ�����
        xx=fft(y1);                         % FFT
        a=2*log(abs(xx)+eps);               % ȡģֵ�Ͷ���
        b=ifft(a);                          % ��ȡ���� 
        [R(k),Lc(k)]=max(b(lmin:lmax));     % ��lmin��lmax������Ѱ�����ֵ
        period(k)=Lc(k)+lmin-1;             % ������������
    end
end