%���׷��������ƺ���
function [Val,Loc,spect]=Formant_Cepst(u,cepstL)
%Val            �����ķ�ֵ
%Loc            ������λ��
%spect          ������
%u                һ֡�����ź�
%cepstL        ��Ƶ���ϴ������Ŀ��
U=fft(u);                                                 % ��ʽ(4-26)����
wlen2=length(u)/2;                                          % ֡��
U_abs=log(abs(U(1:wlen2)));                     
Cepst=ifft(U_abs);                                    % ��ʽ(4-27)����
cepst=zeros(1,wlen2);           
cepst(1:cepstL)=Cepst(1:cepstL);              % ��ʽ(4-28)����
cepst(end-cepstL+2:end)=Cepst(end-cepstL+2:end);
spect=real(fft(cepst));                               % ��ʽ(4-30)����
[Val,Loc]=findpeaks(spect);                      % Ѱ�ҷ�ֵ





