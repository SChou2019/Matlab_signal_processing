%LPC�ڲ巨�Ĺ������ƺ���
function [F,Bw,pp,U]=Formant_Interpolation(u,p,fs)
%pp             �����ķ�ֵ
%U               Ƶ������
%F               �����Ƶ��
%Bw            ��������
%u                һ֡�����ź�
%p                LPC����
%fs                ����Ƶ��
a=lpc(u,p);                                         % ���LPCϵ��
U=lpcar2pf(a,255);                             % ��LPCϵ�����Ƶ������
df=fs/512;                                         % Ƶ�ʷֱ���

[Val,Loc]=findpeaks(U);                     % ��U��Ѱ�ҷ�ֵ
ll=length(Loc);                                  % �м�����ֵ
for k=1 : ll
    m=Loc(k);                                     % ����m-1,m��m+1
    m1=m-1; m2=m+1;
    p=Val(k);                                      % ����H(m-1),H(m)��H(m+1)
    p1=U(m1); p2=U(m2);
    aa=(p1+p2)/2-p;                         
    bb=(p2-p1)/2;
    cc=p;                                           % ��ʽ(4-34)����
    dm=-bb/2/aa;                             % ��ʽ(4-35)����
    pp(k)=-bb*bb/4/aa+cc;              % ��ʽ(4-37)����
    m_new=m+dm;
    bf=-sqrt(bb*bb-4*aa*(cc-pp(k)/2))/aa;      % ��ʽ(4-42)����
    F(k)=(m_new-1)*df;                                  % ��ʽ(4-36)����
    Bw(k)=bf*df;                                            % ��ʽ(4-43)����
end





