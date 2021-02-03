%LPC������Ĺ������ƺ���
function [F,Bw,U]=Formant_Root(u,p,fs,n_frmnt)
%U        Ƶ������
%F               �����Ƶ��
%Bw            ��������
%u                һ֡�����ź�
%p                LPC����
%fs                ����Ƶ��
%n_frmnt      ��������
a=lpc(u,p);                                 % ���LPCϵ��
U=lpcar2pf(a,255);                          % ��LPCϵ���������������
df=fs/512;                                  % Ƶ�ʷֱ���

const=fs/(2*pi);                            % ����  
rts=roots(a);                               % ���
k=1;                                        % ��ʼ��
yf = [];
bandw=[];
for i=1:length(a)-1                     
    re=real(rts(i));                        % ȡ��֮ʵ��
    im=imag(rts(i));                        % ȡ��֮�鲿
    formn=const*atan2(im,re);               
    bw=-2*const*log(abs(rts(i)));           % ��(4-46)�������
    
    if formn>150 & bw <700 & formn<fs/2     % �����������ܳɹ����ʹ���
        yf(k)=formn;
        bandw(k)=bw;
        k=k+1;
    end
end

[y, ind]=sort(yf);                          % ����
bw=bandw(ind);
F = zeros(1,n_frmnt);                      % ��ʼ��
Bw = zeros(1,n_frmnt); 
F(1:min(n_frmnt,length(y))) = y(1:min(n_frmnt,length(y)));   % �������ĸ�
Bw(1:min(n_frmnt,length(y))) = bw(1:min(n_frmnt,length(y))); % �������ĸ�
Bw = Bw(:);