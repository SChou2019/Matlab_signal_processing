%����ط��������ڼ�⺯��
function [vseg,vsl,SF,Ef,period]=pitch_Corr(x,wnd,inc,T1,fs,miniL)
if nargin<6, miniL=10; end
if length(wnd)==1
    wlen=wnd;               % ���֡��
else
    wlen=length(wnd);
end
y  = enframe(x,wnd,inc)';                  % ��֡
[vseg,vsl,SF,Ef]=pitch_vad(x,wnd,inc,T1,miniL);   % �����Ķ˵���
fn=length(SF);
lmin=fix(fs/500);                           % �������ڵ���Сֵ
lmax=fix(fs/60);                            % �������ڵ����ֵ
period=zeros(1,fn);                         % �������ڳ�ʼ��
for i=1 : vsl                             % ֻ���л������ݴ���
    ixb=vseg(i).begin;
    ixe=vseg(i).end;
    ixd=ixe-ixb+1;                        % ��ȡһ���л��ε�֡��
    for k=1 : ixd                         % �Ըö��л������ݴ���
        u=y(:,k+ixb-1);                   % ȡ��һ֡����
        ru= xcorr(u, 'coeff');            % �����һ������غ���
        ru = ru(wlen:end);                % ȡ�ӳ���Ϊ��ֵ�Ĳ���
        [tmax,tloc]=max(ru(lmin:lmax));   % ��Pmin��Pmax��Χ��Ѱ�����ֵ
        period(k+ixb-1)=lmin+tloc-1;      % ������Ӧ���ֵ���ӳ���
    end
end