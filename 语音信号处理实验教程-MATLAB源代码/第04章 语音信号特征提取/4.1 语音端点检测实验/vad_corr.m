%����ط��˵��⺯��
function [voiceseg,vsl,SF,NF,Rum]=vad_corr(y,wnd,inc,NIS,th1,th2)

x=enframe(y,wnd,inc)';             % ��֡
fn=size(x,2);                           % ��֡��
for k=2 : fn                            % ��������غ���
    u=x(:,k);
    ru=xcorr(u);
    Ru(k)=max(ru);
end
Rum=Ru/max(Ru);                       % ��һ��
thredth=max(Rum(1:NIS));                % ������ֵ
T1=th1*thredth;
T2=th2*thredth;
[voiceseg,vsl,SF,NF]=vad_forw(Rum,T1,T2);