%--------------------------------------------------------------
function prob=mixture(mix,x)
%�����������
%���룺mix--��ϸ�˹�ṹ��x--��������,1��SIZE(SIZEΪ������������)
%�����prob--�������

prob=0;
for j=1:mix.M  %MΪÿ��״̬��pdf��
    m=mix.mean(j,:);
    v=mix.var(j,:);
    w=mix.weight(j,:);
%     prob=prob+w*pdf(m,v,x);
    guass=(2*pi*prod(v))^-0.5*exp(-0.5*(x-m)./v*(x-m)');
    prob=prob+w*guass;
end

% ����realmin, �Է�ֹviterbi.m�м���log(prob)ʱ���
if prob==0, prob=realmin; end