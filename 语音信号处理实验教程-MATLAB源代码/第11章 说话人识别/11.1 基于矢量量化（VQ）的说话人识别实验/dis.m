function dis=dis(u,xi)
%DIS:����ŷʽ����
% dis=dis(u,xi):����xi��u�ĸ�����������ŷʽ���룬���ص�dis��
% u������������xi��ά������һ��

if isvector(xi)&&(size(u,1)~=length(xi))
    error('xi����Ϊ������ά���������u������')
end
k=size(u,2);
xi=xi(:);
dis=zeros(1,k);
for i=1:k
    ui=u(:,i);
    dis(i)=sum((xi-ui).^2);
end