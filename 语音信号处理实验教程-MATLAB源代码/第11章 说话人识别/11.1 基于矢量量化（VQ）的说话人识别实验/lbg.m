function v=lbg(x,k)
%lbg�����lbg��ֵ�����㷨
% lbg(x,k) ����������x���ֳ�k�ࡣ���У�xΪrow*col����ÿһ��Ϊһ��������
% ÿ��������row��Ԫ�ء�
% [v1 v2 v3 ...vk]=lbg(...)����k�����࣬����viΪ�ṹ�壬vi.numΪ����
% �к���Ԫ�ظ�����vi.ele(i)Ϊ��i��Ԫ��ֵ��vi.meaΪ��Ӧ���ľ�ֵ

[row,col]=size(x);
%u=zeros(row,k);%ÿһ��Ϊһ������ֵ
epision=0.03;%ѡ��epision����
delta=0.01;
%u2=zeros(row,k);
%LBG�㷨����k������
u=mean(x,2);%��һ���������ģ������ֵ
for i3=1:log2(k)
    u=[u*(1-epision),u*(1+epision)];%˫��
    %time=0;
    D=0;
    DD=1;
    while abs(D-DD)/DD>delta   %sum(abs(u2(:).^2-u(:).^2))>0.5&&(time<=80)   %u2~=u
        DD=D;
        for i=1:2^i3            %��ʼ��
            v(i).num=0;
            v(i).ele=zeros(row,1);
        end
        for i=1:col %��i������
            distance=dis(u,x(:,i));%��i���������������ĵľ���
            [val,pos]=min(distance);
             v(pos).num=v(pos).num+1;%Ԫ�ص�������1
            if v(pos).num==1    %eleΪ��
                v(pos).ele=x(:,i);
            else
                v(pos).ele=[v(pos).ele,x(:,i)];
            end
        end
        for i=1:2^i3 
            u(:,i)=mean(v(i).ele,2);%�µľ�ֵ����
            for m=1:size(v(i).ele,2)
                D=D+sum((v(i).ele(m)-u(:,i)).^2);
            end
        end
    end
end
%u=u;
for i=1:k  %������ֵ
    v(i).mea=u(:,i);
end