function f=zero_pass(x)
%calculate the zero_passing ratio of x 
%ԭ�ͣ�f=zcro(x)
%����˵���� x���������ÿһ��Ϊһ֡����
%          f�������о��󣬵�i����Ϊx�е�i֡�Ĺ�����
[row col]=size(x);
f=zeros(row,1);
for i=1:row
    for j=1:col-1
        if x(i,j)*x(i,j+1)<0;
            f(i)=f(i)+1;
        end
    end
end
