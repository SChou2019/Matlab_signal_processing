clc,clear

N=4;%NΪ����;
M=4;%MΪÿ���˴�ʶ��������
len=5;
load data.mat u;

for iii=1:len
    iii;
for i=1:M
    Dstu=zeros(N,1);
    s=['TX',num2str(iii),'_',num2str(i),'.wav'];
    [x,fs]=audioread(s);
    mel=my_mfcc(x,fs)';%������������
 
    for ii=1:N   %���ii����ƥ��
        for jj=1:size(mel,2) %����������jj����������
            distance=dis(u{ii},mel(:,jj));
            Dstu(ii)=min(distance)+Dstu(ii);
        end
    end
    [val,pos]=min(Dstu);
    if val/size(mel,2)>=81
        fprintf('�����߲���ϵͳ����\n')
    else
        fprintf('������ΪSX%d\n',pos)
    end
end
end
