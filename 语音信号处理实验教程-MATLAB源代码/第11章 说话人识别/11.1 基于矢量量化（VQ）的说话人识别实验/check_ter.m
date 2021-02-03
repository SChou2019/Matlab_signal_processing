function sf=check_ter(sig,f_len,f_inc,len_noise)
%CHECK_TER :check the terminal of the speech
%sig:�����ź�   f_len:֡��  f_inc��֡��
%len_noise:����֡��
%sf�������α��
f=div_frame(sig,f_len,f_inc);%��֡(���δ���
%f=enframe(sig,f_len,f_inc);
row=size(f,1);
sf=zeros(row,1);%�ź�֡���
lmin=3;%��������С����
%��hamming��
for i=1:row
    f(i,:)=add_win(f(i,:),'hamming');
end
%f2=f(1:len_noise,:);
th_ey=(st_energy(f));
th_pz=(zero_pass(f));
%ǰ��������ƽ�������͹�����
n_ey=mean(th_ey(1:len_noise));
n_pz=mean(th_pz(1:len_noise));
th2=105*n_ey;%������ֵT2
th1=108*n_ey;%������ֵT1
%plot(find(th_ey>th1))
%figure

%��������
th3=25*n_pz;%����T3
%�˵���
%flag=0;%flag=1,��ʾ�ο�ʼ��flag=0����ʾ�ν���
for i=1:row
    %j=i-1;
    if th_ey(i)>th1%�϶�������֡
        sf(i)=1;
        j=i-1;
        %������ǰth2�ж�
        while (j>=1)&&(sf(j)~=1)&&(th_ey(j)>th2)
            sf(j)=1;
            j=j-1;
        end            
        %flag=1;
    elseif th_ey(i)>th2
         %�����κ�th2�ж�
        %j=i;
        if (i-1>=1)&&(sf(i-1)==1)
            sf(i)=1;
            %j=j+1;
        end    
        %while (j-1>=1)&&(j<=row)&&(sf(j-1)==1)
         %   sf(j)=1;
          %  j=j+1;
        %end        
    end
end
%�����ж�
for i=1:row
    if th_pz(i)>th3
            %�������ж�
       j=i;
       while (j>0)&&(sf(j)~=1)&&(sf(j+1)==1)&&(th_pz(j)>th3)
             sf(j)=1;
             j=j-1;
       end    
       j=i;
       while (j<=row)&&(sf(j)~=1)&&(sf(j-1)==1)&&(th_pz(j)>th3)
            sf(j)=1;
            j=j+1;
       end        
   end
end     
index=find(sf);
len=length(index);
flag=0;
for i=1:len-1
    if index(i+1)-index(i)<lmin
        sf(index(i)+1:index(i+1)-1)=1;
    end
end
        
    
        
    
