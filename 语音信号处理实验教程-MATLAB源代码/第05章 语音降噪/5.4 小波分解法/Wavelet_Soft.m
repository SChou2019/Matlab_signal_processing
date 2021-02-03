%С������ֵ����
function signal=Wavelet_Soft(s,jN,wname)
[c,l]=wavedec(s,jN,wname);
%��Ƶ����������
first = cumsum(l)+1;
first1=first;
first = first(end-2:-1:1);
ld   = l(end-1:-1:2);
last = first+ld-1;
%--------------------------------------------------------------------------
%����ֵ
cxdsoft=c;
for j=1:jN                                  %j�Ƿֽ�߶�
    flk = first(j):last(j);                 %flk��di��c�е�����
    thr(j)=sqrt(2*log((j+1)/j))*median(c(flk))/0.6745;
    for k=0:(length(flk)-1)             %k��λ�Ƴ߶�
        djk=c(first(j)+k);              %Ϊ�˼򻯳���
        absdjk=abs(djk);
        thr1=thr(j);
        if absdjk<thr1
           djk=0;
       else
           djk=sign(djk)*(absdjk-thr1);   
       end  
       cxdsoft(first(j)+k)=djk;
   end             
end
signal=waverec(cxdsoft,l,wname);