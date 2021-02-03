function trimmed_X = my_vad(x)
%�˵��⣻����Ϊ¼�����������Ϊ�����ź�

Ini = 0.1;          %��ʼ��Ĭʱ��
Ts = 0.01;          %����ʱ��
Tsh = 0.005;        %֡��ʱ��
Fs = 16000;         %����Ƶ��
counter1 = 0;       %�����ĸ���������Ѱ����ʼ��ͽ�����
counter2 = 0;
counter3 = 0;
counter4 = 0;
ZCRCountf = 0;      %���ڴ洢�����ʼ����
ZCRCountb = 0;     
ZTh = 40;           %������ֵ
w_sam = fix(Ts*Fs);                   %���ڳ���
o_sam = fix(Tsh*Fs);                  %֡�Ƴ���
lengthX = length(x);
segs = fix((lengthX-w_sam)/o_sam)+1;  %��֡��
sil = fix((Ini-Ts)/Tsh)+1;            %��Ĭʱ��֡��
win = hamming(w_sam);

Limit = o_sam*(segs-1)+1;             %���һ֡����ʼλ��

FrmIndex = 1:o_sam:Limit;             %ÿһ֡����ʼλ��
ZCR_Vector = zeros(1,segs);           %��¼ÿһ֡�Ĺ������
                                     
%��ʱ�����
for t = 1:segs
    ZCRCounter = 0; 
    nextIndex = (t-1)*o_sam+1;
    for r = nextIndex+1:(nextIndex+w_sam-1)
        if (x(r) >= 0) && (x(r-1) >= 0)
         
        elseif (x(r) > 0) && (x(r-1) < 0)
         ZCRCounter = ZCRCounter + 1;
        elseif (x(r) < 0) && (x(r-1) < 0)
         
        elseif (x(r) < 0) && (x(r-1) > 0)
         ZCRCounter = ZCRCounter + 1;
        end
    end
    ZCR_Vector(t) = ZCRCounter;
end

%��ʱƽ������
Erg_Vector = zeros(1,segs);
for u = 1:segs
    nextIndex = (u-1)*o_sam+1;
    Energy = x(nextIndex:nextIndex+w_sam-1).*win;
    Erg_Vector(u) = sum(abs(Energy));
end

IMN = mean(Erg_Vector(1:sil));  %��Ĭ������ֵ��������ֵ��
IMX = max(Erg_Vector);          %��ʱƽ�����ȵ����ֵ
I1 = 0.03 * (IMX-IMN) + IMN;    %I1��I2Ϊ��ʼ������ֵ
I2 = 4 * IMN;
ITL = 100*min(I1,I2);            %������ֵ����,ǰ��ϵ������ʵ��������ĵõ����ʽ��
ITU = 10* ITL;                  %������ֵ����
IZC = mean(ZCR_Vector(1:sil));  
stdev = std(ZCR_Vector(1:sil)); %��Ĭ�׶ι����ʱ�׼��

IZCT = min(ZTh,IZC+2*stdev);    %��������ֵ
indexi = zeros(1,lengthX);      
indexj = indexi;               
indexk = indexi;
indexl = indexi;

%��Ѱ����������ֵ���޵Ĳ���
for i = 1:length(Erg_Vector)
    if (Erg_Vector(i) > ITU)
        counter1 = counter1 + 1;
        indexi(counter1) = i;
    end
end
ITUs = indexi(1);        %��һ������������ֵ���޵�֡

%��Ѱ���������������޵Ĳ���
for j = ITUs:-1:1
    if (Erg_Vector(j) < ITL)
        counter2 = counter2 + 1;
        indexj(counter2) = j;
    end
end
start = indexj(1)+1;    %��һ���о���ʼ֡

Erg_Vectorf = fliplr(Erg_Vector);%��������������������ҶԳƣ������һ�������൱������ 

%�ظ���������൱���ҽ���֡
for k = 1:length(Erg_Vectorf)
    if (Erg_Vectorf(k) > ITU)
        counter3 = counter3 + 1;
        indexk(counter3) = k;
    end
end
ITUf = indexk(1);

for l = ITUf:-1:1
    if (Erg_Vectorf(l) < ITL)
        counter4 = counter4 + 1;
        indexl(counter4) = l;
    end
end

finish = length(Erg_Vector)-indexl(1)+1;%��һ���о�����֡
    
%�ӵ�һ���о���ʼ֡��ʼ���еڶ��о��������ʣ��˵���
   
BackSearch = min(start,25);
for m = start:-1:start-BackSearch+1
    rate = ZCR_Vector(m);
    if rate > IZCT
        ZCRCountb = ZCRCountb + 1;
        realstart = m;
    end
end
if ZCRCountb > 3
    start = realstart;         
                              
end

FwdSearch = min(length(Erg_Vector)-finish,25);
for n = finish+1:finish+FwdSearch
    rate = ZCR_Vector(n);
    if rate > IZCT
        ZCRCountf = ZCRCountf + 1;
        realfinish = n;
    end
end
if ZCRCountf > 3
    finish = realfinish;        
end

x_start = FrmIndex(start);      %���յ���ʼλ��
x_finish = FrmIndex(finish-1);  %���յĽ���λ��
trimmed_X = x(x_start:x_finish); 