function output=SpectralSubIm(signal,wind,inc,NIS,Gamma,Beta)

nwin=length(wind);
if (nwin == 1)              % �жϴ����Ƿ�Ϊ1����Ϊ1������ʾû���贰����
   W = wind;               % �ǣ�֡��=wind
   wnd=hamming(W);
else
   W = nwin;              % ��֡��=����
   wnd=wind;
end
nfft=W;

y=enframe(signal,W,inc)';
Y=fft(y,nfft);
YPhase=angle(Y(1:fix(end/2)+1,:));          %������������λ
Y=abs(Y(1:fix(end/2)+1,:)).^Gamma;      %������
numberOfFrames=size(Y,2);

N=mean(Y(:,1:NIS)')';                           %��ʼ�������׾�ֵD(k)
NRM=zeros(size(N));                           %�������������ֵ
NoiseCounter=0;
NoiseLength=9;                                  %����ƽ������


YS=Y;                                               %ƽ����ֵ
for i=2:(numberOfFrames-1)
    YS(:,i)=(Y(:,i-1)+Y(:,i)+Y(:,i+1))/3;
end

for i=1:numberOfFrames
    [NoiseFlag, SpeechFlag, NoiseCounter, Dist]=vad_LogSpec(Y(:,i).^(1/Gamma),N.^(1/Gamma),NoiseCounter);       %����Ƶ�׾����VAD���
    if SpeechFlag==0
        N=(NoiseLength*N+Y(:,i))/(NoiseLength+1);                               %���²�ƽ������
        NRM=max(NRM,YS(:,i)-N);                                                          %����������������
        X(:,i)=Beta*Y(:,i);
    else
        D=YS(:,i)-N;                                                                                 %�׼�
        if i>1 && i<numberOfFrames                                                      %��������������            
            for j=1:length(D)
                if D(j)<NRM(j)
                    D(j)=min([D(j) YS(j,i-1)-N(j) YS(j,i+1)-N(j)]);
                end
            end
        end
        X(:,i)=max(D,0);
    end
end
output=OverlapAdd2(X.^(1/Gamma),YPhase,W,inc);



