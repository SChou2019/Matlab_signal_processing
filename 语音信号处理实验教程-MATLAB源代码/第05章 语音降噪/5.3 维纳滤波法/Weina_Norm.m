% --------------------------------------------------------------------------------------------------------------
% ά���˲� enhanced=Weina_Norm(x,wind,inc,NIS,alpha,beta)
% MS��������������,��Ҫ���ƴ����źŹ���Ps
% x:���������ź�
% framesize:֡��
% overlap:֡�ص�����
% NIS������֡֡��
% alpha,beta:���Ʋ���
% ---------------------------------------------------------------------------------------------------------------
 %%
function enhanced=Weina_Norm(x,wind,inc,NIS,alpha,beta)
    nwin=length(wind);           % ȡ����
    if (nwin == 1)              % �жϴ����Ƿ�Ϊ1����Ϊ1������ʾû���贰����
       framesize= wind;               % �ǣ�֡��=win
       wnd=hamming(framesize);                      % ���ô�����
    else
       framesize = nwin;              % ��֡��=����
       wnd=wind;
    end
    
    y=enframe(x,wnd,inc)';             % ��֡
    framenum=size(y,2);                           % ��֡��
    y_fft = fft(y);                         % FFT
    y_a = abs(y_fft);                       % ��ȡ��ֵ
    y_phase=angle(y_fft);                   % ��ȡ��λ��
    y_a2=y_a.^2;                            % ������
    noise=mean(y_a2(:,1:NIS),2);               % ����������ƽ������
    signal=zeros(framesize,1);
    for i=1:framenum
         frame=y(:,i);                                 %ȡһ֡����
         y_fft=fft(frame);                 %���ź�֡y_ham���ж�ʱ����Ҷ�任,�õ�Ƶ���ź�y_fft
         y_fft2=abs(y_fft).^2;     %����Ƶ���ź�y_fftÿ֡�Ĺ�����y_w
     
         %���������׼�ȥ������
         for k=1:framesize
                if   abs( y_fft2(k) )  >=alpha*noise(k)%(k,i)
                      signal(k)=y_fft2(k)-alpha*noise(k);%(k,i);
                      if signal(k)<0
                          signal(k)=0;
                      end
                else
                      signal(k)=beta*noise(k);%*0.01;
                end
                
         end
         %����H(W)
         Hw=( signal./(signal+1*noise) ).^1 ;
         %ά���˲������
         yw(:,i)=Hw.*y_fft;
         yt(:,i)=ifft(yw(:,i));
    end
  %������λ����������ȵ�
    enhanced=filpframe(yt',wnd,inc);



  
  
  
  
  
  