% --------------------------------------------------------------------------------------------------------------
% ά���˲�enhancement=Weina_Im(y_fft,framesize,framenum,length);
% MS��������������, D-D�������������, ����Ҫ���ƴ����źŹ���
% y_fft����֡������FFT�任��
% framesize��֡��
% framenum����֡����
% length����������
% enhancement����ǿ�������
% --------------------------------------------------------------------------------------------------------------

function enhancement=Weina_Im(x,wind,inc,NIS,alpha)
    Length=length(x);
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
    y_fft2=y_a.^2;                            % ������
    noise=mean(y_fft2(:,1:NIS),2);               % ����������ƽ������
   
    snr_x_q=0.96;                                    %ǰһ֡��������ȣ���ʼֵ��Ϊ0.96
    for i=1:framenum
         Mag_y=y_a(:,i);                          
         snr_h=y_fft2(:,i)./noise;%(:,i);                   %������������
         snr_x=alpha.*snr_x_q+(1-alpha).*max(snr_h-1,0);        %���������,����"D-D"��  ,framesize*1
         Hw=snr_x./(1+snr_x);                             %ά���˲�
         M=Mag_y.*Hw;                                     %ά�ɺ�ķ���ֵ
         Mn=M.*exp(1i.*y_phase(:,i));                       %������λ
         snr_x_q=M.^2./noise;%(:,i);                        %���¹��Ƶ�ǰһ֡���������
         signal(:,i)=real(ifft(Mn));
    end
enhancement=filpframe(signal',wnd,inc);

    

    
    
    
    
    
     