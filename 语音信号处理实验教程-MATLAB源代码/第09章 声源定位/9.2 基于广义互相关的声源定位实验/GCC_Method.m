function [G] = GCC_Method(m,s1,s2,wnd,inc)
% m            Ԥ�׻��˲������ͣ�'standard','roth','scot','phat','ml'
%s1,s2         ���������ź�
% Fs           ����Ƶ�� (Hz)
% wnd        ��������֡��
% inc           ֡��
%
% G            ���Ƶ�ʱ��ֵ

N=wnd;
wnd=hamming(N);
x=enframe(s1,wnd,inc);
y=enframe(s2,wnd,inc);
n_frame=size(x,1);

switch lower(m)
    case 'standard'
        % ��׼GCC
    for i=1:n_frame
        x = s1(i:i+N);
        y = s2(i:i+N);
        X=fft(x,2*N-1);
        Y=fft(y,2*N-1);
        Sxy=X.*conj(Y);
        gain=1;
        Cxy=fftshift(ifft(Sxy.*gain));
        [Gvalue(i),G(i)]=max(Cxy);%������ֵmax,�����ֵ���ڵ�λ�ã��ڼ��У�location;
    end;
    
    case 'roth'
        % Roth filter
    for i=1:n_frame
        x = s1(i:i+N);
        y = s2(i:i+N);
        X=fft(x,2*N-1);
        Y=fft(y,2*N-1);
        Sxy=X.*conj(Y);
        Sxx=X.*conj(X);
        gain=1./abs(Sxx);
        Cxy=fftshift(ifft(Sxy.*gain));
        [Gvalue(i),G(i)]=max(Cxy);%������ֵmax,�����ֵ���ڵ�λ�ã��ڼ��У�location;
    end;

    case 'scot'
        % Smoothed Coherence Transform (SCOT)
    for i=1:n_frame
        x = s1(i:i+N);
        y = s2(i:i+N);
        X=fft(x,2*N-1);
        Y=fft(y,2*N-1);
        Sxy=X.*conj(Y);
        Sxx=X.*conj(X);
        Syy=Y.*conj(Y);
        gain=1./sqrt(Sxx.*Syy);
        Cxy=fftshift(ifft(Sxy.*gain));
        [Gvalue(i),G(i)]=max(Cxy);%������ֵmax,�����ֵ���ڵ�λ�ã��ڼ��У�location;
    end;
    
    case 'phat'
        % Phase Transform (PHAT)
    for i=1:n_frame
        x = s1(i:i+N);
        y = s2(i:i+N);
        X=fft(x,2*N-1);
        Y=fft(y,2*N-1);
        Sxy=X.*conj(Y);
        gain=1./abs(Sxy);
        Cxy=fftshift(ifft(Sxy.*gain));
        [Gvalue(i),G(i)]=max(Cxy);%������ֵmax,�����ֵ���ڵ�λ�ã��ڼ��У�location;
    end;

    case 'ml'
        % �����Ȼ��Ȩ����
    for i=1:n_frame
        x = s1(i:i+N);
        y = s2(i:i+N);
        X=fft(x,2*N-1);
        Y=fft(y,2*N-1);
        Sxy=X.*conj(Y);
        Sxx=X.*conj(X);
        Syy=Y.*conj(Y);
        Zxy=(Sxy.*Sxy)/(Sxx.*Syy);
        gain=(1./abs(Sxy)).*((Zxy.^2)./(1-Zxy.^2));
        Cxy=fftshift(ifft(Sxy.*gain));
        [Gvalue(i),G(i)]=max(Cxy);%������ֵmax,�����ֵ���ڵ�λ�ã��ڼ��У�location;
    end;        
    otherwise error('Method not defined...');
end


