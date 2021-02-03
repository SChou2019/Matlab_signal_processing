function [AngleM] = Spectrum_Method(m,s1,s2,wlen,inc,range)
% m                     Ԥ�׻��˲������ͣ�'capon','music','esprit'
%s1,s2                 ���������ź�
% wlen                ������
% inc                   ֡�� 
%range                 �Ƕȷ�Χ
% AngleM           ���ƵĽǶ�
fs=8000;
M=2;                        %��Ԫ��Ŀ
p=1;                         %��Դ�� 
d=0.5;
a=[0:1:M-1];
t=[-pi:1/fs:pi]';         %����ʱ��

wnd=hamming(wlen);
x=enframe(s1,wnd,inc);
y=enframe(s2,wnd,inc);
n_frame=size(x,1);
Rangle=-range:.1:range;

switch lower(m)
    case 'capon'
        for z=1:n_frame
            X(1,:)=x(z,:);
            X(2,:)=y(z,:);
            Rx=X*X'/length(t); 
            [E,D,V]=svd(Rx);
            Rinv=inv(Rx);                                         %Capon
            i=1;
            for  theta=-range:.1:range;
                aaa=exp(-j*pi*sin(theta*pi/180)*a);
                bbb=1/sum((abs(aaa*Rinv)).^2);
                S_theta(1,i)=bbb;  
                i=i+1;
            end
            S=10*log10(S_theta/max(S_theta));
            [xa,ya]=max(S);
            AngleM(z)=Rangle(ya);
        end
        
    case 'music'
        for z=1:n_frame
            X(1,:)=x(z,:);
            X(2,:)=y(z,:);
            Rx=X*X'/length(t);
            [E,D,V]=svd(Rx);
        	Nn=E(:,p+1:M)*E(:,p+1:M)';                  %MUSIC
            i=1;
            for  theta=-range:.1:range;
                aaa=exp(-j*pi*sin(theta*pi/180)*a);
                bbb=1/sum((abs(aaa*Nn)).^2);
                S_theta(1,i)=bbb;  
                i=i+1;
            end
            S=10*log10(S_theta/max(S_theta));
            [xa,ya]=max(S);
            AngleM(z)=Rangle(ya);
        end

    case 'esprit'                                                      %esprit
    for z=1:n_frame
        X(1,:)=x(z,:);
        X(2,:)=y(z,:);
        Rx=X*X'/length(t);
        [U,D,V]=svd(Rx);
        S=U(:,1:p);
        phi = S(1:M-1,:)\S(2:M,:);
        w=angle(eig(phi));
        AngleM(z)=asin(w/d/pi/2)*180/pi;
    end;
    
    otherwise error('Method not defined...');

end


