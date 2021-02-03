function [h]=rir(fs, mic, n, r, rm, src);
%     ����������Ӧ
%     fs    ����Ƶ��
%     mic ��˷����꣨��������  
%     n     ������Դ���� (2*n+1)^3 
%     r 	ǽ�ڷ���ϵ����-1<R<1��
%     rm   ����ߴ磨��������
%     src   ��Դ���꣨��������
%
%     h     ����������Ӧ

nn=[-n:1:n];                                          
rms=nn+0.5-0.5*(-1).^nn;                    
srcs=(-1).^(nn);                                    
xi=[srcs*src(1)+rms*rm(1)-mic(1)];      % ʽ��9-2��
yj=[srcs*src(2)+rms*rm(2)-mic(2)];      % ʽ��9-3��
zk=[srcs*src(3)+rms*rm(3)-mic(3)];      %ʽ��9-4��

[i,j,k]=meshgrid(xi,yj,zk);                         % convert vectors to 3D matrices
d=sqrt(i.^2+j.^2+k.^2);                         % ʽ��9-5��
time=round(fs*d/343)+1;                     % ʽ��9-6��
              
[e,f,g]=meshgrid(nn, nn, nn);                  % convert vectors to 3D matrices
c=r.^(abs(e)+abs(f)+abs(g));                    % ʽ��9-9��
e=c./d;                                                   % ʽ��9-10��

h=full(sparse(time(:),1,e(:)));                     % ʽ��9-11��

