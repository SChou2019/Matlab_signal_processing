%**************************************************************************
%���֡��ȡ�����źŵĹ����
%����sΪһ������,�洢��ÿ֡�����ź�����Ӧ�Ĺ����ֵ
%**************************************************************************
function [s,s_bw,s_amp]=pick_peak(sound,samplerate)

%sex��1--������2--Ů��
% sex=1;
%�������źŽ���Ԥ�����õ���֡����ź�
Sw=pre(sound,samplerate);
[frame_number,point]=size(Sw);

%��֡����,EST�ĳ�ʼֵΪÿ֡��ʼֵ��ƽ��
for t=1:1:frame_number
    [ppp,bw,amp]=find_peak(Sw(t,:),samplerate,1);
    pp(t,:)=[ppp,zeros(1,4-length(ppp))];
end
for t=1:1:4
    a=find(pp(:,t)==0);
    if isempty(a)
        l=0;
    else
        l=length(a);
    end
    EST(t)=sum(pp(:,t)/(frame_number-l));
end
%     if sex==1
%         EST(1)=320;
%         EST(2)=1440;
%         EST(3)=2760;
%         EST(4)=3200;
%     elseif sex==2
%         EST(1)=480;
%         EST(2)=1760;
%         EST(3)=3200;
%         EST(4)=3520;
%     end

% EST=mean(p,1);
%s������¼��ǰ֡�Ĺ�������ֵ
s=zeros(frame_number,4);
s_bw=zeros(frame_number,4);
s_amp=zeros(frame_number,4);
%��ʼ��֡����
for t=1:1:frame_number
    r=1;
    while(s(t,1)==0 || s(t,2)==0 || s(t,3)==0)
         [p,bw,amp]=find_peak(Sw(t,:),samplerate,r);
         if r>0.88
             r=r-0.02;
         else
             break;
         end
         peak_number=length(p);
        %p_mark������ע������ķ�ֵ���Ѿ������Ϊ0��δ�������Ϊ1��
        p_mark=ones(1,peak_number);
        
        %e������ŷ�ֵ��Ԥ��ֵ�Ĳ�
        e=zeros(1,peak_number);
        %����ȡ���ķ�ֵ��Ӧ���ÿ��������λ���ϣ�Ҫʹ�ö�Ӧ�Ĺ����ֵ��ESTֵ��Ϊ�ӽ�
         for x=1:1:4
             for y=1:1:peak_number
                 e(y)=abs(p(y)-EST(x));
             end
            [a,b]=min(e);
            s(t,x)=p(b);
            s_bw(t,x)=bw(b);
            s_amp(t,x)=amp(b);
            p_mark(b)=0;
         end
         %����������������ķ�ֵ��ͬ���򽫸�ֵ��������ESTֵ�����С��λ�ã���������һ�������ֵ��ȥ
         for x=1:1:3
             for y=(x+1):1:4
                 if s(t,x)==s(t,y)
                     e(x)=abs(s(t,x)-EST(x));
                     e(y)=abs(s(t,y)-EST(y));
                     if e(x)>e(y)
                         s(t,x)=0;
                     else 
                         s(t,y)=0;
                     end
                 end
             end
         end
    
        %����û�з���ķ�ֵ
        for k=1:1:peak_number
            if p_mark(k)~=0
                if s(t,k)==0
                    s(t,k)=p(k);
                    p_mark(k)=0;
                elseif ((amp(k)>=s_amp(t,k)*0.5) && k~=4)
                    if (s(t,k+1)==0)
                        s(t,k+1)=s(t,k);
                        s_bw(t,k+1)=s_bw(t,k);
                        s_amp(t,k+1)=s_amp(t,k);
                        s(t,k)=p(k);
                        s_bw(t,k)=bw(k);
                        s_amp(t,k)=amp(k);
                        p_mark(k)=0;
                    end
                elseif (amp(k)>=s_amp(t,k)*0.5) && (k~=1)
                    if (s(t,k-1)==0)
                        s(t,k-1)=s(t,k);
                        s_bw(t,k-1)=s_bw(t,k);
                        s_amp(t,k-1)=s_amp(t,k);
                        s(t,k)=p(k);
                        s_bw(t,k)=bw(k);
                        s_amp(t,k)=amp(k);
                        p_mark(k)=0;
                    end
                end
            end
        end
        for k=1:1:4
            if (s(t,k)~=0)
                EST(k)=s(t,k);
            end
        end
    end
end