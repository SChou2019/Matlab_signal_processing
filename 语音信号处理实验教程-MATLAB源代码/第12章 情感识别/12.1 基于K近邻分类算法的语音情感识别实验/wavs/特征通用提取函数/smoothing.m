%**************************************************************************
%��ÿ֡�Ĺ����ֵ����ƽ��
%**************************************************************************
function [f,bw]=smoothing(s,s_bw)
[frame_number,peak_number]=size(s);
f=zeros(frame_number,peak_number);
for t=1:1:peak_number
    ss=s(:,t);
    ss_bw=s_bw(:,t);
    a=find(ss==0);     %�ж���û��֡©����t�������
    l=length(a);
    if l~=0
        ss(a)=sum(ss)/(frame_number-l);
        ss_bw(a)=sum(ss_bw)/(frame_number-l);
    end
% %     %�ж�©����֡�ǵ���֡��������֡
% %     %���©�����ǵ���֡������������֡��ƽ��ֵ�����֡��ֵ
% %     if l==1
% %         ss(a)=sum(ss)/(frame_number-l);
% %     elseif l>1
% %         for j=1:1:l-1
% %             b(j)=a(j+1)-a(j);
% %         end
% %         %��һ֡�����һ֡��������
% %         if b(1)~=1
% %             ss(a(1))=sum(ss)/(frame_number-l);
% %             ss_amp(a(1))=sum(ss_amp)/(frame_number-l);
% %         end
% %         if b(l-1)~=1
% %             ss(a(l))=sum(ss)/(frame_number-l);
% %             ss_amp(a(l))=sum(ss_amp)/(frame_number-l);
% %         end
% %     
% %         %�м��֡����ͷ�жϣ���������ʱ����Ϊ����֡
% %         for j=2:1:l-2
% %             if b(j)~=1 && b(j+1)~=1
% %                 ss(a(j))=sum(ss)/(frame_number-l);
% %                 ss_amp(a(j))=sum(ss_amp)/(frame_number-l);
% %             end
% %         end
% %     end
% %     clear b;        
    %�ж��Ƿ���֡����ͻ��ֵ
    %dΪ����֡������ľ��룬d(j)=ss(j+1)-ss(j);
    d=zeros(1,frame_number-1);
    for j=1:1:frame_number-1
    %      if ss(j+1)~=0 && ss(j)~=0
            d(j)=ss(j+1)-ss(j);
%         end
    end
    
    %thetaΪ����ֵ����Ϊ240hz
    theta=240;
    %һ֡����ͻ������
    for j=3:1:frame_number-5
        if abs(d(j))>theta
            if abs(d(j-1))<theta && abs(d(j+1)+d(j-1))<theta && abs(d(j+2))<theta
                ss(j+1)=(ss(j+2)+ss(j))/2;
                ss_bw(j+1)=(ss_bw(j+2)+ss_bw(j))/2;
            elseif abs(d(j-1))<theta && abs(d(j+2)+d(j+1)+d(j))<theta && abs(d(j+3))<theta
                ss(j+1)=(ss(j+3)+ss(j))/2;
                ss_bw(j+1)=(ss_bw(j+3)+ss_bw(j))/2;
            elseif abs(d(j-1))<theta && abs(d(j+3)+d(j+2)+d(j+1)+d(j))<theta && abs(d(j+4))<theta
                ss(j+1)=(ss(j+4)+ss(j))/2;
                ss_bw(j+1)=(ss_bw(j+4)+ss_bw(j))/2;
            end
        end
    end

    %���������ƽ��    
    for j=2:1:frame_number-1
        ss(j)=0.25*ss(j-1)+0.5*ss(j)+0.25*ss(j+1);
        ss_bw(j)=0.25*ss_bw(j-1)+0.5*ss_bw(j)+0.25*ss_bw(j+1);
    end

    for j=2:1:frame_number-1
        ss(j)=0.25*ss(j-1)+0.5*ss(j)+0.25*ss(j+1);
        ss_bw(j)=0.25*ss_bw(j-1)+0.5*ss_bw(j)+0.25*ss_bw(j+1);
    end

    f(:,t)=ss;
    bw(:,t)=ss_bw;
end
         
    
    
     