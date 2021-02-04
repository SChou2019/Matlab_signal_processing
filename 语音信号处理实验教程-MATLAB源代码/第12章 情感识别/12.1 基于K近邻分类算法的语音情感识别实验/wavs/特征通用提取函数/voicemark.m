% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% �������źŽ��ж˵���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [begin,last]=voicemark(speech)
%��һ������ ˫����������ȡһ������
%ȷ��֡��framelength
%framelength=0.025;
point=400;
c=max(abs(speech));
speech=speech/c;

%����֡��
n=floor(length(speech)/point);  
enframe=zeros(point,n);

%��֡��ÿһ֡400����
for i=1:n;
 for j=1:point;
   enframe(j,i)=speech(j+(i-1)*point);
 end
end

%����ÿһ֡�Ķ�ʱ����
for i=1:n;   
  energy_short(i)=sum(enframe(:,i).^2);
end

%��ʱ������һ��
a=max(energy_short);
energy_short=energy_short/a;

%Ѱ���������
for i=1:n;   
  if energy_short(i)>0.2 break;  
  end
end
frame_begin=i;
begin=(i-1)*point+1;

%Ѱ�������յ�
for i=n:-1:1;   
  if energy_short(i)>0.2 break; 
  end
end

last=i*point;
% if last>length(speech)
%     last=length(speech);
% end