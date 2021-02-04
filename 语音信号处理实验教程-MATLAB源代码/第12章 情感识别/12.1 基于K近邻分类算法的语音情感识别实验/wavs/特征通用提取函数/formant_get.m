%-----------------------------------------------------------------------
%��֡��ȡ�����
%-----------------------------------------------------------------------
function [formant_frame,bw_frame]=formant_get(speech,samplerate)
M=160;%֡��
%[speech,samplerate]=wavread('e.wav');
[begin,last]=voicemark(speech);
%��ȡ�����Ƶ�ʺʹ���
[s,s_bw,s_amp]=pick_peak(speech,samplerate);
%��Ƶ�ʺʹ������ƽ��
[formant_frame,bw_frame]=smoothing(s,s_bw);

%%��ʼ�Թ������д�����������0
%------------�����жϷ�֡����Щ֡�ǷǾ�����
frame_num=size(formant_frame,1);
frame_begin=0;%��ʼ֡
frame_last=0;%����֡
%%%%%ȷ����ʼ֡
for i=1:frame_num
    if begin<=(i-1)*M+1
        break;
    end
    frame_begin=frame_begin+1;
end
%%%%%ȷ������֡
for i=1:frame_num
    if last<=(i-1)*M+1
        break;
    end
    frame_last=frame_last+1;
end
%������ʼ֡�ͽ���֡�Թ�������ݽ��д���
for i=1:frame_begin-1
    formant_frame(i,:)=zeros(1,4);
end
for i=frame_last:frame_num
    formant_frame(i,:)=zeros(1,4);
end
 frame_last=frame_last-1;
 formant_frame = formant_frame(frame_begin:frame_last,:);
%------------------------------------------------------------