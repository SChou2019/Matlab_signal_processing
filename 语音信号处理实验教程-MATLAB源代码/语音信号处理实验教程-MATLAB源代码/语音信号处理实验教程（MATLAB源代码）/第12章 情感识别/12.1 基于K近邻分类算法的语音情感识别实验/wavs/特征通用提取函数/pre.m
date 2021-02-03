%**************************************************************************
%�������źŽ���Ԥ����,�õ���֡����ź�
%**************************************************************************
function Sw=pre(sound,samplerate)
% addpath('E:\xin\speechstuff\sysbatch\Voicebox');

%�ö�ʱ�����������źŽ��ж˵���
% % % % % % % % % [begin,last]=voicemark(sound)
% % % % % % % % % sound=sound(fix((begin+(last-begin)/4)):fix((last-(last-begin)/4)));
% r=1;
%p:��ֵ��Ӧ��Ƶ��
%amp:��ֵ��Ӧ�ķ���

%�������źŽ���Ԥ����
%Y(z)/X(z)=H(z)=1-0.98z^(-1)---Y(z)/X(z)=H(z)=1-0.94z^(-1)
%Ԥ���غ�Ľ������voice��
voice(1)=sound(1);
for t=2:1:length(sound)
    voice(t)=sound(t)-0.98*sound(t-1);
end

%�������źŽ��з�֡�������ݶ�ʱ���ԣ�֡����Ϊ25ms,���ڸ�ʵ�������Fs=16kHz
%�����Ӧ��ÿ֡��400���ź���ֵ
frame_length=0.025;                                         %%���� 
point=samplerate*frame_length;
%Ϊ��ʹ֡��֮֡��ƽ�����ɣ����ý����ֶεķ���(֡��)
%�ݶ�֡����֡���ı�ֵΪ10ms/25ms=0.4
frameinc=round(point*0.4);
%�����źŵ�֡��
Sw=enframe(voice,point,frameinc);

%hamming��
for t=1:point
    hamming_window(t)=0.54-0.46*cos(2*pi*(t-1)/(point-1));
end

[x,y]=size(Sw);
for t=1:1:x
    for j=1:1:y
    Sw(t,j)=Sw(t,j)*hamming_window(j);
    end
end