function f=mel2freq(mel)
%MEL2FREQ:HzΪ��λ��Ƶ�ʵ�MelΪ��λ��Ƶ��ת��
%��mel2freq

f=700*(exp(mel/1125)-1);