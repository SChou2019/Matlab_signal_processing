function fmel=freq2mel(freq)
%FREA2MEL:HzΪ��λ��Ƶ�ʵ�MelΪ��λ��Ƶ��ת��
%��mel2freq

fmel=1125*log(1+freq/700);