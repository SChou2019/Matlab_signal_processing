%pr1_1.m�����ɷ�����ʾ
fs=16000;                                   %����Ƶ��
duration=2;                                %ʱ�䳤��
n=duration*fs;
t=(1:n)/fs;
fprintf('Begin by pressing any key %gseconds:\n',duration);pause
fprintf('recording...\n');
y=wavrecord(n,fs);
fprintf('Finish\n');
fprintf('Press any key to play audio:\n');pause
wavplay(y,fs);
wavwrite(y,fs,'9');