% ʵ��Ҫ��LSB������Ϣ����ʵ��
clear all;
clc;
% ��ȡ��������
[x_org,fs,bits]=wavread('C8_1_y.wav');
% �޷��Ż���������
if (bits==16)
   x=uint16((x_org+1)*2^(bits-1));
elseif (bits==8)
   x=uint8((x_org+1)*2^(bits-1));
end
% ������������ 
load  'C8_1_y.DAT' -mat;
nBits=1;
% Ƕ��������Ϣ
[x_embed,m_len]=hide_LSBEmbed(x,message,nBits);

% ��ȡ������Ϣ
[message_rec]=hide_LSBExtract(x_embed,m_len,nBits);

% ���������Ա�
% Step 1  ���η���
figure(1);
subplot(311);
plot(x);title('ԭʼ����');
xlabel('������')
ylabel('����')
subplot(312);
plot(x_embed);title('Ƕ��������Ϣ����');
xlabel('������')
ylabel('����')
subplot(313);
plot(x-x_embed);title('����֮��');
xlabel('������')
ylabel('����')
ylim([-10 10]);

% Step 2  �ָ��ʷ���
figure(2);
subplot(211);
imshow(reshape(message,m_mess,n_mess),[0 1]);
title('ԭʼ������Ϣ');
subplot(212);
len_mess=length(message);
message_rec=message(1:len_mess);
err_rate=sum(abs(message-message_rec))/len_mess;
imshow(reshape(message_rec,m_mess,n_mess),[0 1]);
s_title=sprintf('�ָ���������Ϣ��������%.2f %%',err_rate);
title(s_title);


