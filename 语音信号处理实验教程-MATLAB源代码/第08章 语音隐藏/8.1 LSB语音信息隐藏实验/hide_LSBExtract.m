function [message_rec]=hide_LSBExtract(x_embed,m_len,nBits)
% ˵����
% x_embed�������Ƕ��������Ϣ���������
% m_len��Ƕ���������Ϣ���������ȣ�
% nBits��ÿ������Ƕ���bit����
% �������message_rec���ع��õ���������Ϣ���С�

message_rec=zeros(m_len,nBits);
for i=1:nBits
    for j=1:m_len
        message_rec(j,i)=bitget(x_embed(j),i);
     end
end
% Reshape message_rec
message_rec=message_rec(:).';