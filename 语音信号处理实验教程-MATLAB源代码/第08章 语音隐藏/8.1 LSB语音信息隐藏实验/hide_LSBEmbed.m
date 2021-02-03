function [x_embed,m_len]=hide_LSBEmbed(x,message,nBits)
% ˵�����������x��������������ݣ�message�Ǵ�Ƕ���������Ϣ��
% nBits��ÿ������Ƕ���bit�����������x_embed��Ƕ��������Ϣ���������
% m_len�Ƿ���Ƕ��������Ϣ���������ȡ�

% Step 1�� ȷ��Ƕ��������Ϣ����������
% ��ȡmessage����,
len=length(message);
% ����nBits�����¹���message
pads=mod(len,nBits);
if( pads ) 
    len=len+nBits-pads;
    message=[message,zeros(1,nBits-pads)];
end
m_len=len/nBits;
mess_n=reshape(message,m_len,nBits);

% Step 2�� ����������Ƕ��������Ϣ
for i=1:nBits
    for j=1:m_len
        % �������ĵ�iλǶ����Ϣ
        if(mess_n(j,i))
            x(j)=bitset(x(j),i);
        else
            x(j)=bitset(x(j),i,0);
        end
    end
end
x_embed=x;
    