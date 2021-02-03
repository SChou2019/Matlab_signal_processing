function f=add_win(input_frame,win_stype)
%������f=add_win��input_frame,win��
%����˵����input_frame:����֡
%          win:�����ͣ�rect,hamming,hanning
%          f�����ؼӴ�֡��Ϊ������
%�������ܣ��������źżӴ�

fra=input_frame(:);
len=length(fra);
switch win_stype
    case 'rect'
        win=ones(len,1);
    case 'hamming'
        win=hamming(len);
    case 'hanning'
        win=hanning(len);
    otherwise
        error('wrong win_type')
        return
end
f=(fra.*win)';
        
        