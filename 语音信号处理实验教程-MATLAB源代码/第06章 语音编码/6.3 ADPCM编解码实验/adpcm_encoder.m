% APDCM���뺯��
function code=adpcm_encoder(x,sign_bit)
x = x*128; 
len = length(x); 

% ���ɲ������ұ� 
index = [-1 4]; 
currentIndex = 2; 
startval = 1;                           % ���� 1
endval = 127;                       % ����startval��С��128- 1 
base = exp( log(2)/8 );         % �����(1,2)�� 

 % ���Ʋ��� 
const = startval/base; 
numSteps = round( log(endval/const) / log(base) ); 
n = 1:numSteps; 
base = exp( log(endval/startval) / (numSteps-1) ); 
const = startval/base; 
table2 = round( const*base.^n ); 

ss = zeros(1,len);
ss(1) = table2(1); 
z = zeros(1,len);
code = zeros(1,len); 
neg = 0; 
  
for n = 2:len 
    d(n) = x(n) - z(n-1);

    if (d(n) < 0) 
        neg = 1; 
        code(n) = code(n) + sign_bit; 
        d(n) = -d(n); 
    else 
        neg = 0; 
    end 

    if (d(n) >= ss(n-1))
        code(n) = code(n) + 1; 
    end 
    %������������ 
    if (neg) 
        temp = code(n) - sign_bit; 
    else
        temp = code(n); 
    end 

    temp2 = (temp+.5)*ss(n-1); 
    if (neg) 
        temp2 = -temp2; 
    end 
    z(n) = z(n-1) + temp2;
    if (z(n) > 127) 
        z(n) = 127; 
    elseif (z(n) < -127)
        z(n) = -127; 
    end 

    % �����µĲ���
    temp = temp + 1;
    currentIndex = currentIndex + index(temp);
    if (currentIndex < 1) 
        currentIndex = 1; 
    elseif (currentIndex > numSteps) 
        currentIndex = numSteps; 
    end 
    ss(n) = table2(currentIndex); 
end 
