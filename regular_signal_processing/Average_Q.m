A = [1 2 3;4 5 6];

[m,n] = size(A);
Aver_Q = zeros(m,1);
for i = 1:m
    sum = 0;
    for j = 1:n
        sum = sum + power(A(i,j),2);
    end
    Aver_Q(i) = sqrt(sum / n);
end  
