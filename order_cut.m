%阶次切片
clc
len = 141;
for i = 40:len
    %[row,col] = find(f>bandfre(1,i)&f<bandfre(2,i));
    [row,col] = find(f>bandfre(1,i)&f<bandfre(2,i));
    if isempty(col) 
        order_2th(i) = 0;
        print("数据不足")
    elseif length(col) == 1
        order_2th(i) = order_2th(i) + colormap(i,col);
    else
        for j = 1:len(col)-1
            order_2th(i) = order_2th(i) + colormap(i,col(j));
        end
    end
end