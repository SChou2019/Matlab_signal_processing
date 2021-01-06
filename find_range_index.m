f = 1:1:20;
index = find(f >5 & f < 10);
a = 0;
for j = 1:length(index)
    a = a + f(index(j));
end


amp = 5;
degree = 30;
% real = amp*sin(degree*pi/180);
% imag = amp*cos(degree*pi/180);
real = amp*sind(30);
imag = amp*cosd(30);
imag = amp*cos(deg2rad(30))