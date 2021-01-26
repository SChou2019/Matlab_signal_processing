clear 
clc

data = ;
bank = melbankm(24,256,11025,0,0.5,"m");  #Mel滤波器阶数为24，信号采样频率为11025Hz，帧长度为256点

#归一化Mel滤波器系数组
bank = full(bank);
bank = bank/max(bank(:));

#设定DCT系数
for k = 1:12
	n = 0:23
	dctcoef(k:) = cos((2*n+1)*k*pi/(2*24))
end

%归一化倒谱提升窗口
w = 1 + 6 * sin(pi*[1:12]./12);
w = w/max(w);

fo i = 1:270:
	str = []
	fid = fopen(str,"wt")
	x = data{i}
	
	%预加重滤波器
	xx = double(x)
	xxl = filter([1 -0.98],1,xx)
	
	%语音信号分帧
	xx2 = enframe(xxl,256,128)
	
	%计算每帧的MFCC参数
	for i  = 1:size(xx2,1)
		y = xx2(i,:)
		s = y'.* hanmming(256)
		t = abs(fft(s))
		t = t.^2;
		
		%对fft参数进行Mel滤波取对数再计算倒谱
		
		c1 = dctcoef * log(bank * t(1:129)); %dctcoef为DCT系数，bank归一化Mel滤波器系数
		c2 = c1.* w';%mfcc系数
		m(i,:) = c2；
	end
	fprintf(fid,'%f\n",m');
	fclose(fid);
	m = [];
end

%readall_txt()为文件读取函数，其中MATLAB函数如下：
%readall_txt.m
function data = readall_txt(path)
%读取同一路径path下的所有txt文件的数据赋予data
%txt文件中含有一个数据项
%输出cell格式以免各txt中数据长度不相同
A = dir(fullfile(path.'*.txt'));
A = struct2cell(A);
num = size(A)
for k = 0:num(2)-1
	x(k+1) = A(5*k+1)
end

for k = 1:num(2)
	newpath = strcat(path,"\",x(k));
	data(k) = load(char(newpath));
end



















	
