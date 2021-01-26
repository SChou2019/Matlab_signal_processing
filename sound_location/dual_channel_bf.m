%总共有四个函数分别如下，导入采集到的音频文件，设置好代码读取路径就可以实现了。
 %第一个函数：采集系统
 % In this program, the signals are divided in pieces of lengData samples
 % and process separately. This program calls two functions:
 % f_adap,频域自适应
 % get_angle，计算角度
 % At the end it returns the direction in degrees.
 % CONSTANT VALUES
 i=0;
 fs=44100; % Sampling frequency (Hz)
 d_micro=0.1; % Distance between microphones (m)
 c=340; % Speed of sound (m/s)
 muestrasMAX=ceil(d_micro*fs/c); % Maximum number of samples Nmax
 DESP=ceil(muestrasMAX*1.5); % Delay we insert in the micro 2 ？作用？
 % We leave 50% of margin of error.
 lengData=44100*0.5; % Number of samples in which the voice
 % signals are divided to be processed.
 % All values depend on fs and d_micro in
 % case it was necessary to change them.
 % Importing the file to process.
 signal1=audioread('60.wav')';% LOOP WHERE THE DIFFERENT PARTS OF THE IMPORTED FILE ARE PROCESSED
 for k=lengData:lengData:length(signal1)
 tic; % Measure of time tic;toc;
 signal=signal1(1,k-(lengData-1):k); % Signal in MIC B
 d=signal1(2,k-(lengData-1):k); % Signal in MIC A
 % NORMALIZATION PROCESS
 M1=max(abs(signal)); % Maximum of channel 1
 M2=max(abs(d)); % Maximum of channel 2
 M3=max(M1,M2); % Normalization value
 %signal=signal/M32; % Normalizing
 signal=signal/M3; % Normalizing
 %d=d/M32;
 d=d/M3;
 % LMS ALGORITHM
 hDESP=[zeros(1,DESP) 1]; % Filter to delay the signal DESP samples.
 d1=conv(hDESP,d);
 P=50; % Parameters of the algorithm
 mu=0.0117;
 h0=zeros(1,P); h0(1)=0; % Initialazing the adaptative filter
 [h ,y,e]=f_adap(signal,d1,h0,mu); % Recursive function calculating the
 % coefficients of the filter h(n)
 % PROCESSING THE FILTER BEFORE THE FREQUENCY ANALYSIS.
 h1=[zeros(1,DESP-muestrasMAX-3),h(DESP-muestrasMAX-2:length(h))];
 h1(DESP+muestrasMAX+2:length(h1))=0;
 h1(DESP+1)=h1(DESP+1)/2;
 [B,I]=sort(h1,'descend');
 H1=[zeros(1,I(1)-3),h(I(1)-2:I(1)+2),zeros(1,length(h)-(I(1)+2))];
 % FREQUENCY ANALYSIS TO OBTAIN THE DELAY (IN SAMPLES)
 % 1-FFT
 lh=128; % Length of the FFT
 H=fft(h1,lh); % FFT of the filter h(n)
 % 2-ANGLE(+UNWRAP)
 alpha=angle(fftshift(H)); % Obtaining the phase
 q=unwrap(angle(fftshift(H)));
 % 3-SLOPE
 M=diff(q); % Obtaining the slope of the phase
 % 4-SLOPE’S AVERAGE
 lM=length(M)+2; % The slope M1 is not a unique value,
 p1=floor(lM/2-4); % it’s an array. So we calculate the
 p2=ceil(lM/2+4); % average of the values, K.
 K=mean(M(p1:p2));
 Nprime=(-Klh/(2*pi)); % Number of samples before
 % substracting DESP.
 % 5-SAMPLES
 if Nprime<0 % Two possible cases: negative or positive
 N=Nprime+lh;
 N=N-DESP;
 else
 N=Nprime;
 N=N-DESP;
 end
 % CALLING THE FUNCTION WHICH RETURNS THE ANGLE
 angleGRAD1=get_angle(N,fs,d_micro);
 if isreal(angleGRAD1)==1 % Security measures in case
 angleGRAD1; % the number is complex
 i=i+1;
 else
 angleGRAD1=real(angleGRAD1);
 i=i+1;
 end
 timeElapsed=toc; % Time is kept in variable timeElapsed
 timeElapsed;
 i;
 end
%第二个函数：适应度函数
 % PERFORMS THE CALCULATION OF THE COEFFICIENTS OF THE FILTER h(n).
 % Inputs: - x = Signal in MIC A
 % - d = Signal in MIC B
 % - h0 = Initial filter (equals to 0)
 % - mu = Step-size
 % Outputs: - h = Desired filter
 % - y = Convolution between x and h
 % - e = Error function
 function [h,y,e] = f_adap(x,d,h0,mu)
 % Implements the LMS algorithm.
 %Inputs: x(n) Original signal
 % d(n) Delayed signal
 % h0 Original filter
 % mu Constant value
 %Outputs: h(n) Filter
 % y(n)= x(n)h(n)
 % e(n) Error function (must be zero)
 h=h0; P=length(h);
 N=length(x);
 y=zeros(1,N); e=y; % Reserve space for y[] y e[]
 rP=0-P+1;%rP的作用？整数，某个长度
 for k=P:N
 xx=x(k+rP); % Last P inputs x[k], x[k-1], … x[k-P]
% y(k)=xxh'; % Filter output: xh Convolution？？
%此处无卷积运算
 %注意查看尺寸
 y(k) = conv(xx,h');%自己添加
 e(k)=d(k)-y(k); % Error
 h=h+mu(k)*xx; % We update the filter coefficients.
 end
 end
%第三个函数：计算模型角度的
 % OBTAINS THE ANGLE BY PERFORMING A CERTAIN NUMBER OF TRIGONOMETRIC
 % CALCULATIONS. IT CALLS THE FUNCTION:
 % hiper
 % Inputs: - N = Number of samples
 % - fs = Sampling Frequency
 % - d_micro = Distance between microphones
 % Outputs: - angle = Angle in degrees
 function[angle]= get_angle(N,fs,d_micro)
 if N~=0
 j=0.1; % Steps
 x=-20:j:20; % x axis
 [y1]=hiper(N,x,-d_micro/2,fs);% Calling function hiper
 x1=round(length(x)/4);x2=round(length(x)/8);
 pendiente=(y1(x1)-y1(x2))/(j*(x1-x2)); % Slope
 if N>0
 angulorad=atan(pendiente);
 angulo1=angulorad180/pi;
 angle=-90-angulo1;
 else
 angulorad=atan(-pendiente);
 angulo1=angulorad180/pi;
 angle=90-angulo1;
 end
 else
 angle=0;
 end
 end
%第四个函数：
 % OBTAIN THE COORDINATES OF THE POSITIONS WHERE THE SPEAKER CAN BE.
 % Inputs: - muestras = Delay between signals in samples
 % - x = Values of the x-axis
 % - xA = x coordinate of MIC A
 % - fs = Sampling frequency
 % Outputs: - y1 = Values of the y coordinate of the speaker
 function [y1]= hiper(muestras,x,xA,fs)
 c=340; % speed of sound (m/sec)
 pot=2*ones(1,length(x));
 dist=muestras/fs; % distance to B prime
 y1=sqrt(dist2/4-xA2+(4*xA2/dist2-1)*x.^pot); % formula with following
 % requisitions:
 % xA=-xB; yA=yB=0
 end
