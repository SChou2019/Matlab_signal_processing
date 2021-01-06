%% SRP Estimate of Direction of Arrival at Microphone Array
% Estimate the direction of arrival of a signal using the SRP-PHAT
% algorithm. 
%%
%close all
fs=16000;        % sampling frequency (arbitrary)
D=2;            % duration in seconds

L = ceil(fs*D)+1; % signal duration (samples)
n = 0:L-1;        % discrete-time axis (samples)
t = n/fs;         % discrete-time axis (sec)
x = chirp(t,0,D,fs/2)';   % sine sweep from 0 Hz to fs/2 Hz，信号用不着

%% load recorded office noise audio
noisepath = '../../noise/';
[noise2,fs] = audioread([noisepath,'音轨-2.wav']);
noise0 = audioread([noisepath,'音轨.wav']);
noise5 = audioread([noisepath,'音轨-5.wav']);
noise = [noise2,noise0,noise5];

%use a clean speech audio as desired signal
pathname = '../';
[speech ,fs] = audioread([pathname,'speech.wav']);
%scale source signal to obtain 0 dB input SNR    
speech = speech(1:length(noise0))/2;   


% x = filter(Num,1,x0);
c = 340.0;
%%
% Create the 5-by-5 microphone URA.
d = 0.042;
N = 3;
mic = phased.OmnidirectionalMicrophoneElement;
array = phased.ULA(N,d,'Element',mic);
% array = phased.URA([N,N],[d,d],'Element',mic);

%%
% Simulate the incoming signal using the |WidebandCollector| System
% object(TM).
arrivalAng = [15;0];
collector = phased.WidebandCollector('Sensor',array,'PropagationSpeed',c,...
    'SampleRate',fs,'ModulatedInput',false);
s = collector(speech,arrivalAng);
% signal = signal(1:4800,:);
signal = s+noise;

%%
path = '../../TestAudio/num3_MIC5/';
[noise7,fs] = audioread([path,'音轨-7.wav']);
noise0 = audioread([path,'音轨.wav']);
noise4 = audioread([path,'音轨-4.wav']);
signal = [noise7,noise0,noise4];
%%
t = 0;
P = zeros(1,length(-90:step:90-step));
step = 1;
tic
for i = -90:step:90-step
    [ DS, x1] = phaseshift(signal,fs,256,256,128,d,i/180*pi);
    t = t+1;
    P(t) = DS'*DS;
end
toc
[m,index] = max(P);
figure,plot(-90:step:90-step,P/max(P))
(index-90+step)*step



function [ DS, x1] = phaseshift( x,fs,N,frameLength,inc,d,angle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%frequency-domain delay-sum beamformer
%   
%      input :
%          x : input signal ,samples * channel
%          fs: sample rate
%          N : fft length,frequency bin number
%frameLength : frame length,usually same as N
%        inc : step increment
%          d : array element spacing
%      angle : incident angle
%
%     output :
%         DS : delay-sum output
%         x1 : presteered signal,same size as x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 256;
% inc = 32;
% frameLength = 256;
c = 340;
Nele = size(x,2);
omega = zeros(frameLength,1);
H = zeros(N/2+1,Nele);
% tao = d*sin(angle(1))*cos(angle(2))*[0:Nele-1]/c;     %方位角 -90 < theta <90
tao = d*sin(angle(1))*[0:Nele-1]/c;     %方位角 -90 < theta <90
yds = zeros(length(x(:,1)),1);
x1 = zeros(size(x));
for i = 1:inc:length(x(:,1))-frameLength
    for k = 2:N/2+1
        omega(k) = 2*pi*(k-1)*fs/N;        
%         H(k,:) = [1;exp(-1j*omega(k)*tao);exp(-1j*omega(k)*2*tao)];
        %对齐向量，以第一个阵元为参考，
        %例如若第一个阵元最慢(theata>0),则将第2、3、....个阵元分别延迟exp(-j*w*m*tao)
%         H(k,:) = [1;exp(-1j*omega(k)*tao);exp(-1j*omega(k)*2*tao);];
        H(k,:) = exp(-1j*omega(k)*tao);
    end
    d = fft(bsxfun(@times, x(i:i+frameLength-1,:),hamming(frameLength)'));%通过广播机制批处理
%     d = fft(x(i:i+frameLength-1,:).*hamming(frameLength)');
%     x_fft = d(1:129,:).*H;%./abs(d(1:129,:));
    x_fft=bsxfun(@times, d(1:N/2+1,:),H);
    
    % phase transformed
    x_fft = bsxfun(@rdivide, x_fft,abs(d(1:N/2+1,:)));
    yf = sum(x_fft,2);
    Cf = [yf;conj(flipud(yf(1:N/2-1)))];
    
    % 恢复延时累加的信号
    yds(i:i+frameLength-1) = yds(i:i+frameLength-1)+(ifft(Cf));
    
    
    % 恢复各路对齐后的信号
    xf  = [x_fft;conj(flipud(x_fft(1:N/2-1,:)))];
    x1(i:i+frameLength-1,:) = x1(i:i+frameLength-1,:)+(ifft(xf));
end
DS = yds/Nele;  


end
