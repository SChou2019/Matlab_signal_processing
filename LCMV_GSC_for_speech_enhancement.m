% LCMV-GSC for speech enhancement
% author : Xu Changlai,6/2,2019

clear all
close all

[speech , fs ] = audioread('male_female_pure_mixture.wav');
speech = speech';
[Nch,Nz] = size(speech);
Nfft =floor( fs*64/1000); % 64 ms per frame
Nbin = floor(Nfft/2+1);%ferquency point number
Nfrm = floor(Nz/Nbin)-1;%resolution 
win = sqrt(hanning(Nfft))';%why sqrt?

yout = zeros(1,Nz);
Ybin_nonclosed = zeros(1,Nbin);
 
q = zeros(Nch+1, Nbin);
pest = zeros(Nbin,1);
mu = 0.05;
alphaP = 0.9;
Yfbf = zeros(Nch,Nfft);
phi_x = zeros(Nbin);
phi_n = zeros(Nbin);
PhiN = zeros(Nch+1, Nch+1, Nbin);
PhiS = zeros(Nch+1, Nch+1, Nbin);
PhiSN = zeros(Nch+1,Nbin);

% processing
[RTF,SPP,Mark] = segmentation (speech,fs,75,64);
% resampling C 
C1 = shiftdim(RTF,2);%shitf demision
for nsr = 1 : size(C1,3)
    C2(:,:,nsr) = resample(C1(:,:,nsr),64,64);
end
C = shiftdim(C2,1);

nsrce = size(C,2);
g = [1;zeros(nsrce-1,1)];
g = flipud(g); % change the enhanced person in the case of 2 speakers
enhansp = find(1 == g); 
 for frm = 1 : Nfrm  
    %STFT
    for ch = 1 : Nch
         Y(ch ,:) = fft(win .* speech(ch ,(frm-1)*Nbin+1:(frm-1)*Nbin+Nfft),Nfft);
    end

    for bin=1:Nbin  
        w0(:,bin) = C(:,:,bin)/(C(:,:,bin)'* C(:,:,bin)) * g;   
        B(:,:,bin) = eye(Nch,Nch) - C(:,:,bin) /(C(:,:,bin)'*C(:,:,bin))*C(:,:,bin)';
    % processing      
        % FBF filtering
        Yfbf(bin) = w0(:,bin)' * Y(:,bin)/norm(w0(:,bin));    
        % BM filtering
        u(:,bin) = B(:,:,bin) * Y(:,bin);
        
        Yout(bin) = Yfbf(bin) - q(:,bin)'* [Yfbf(bin);u(:,bin)];

        % SDW-MWF
        S(:,bin) = [Yout(bin);1e-10 * ones(Nch,1)];
        N(:,bin) = [Yfbf(bin)-Yout(bin); u(:,bin)];
        PhiS(:,:,bin) = 0.98 * PhiS(:,:,bin) + 0.02 * S(:,bin) * S(:,bin)';
        PhiN(:,:,bin) = 0.98 * PhiN(:,:,bin) + 0.02 * N(:,bin) * N(:,bin)';
        PhiSN(:,bin) = 0.98 * PhiSN(:,bin) + 0.02 * N(:,bin) * (Yfbf(bin)-Yout(bin))';
        q(:,bin) = MWF(PhiS(:,:,bin), PhiN(:,:,bin),PhiSN(:,bin),1.11); 
    end
    % load all stuff
    yout((frm-1)*Nfft/2+1:(frm-1)*Nfft/2+Nfft) = yout((frm-1)*Nfft/2+1:(frm-1)*Nfft/2+Nfft) + win .* real(ifft([Yout,conj(Yout(end-1:-1:2))]));   
end
audiowrite('RTF.wav', yout,fs);

% plot
figure(2);
subplot(3,1,1);
plot(audioread('male.wav'));
hold on 
plot(Mark(1,:));
subplot(3,1,2);
plot(audioread('female.wav'));
hold on 
plot(Mark(2,:));
subplot(3,1,3);
plot(yout);

% [ scoresbefore ] = pesq( 'male.wav', 'male_female_pure_mixture.wav' );
% [ scoresafter ] = pesq( 'male.wav', 'RTF.wav' );
% [ scoresideal ] = pesq( 'male.wav', 'male.wav' );
% fprintf('scorebefore: %f\n',scoresbefore);
% fprintf('scoreafter: %f\n',scoresafter);
% fprintf('scoreideal: %f\n',scoresideal);
% fprintf(['improved PESQ socre : %f\n'],scoresafter-scoresbefore);

% the MWF beamformer
% SDWMWF:   h = (PhiX + mu * PhiN)^-1 * PhiSN(:,bin)
% input:    PhiN, (Nch, Nch, Nbin) the noise covariance matrix
%           PhiX, (Nch, Nch, Nbin) the speech covariance matrix
%           PhiSN,(Nch, Nbin) the noise&speech across-covariance matrix
%           mu, the speech distortion/noise reduction trade-off parameter
% output:   h, (Nch, Nbin)  the beamformer coefficients
% author : Xu Changlai,6/1,2019

function h = MWF(PhiX,PhiN,PhiSN,mu)
if nargin < 3
    mu = 1;                 % typical value {0, 1}
end

[Nch, ~, Nbin] = size(PhiX);
h = zeros(Nch, Nbin);

for bin = 1:Nbin
    if rcond(PhiN(:,:,bin)) < eps
     %   disp(['bin ' num2str(bin) ': Noise covariance ill-conditioned.']);
        PhiN(:,:,bin) = PhiN(:,:,bin) + 1e-10 * eye(Nch);
    end
    h(:,bin)  = (PhiX(:,:,bin) + mu * PhiN(:,:,bin)) \ PhiSN(:,bin);
end
end

function [RTF,SPP,Mark] = segmentation (speech,fs,ov,t_p_frm)

% reference : Data-Driven Source Separation Based on Simplex Analysis,2018,Bracha
%******input
% speech : source in time domain
% fs :sample rate
% t_p_frm : time delay in per frame (ms)
% ov :  overlap ov%

%******output
% RTF: Nch x nsrce x Nbin   relative transmission fuction
% nsrce : source number
% Mark: nsrce x (length of speech)    to mark segments

%author : Xu Changlai,6/2,2019

[Nch,Nz] = size(speech);
Nfft =floor( fs*t_p_frm/1000); % 64 ms per frame
Nbin = floor(Nfft/2+1);
Nov = 100 / (100-ov) ;
Lbin = floor( Nfft/Nov ); 
Nfrm = floor(Nz/Lbin)-(Nov-1);
win = sqrt(hanning(Nfft))';

for frm = 1 : Nfrm          
    %STFT
    for ch = 1 : Nch
         Y(ch ,frm,:) = fft(win .* speech(ch ,(frm-1)*Lbin+1:(frm-1)*Lbin+Nfft),Nfft);
    end
end

upbin = floor(4.8*1000/fs*Nfft);
lowbin = floor(0.3*1000/fs*Nfft+1);
for frm = 1: Nfrm
  for bin = lowbin : upbin   %  0.2 ~ 4.8 KHz
    if frm == 1
          Am(:,frm,bin-lowbin+1) = (Y(2:end,frm,bin) * Y(1,frm,bin)' + Y(2:end,frm+1,bin) * Y(1,frm+1,bin)') ./...
              (Y(1,frm,bin) * Y(1,frm,bin)' + Y(1,frm+1,bin) * Y(1,frm+1,bin)'); 
    
    elseif frm == Nfrm    
          Am(:,frm,bin-lowbin+1) = (Y(2:end,frm-1,bin) * Y(1,frm-1,bin)' + Y(2:end,frm,bin) * Y(1,frm,bin)') ./...
              (Y(1,frm-1,bin) * Y(1,frm-1,bin)' + Y(1,frm,bin) * Y(1,frm,bin)'); 
    else
          Am(:,frm,bin-lowbin+1) = (Y(2:end,frm-1,bin) * Y(1,frm-1,bin)' + Y(2:end,frm,bin) * Y(1,frm,bin)'+ Y(2:end,frm+1,bin) * Y(1,frm+1,bin)') ./...
              (Y(1,frm-1,bin) * Y(1,frm-1,bin)' + Y(1,frm,bin) * Y(1,frm,bin)'+ Y(1,frm+1,bin) * Y(1,frm+1,bin)'); 
    end
  end
  ac(:,frm) = reshape((reshape(Am(:,frm,:),Nch-1,upbin-lowbin+1))',(Nch-1)*(upbin-lowbin+1),1);
  a(:,frm) = [real(ac(:,frm));imag(ac(:,frm))];
  for n = 1 : frm
     W(frm,n) = a(:,frm)'* a(:,n)/(2*(Nch-1)*(upbin-lowbin+1));
     W(n,frm) = W(frm,n);
  end
end

% EVD on W
[U,D] = eig(W);
norm_eigv = diag(D)/D(end,end);
V = [];
for cnt = length(norm_eigv):-1:1
  if norm_eigv(cnt) < .119  %  0.11 ~ 0.128
      nsrce = length(norm_eigv) - cnt;
      break;
  end
  V = [V,U(:,cnt)];
end
% find probability vector
[~,I1] = max(sum(V.^2,2));
e(1,:) = V(I1,:); 
V1 = V - repmat(e(1,:),size(V,1),1);
[~,I2] = max(sum(V1.^2,2));
e(2,:) = V(I2,:);
if nsrce > 2
  Er = [];
  for r = 3:nsrce      
      er = e(r-1,:)-e(1,:); 
      Er = [Er er'];
      temp = pinv(Er' * Er);%pseudoantique
      P = eye(nsrce)- Er * temp * Er';
      [~,I] = max(sum((P * V1').^2,1));
      e(r,:) = V(I,:);
  end  
 end
Q = e';
SPP = (Q \ V')'; % source present probality
% clustering and estimate the RTF 
for n = 1 : nsrce
    Ydom = zeros(Nch,1,Nbin);
    Yref = zeros(1,1,Nbin);
    mark = zeros(1,Nz);
    L = find (SPP(:,n) > .96);%classic .95
    disp(L);
    for i = 1 : length(L)
      Ydom  = Ydom + Y(: ,L(i),1:Nbin) .* repmat(conj(Y(1 ,L(i),1:Nbin)),Nch,1);
      Yref  = Yref + Y(1 ,L(i),1:Nbin) .* conj(Y(1 ,L(i),1:Nbin));
      mark((L(i)-1)*Lbin+1:(L(i)-1)*Lbin+Nfft) = ones(1,Nfft);
    end 
    RTF(:,n,:) = Ydom ./ repmat(Yref,Nch,1,1);
    % making marks
    Mark(n,:) = mark/10;
end

figure(1);
plot(speech(1,:));
hold on
plot(Mark(1,:));
hold on
plot(Mark(2,:));
end