%%
%Code starts from line 87
%Load attached workspace before running the code

%%

% % Recording, Playing and Writing Audio File
% % clc;close all;clear all;
% % warning off
% % recObj = audiorecorder;% audiorecorder creates an 8000 Hz, 8-bit, 1 channel audiorecorder object.
% % audiorecorder(Fs, NBITS, NCHANS) creates an audiorecorder object with 
% %     sample rate Fs in Hertz, number of bits NBITS, and number of channels NCHANS. 
% %     Common sample rates are 8000, 11025, 22050, 44100, 48000, and 96000 Hz.
% %     The number of bits must be 8, 16, or 24. The number of channels must
% %     be 1 or 2 (mono or stereo).
% % audiorecorder(Fs, NBITS, NCHANS, ID) creates an audiorecorder object using 
% %     audio device identifier ID for input.  If ID equals -1 the default input 
% %     device will be used.
% %     
% % Fs =8000 ; % Sampling frequency in hertz8000, 11025, 22050, 44100, 48000, and 96000 Hz.
% % nBits = 16 ;% 8, 16, or 24
% % nChannels = 1 ; %Number of channels--2 options--1 (mono) or 2 (stereo)
% % ID = -1; % default audio input device like Microphone
% % recObj = audiorecorder(Fs,nBits,nChannels,ID);
% %  
% % disp('Start speaking.')
% % recordblocking(recObj,0.5); 
% % recordblocking(OBJ, T) records for length of time, T, in seconds;
% %              does not return until recording is finished.
% % disp('End of Recording.');
% % play(recObj);
% % mySpeech = getaudiodata(recObj); % returns the recorded audio data as a double array
% % getaudiodata(OBJ, DATATYPE) returns the recorded audio data in
% %      the data type as requested in string DATATYPE.  Valid data types
% %      are 'double', 'single', 'int16', 'uint8', and 'int8'.
% % Write audio file
% % audiowrite('test.wav',mySpeech,Fs);
% %  

%Playing Recorded Audio file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Reading and Plotting Audio Signal with Noise %%%%%%%%%%%%%%%%%%%%%%%%%
% clc;clear all; close all
% [name,Fs] = audioread('myname.wav');
% [rollno,Fs] = audioread('myrollno.wav');
% [sec,Fs] = audioread('mysec.wav');
% [stream,Fs] = audioread('mystream.wav');
% samples=length(name);
% name=name/max(name);
% rollno=rollno/max(rollno);
% sec=sec/max(sec);
% stream=stream/max(stream);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sound(name,Fs)
% pause(2.5)
% sound(rollno,Fs)
% pause(2.5)
% sound(sec,Fs)
% pause(2.5)
% sound(stream,Fs)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% n_samples = name;
% r_samples = rollno;
% s_samples = sec;
% st_samples = stream;

%%%%%
%starting reconstruction process
% % % name_recon=0;
% % % roll_recon=0;
% % % sec_recon=0;
% % % stream_recon=0;
% % % for k=0:length(n_samples)-1
% % %     l = linspace(k,-16000+k,4*2228000); %Fs/2=2228000, our sound
% sample is 2 seconds
% % %     name_recon=name_recon+n_samples(k+1)*sinc(l);
% % %     roll_recon=roll_recon+r_samples(k+1)*sinc(l);
% % %     sec_recon=sec_recon+s_samples(k+1)*sinc(l);
% % %     stream_recon=stream_recon+st_samples(k+1)*sinc(l);
% % % end

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%% START FROM HERE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 load WorkSpace_889K_x2
%
%  Fs = length(name_recon)/2; %each of our message is 2 seconds long
% sound(name_recon(1:557:length(name_recon)),8000) % 557 =  8912000 / 16000
% pause(2.5)
% sound(roll_recon(1:557:length(roll_recon)),8000)
% pause(2.5)
% sound(sec_recon(1:557:length(sec_recon)),8000)
% pause(2.5)
% sound(stream_recon(1:557:length(stream_recon)),8000)
% Message 1 figures 
figure(1)
subplot(211)
plot(t,name_recon)
xlabel('t (sec)','fontweight','bold','fontsize',14);
ylabel('m1(t)','fontweight','bold','fontsize',14);
grid on;
title(' Name Reconstructed In Time Doman')

m1k=fft(name_recon)/Fs*2;
m1k=fftshift(abs(m1k));
fx=-Fs/2:Fs/length(m1k):((Fs/2)-(Fs/length(m1k)));

subplot(212)
stem(fx,m1k,'b','LineWidth',0.1);
 xlim([-4000 4000]);
xlabel('F(Hz)','fontweight','bold','fontsize',14);
ylabel('m1k (f)','fontweight','bold','fontsize',14);

grid on;
title('Name Reconstructed In Freq Doman')


% Message 2 figures 
figure(2)
subplot(211)
plot(t,roll_recon)
xlabel('t (sec)','fontweight','bold','fontsize',14);
ylabel('m2(t)','fontweight','bold','fontsize',14);
grid on;
title(' Roll no Reconstructed In Time Doman')

m2k=fft(roll_recon)/Fs*2;
m2k=fftshift(abs(m2k));
fx=-Fs/2:Fs/length(m2k):((Fs/2)-(Fs/length(m2k)));

subplot(212)
stem(fx,m2k,'b','LineWidth',0.1);
xlabel('F(Hz)','fontweight','bold','fontsize',14);
ylabel('m2k(f)','fontweight','bold','fontsize',14);
 xlim([-4000 4000]);
grid on;
title('Roll no Reconstructed In Freq Doman')


% Message 3 figures 
figure(3)
subplot(211)
plot(t,sec_recon)
xlabel('t (sec)','fontweight','bold','fontsize',14);
ylabel('m3(t)','fontweight','bold','fontsize',14);
grid on;
title(' Section Reconstructed In Time Doman')

m3k=fft(sec_recon)/Fs*2;
m3k=fftshift(abs(m3k));
fx=-Fs/2:Fs/length(m3k):((Fs/2)-(Fs/length(m3k)));

subplot(212)
stem(fx,m3k,'b','LineWidth',0.1);
xlabel('F(Hz)','fontweight','bold','fontsize',14);
ylabel('m3k(f)','fontweight','bold','fontsize',14);
 xlim([-4000 4000]);
grid on;
title('Section Reconstructed In Freq Doman')

% Message 4 figures 
figure(4)
subplot(211)
plot(t,stream_recon)
xlabel('t (sec)','fontweight','bold','fontsize',14);
ylabel('m4(t)','fontweight','bold','fontsize',14);
grid on;
title(' Stream Reconstructed In Time Doman')

m4k=fft(stream_recon)/Fs*2;
m4k=fftshift(abs(m4k));
fx=-Fs/2:Fs/length(m4k):((Fs/2)-(Fs/length(m4k)));

subplot(212)
stem(fx,m4k,'b','LineWidth',0.1);
xlabel('F(Hz)','fontweight','bold','fontsize',14);
ylabel('m4k(f)','fontweight','bold','fontsize',14);
 xlim([-4000 4000]);
grid on;
title('Stream Reconstructed In Freq Doman')

%% 
n=0:2*Fs-1;
t=n*(1/Fs);
B=5;

kf1 = 3300*B/max(name_recon);
kf2 = 3300*B/max(roll_recon);
kf3 = 3300*B/max(sec_recon);
kf4 = 3300*B/max(stream_recon);

freqdev1=kf1*max(name_recon);
freqdev2=kf2*max(roll_recon);
freqdev3=kf3*max(sec_recon);
freqdev4=kf4*max(stream_recon);
fc4=889000; %selected according to roll number
fc3=829000; %left a space of 60k between each fc
fc2=769000;
fc1=709000;
%Individual modulation
x1=cos(2*pi*fc1*t + 2*pi*kf1*cumtrapz(name_recon)/length(name_recon));
x2=cos(2*pi*fc2*t + 2*pi*kf2*cumtrapz(roll_recon)/length(roll_recon));
x3=cos(2*pi*fc3*t + 2*pi*kf3*cumtrapz(sec_recon)/length(sec_recon));
x4=cos(2*pi*fc4*t + 2*pi*kf4*cumtrapz(stream_recon)/length(stream_recon));
%%   Multiplexing  and Transimission 
%Frequency multiplexing
u=x1+x2+x3+x4;
% Local oscillator Freqs at receiver end
fif = 450000; 
flo1 = fc1 + fif;
flo2 = fc2 + fif;
flo3 = fc3 + fif;
flo4 = fc4 + fif;

%%    User interface 
string msg;
msgin=inputdlg('Please Enter your choice: 1,2,3 or 4                                            .','Input');
%msg = '1';
ch1 = strcmpi(msgin, '1') ;
ch2 = strcmpi(msgin, '2') ;
ch3 = strcmpi(msgin, '3') ;
ch4 = strcmpi(msgin, '4') ;
while ~(ch1 | ch2 | ch3 |ch4 )
    msgin=inputdlg('Please Enter a valid choice: 1,2,3 or 4                                      .','Warning');
    ch1 = strcmpi(msgin, '1') ;
    ch2 = strcmpi(msgin, '2') ;
    ch3 = strcmpi(msgin, '3') ;
    ch4 = strcmpi(msgin, '4') ;
end

if ch1==1
   flox = flo1;
    freqdev=freqdev1;
else if ch2==1
        flox = flo2;
         freqdev=freqdev2;
else if ch3==1
        flox = flo3;
         freqdev=freqdev3;
    else
    flox = flo4;
     freqdev=freqdev4;
    end
end

end
%%      Shifting to Fif, Bandpass filtering and  Demodulation 
s = u.*cos(2*pi*flox*t);
i = filter(inter,s); %'inter' filter created by filterDesign Tool
                     % passband BW = 65K Hz, filter is not ideal
                   % Low Fpass = 420K Hz & high Fpass = 485K Hz
% y=gradient(i);
% k=find(y<0);
% y(k)=0;
% [b,a]=butter(2,(3500/Fs/2));
% r=filter(b,a,y);
% r=50*(r-mean(r));
% r=100*gradient(r);
r=fmdemod(i,fif,Fs,freqdev);
r=r*2;


sound(r(1:557:length(r)),8000);
%%   Attempt at analog butter filter 
%  
% % % h_band=gaurd/2 -7500;
% % % h=(fc4-h_band)/Fs/2;
% % % l=(fc4+h_band)/Fs/2;
% % % [bb,aa]=butter(5,[h l],'s');
% % % [z p k] = cheb2ap(4,20);
% % % [bp ap] = zp2tf(z,p,k);
% % [bp ap] = butter(10,fc1/(Fs/2)); 
% % [b,a]=lp2bp(bp,ap,(2*pi*fc1),2*pi*2*h_band);%wn=cuttofffreq/(Fs/2)
% % % [bb,aa]=fir1(48,[h l]);
% % freqz(bp,ap);
% % title('bp ap')
% % y4=filter(b,a,x4);
% % 
% % r=fmdemod(u,fc1,Fs,freqdev);
% % sound(r(1:278:length(r)),8000)
% 
% 
%% Modulation and Multiplexing
%%%  PLOT of U %%%
uk=fft(u)/Fs*2;
uk=fftshift(abs(uk));
fx=-Fs/2:Fs/length(uk):((Fs/2)-(Fs/length(uk)));
figure(5)
stem(fx,uk,'b','LineWidth',0.1);
title('Multiplexed Signal')
xlabel('F(Hz)','fontweight','bold','fontsize',14);
ylabel('U(f)','fontweight','bold','fontsize',14);
 xlim([-1000000 1000000]);
grid on;
title('uk')


%% Shifting to Intermediate Frequency
%%%  PLOT of S %%%
sk=fft(s)/Fs*2;
sk=fftshift(abs(sk));
fx=-Fs/2:Fs/length(sk):((Fs/2)-(Fs/length(sk)));
figure(6)
stem(fx,sk,'b','LineWidth',0.1);
title('Output of IF mixer')
xlabel('F(Hz)','fontweight','bold','fontsize',14);
ylabel('S(f)','fontweight','bold','fontsize',14);
% xlim([28000 50000]);
grid on;
title('sk')

%%  Bandpass Filtering at Fif
%%%  PLOT of I %%%
ik=fft(i)/Fs*2;
ik=fftshift(abs(ik));
fx=-Fs/2:Fs/length(ik):((Fs/2)-(Fs/length(ik)));
figure(7)
stem(fx,ik,'b','LineWidth',0.1);
title('Output of bandpass centered at Fif')
xlabel('F(Hz)','fontweight','bold','fontsize',14);
ylabel('I(f)','fontweight','bold','fontsize',14);
 xlim([-750000 750000]);
grid on;
title('ik')

%%  
%%%  PLOT of M %%%
rk=fft(r)/Fs*2;
rk=fftshift(abs(rk));
fx=-Fs/2:Fs/length(rk):((Fs/2)-(Fs/length(rk)));
figure(8)
stem(fx,rk,'b','LineWidth',0.1);
title('Message signal Spectrum')
xlabel('F(Hz)','fontweight','bold','fontsize',14);
ylabel('M(f)','fontweight','bold','fontsize',14);
 xlim([-4000 4000]);
grid on;
title('rk')
