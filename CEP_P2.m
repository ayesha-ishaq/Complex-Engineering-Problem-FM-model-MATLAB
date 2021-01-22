%% Reading Text file and binary representation 
fidi = fopen('Text.txt','r');  %open file and get File ID
data=fread(fidi);              %Read data in double format
fclose(fidi);
cdata=char(data);              %Convert data into char (checking if data is read correctly)
bdata1=dec2bin(data,8)-'0';     %Convert double to binary(double format)
%check=bi2de(bdata1,'left-msb');
bdata=transpose(bdata1);        %Transpose before reshape
x=reshape(bdata,1,length(cdata)*8);     %Vector of binary values


%% Transimitter 
f =10000000; % f in Hz 10MHz
Fs =2.2*f; % samples per symbol 22MHz
Tb =0.000001; % Bit duration
Ab =1; % Bit Amplitude.
%Code to define time t that contains N samples for Tb second(s)
Ts=1/Fs;
N=length(x)*Tb*Fs;
n=0:N-1;   %Sampling Index 
t=n*Ts;    %time vector
w = sqrt(2/Tb)*cos(2*pi*f*t); % Defining the Basis Function/ Carrier Signal

%Code to generate and plot BPSK_Pulse
u=ones(1,length(t));
Fss=Fs*Tb;

for i=1:length(x)
    if x(i) == 0
 
        if i ==1
            u(1:1:Fss)= u(1:1:Fss)*(-1);
        else if i==length(x)
         u(( (i-1)*Fss+1):1:i*Fss)=u(( (i-1)*Fss+1):1:i*Fss).*(-1);       
        else
        u(( (i-1)*Fss+1):1:i*Fss)=u(( (i-1)*Fss+1):1:i*Fss).*(-1);
        
    end
        end  
    end
end

tt=t(1:1:length(t)/length(data));
figure(1) %Binary PAM Plot
plot(tt,u(1:1:length(tt)),'b','LineWidth',1)
grid on;
ylim( [-1.5 1.5])
title('Binary PAM')
xlabel('Time (sec)','fontweight','bold','fontsize',14);
ylabel('Amplitude (V)','fontweight','bold','fontsize',14);

u=u.*w;   %BPSK

figure(2)
plot(tt,u(1:1:length(tt)),'b','LineWidth',1)
grid on;
xlim( [0 3/1000000])

title('BPSK')
xlabel('Time (sec)','fontweight','bold','fontsize',14);
ylabel('Amplitude (V)','fontweight','bold','fontsize',14);

%% Noise

variance = input('Enter vareince:');

r=u+sqrt(variance)*randn(1,length(u)); % Addition of Noise in the Channel or use awgn (r, SNR) command

figure(3)
plot(tt,r(1:1:length(tt)),'b','LineWidth',1)
grid on;
xlim( [0 3/1000000])

title('BPSK with noise')
xlabel('Time (sec)','fontweight','bold','fontsize',14);
ylabel('Amplitude (V)','fontweight','bold','fontsize',14);
%% Correlator Receiver and Detection Block

rec=0;
for i=1:length(x)
    
 
        if i ==1
            ss(1:1:Fss)= r(1:1:Fss).*w(1:1:Fss);
            s(i)=sum(ss(1:1:Fss));
        else if i==length(x)
         ss(( (i-1)*Fss+1):1:i*Fss)=r(( (i-1)*Fss+1):1:i*Fss).*w(( (i-1)*Fss+1):1:i*Fss); 
         s(i)=sum(ss(( (i-1)*Fss+1):1:i*Fss));
        else
        ss(( (i-1)*Fss+1):1:i*Fss)=r(( (i-1)*Fss+1):1:i*Fss).*w(( (i-1)*Fss+1):1:i*Fss);
        s(i)=sum(ss(( (i-1)*Fss+1):1:i*Fss));
        
        
    end
        end
        
        if s(i)>0
            rec(i)=1;
        else
            rec(i)=0;
     
        end
end



%% Error checking 

err=0;
index=1;
for i=1:length(rec)
    
    if x(i) ~= rec(i)
        err(index)=i;
        index = index+1;
    end
end
    
    if length(err)==1& err(1)==0
        data_bits = transpose(x);
        Rcvd_Bits = transpose(rec);
        disp('Congragulations! Data received correctly')
        
    else 
        
        data_bits = transpose(x);
        Rcvd_Bits = transpose(rec);
        
       
        disp('Warning! errors detected.')
        %Total_Bits = length(
        Number_of_Errors = length(err)
        
        Position = err
    end
    
    %% Convert binary back to char
y = reshape(rec,8,length(cdata));
y = transpose(y);
y = bi2de(y,'left-msb');
 out = char(y);
 disp(out')
%% Write to Output Text File
fid2 = fopen('Out.txt','w');  %open file and get File ID
fprintf(fid2,'%c',out);
%print('%c',transpose(out));
%output = transpose(out)
fclose(fid2);

