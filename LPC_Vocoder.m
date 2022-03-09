%Author: Yunyang Zeng

%% a) LPC vocoder
%% Read audio
clear 
close all
[welc,Fs]=audioread("welcome16k.wav");
welc=welc(1:end-mod(length(welc),160));
N=320;
window=hamming(N);
a_list=[];
g_list=[];
%white_noise=wgn(N,1,0);
seq_odd=[];
resyn_even=[];
resyn_odd=[];
seq_even=[];
P=14;

for i=0:(length(welc)-N)/160
    seg=welc(160*i+1:160*i+N);
    wind_seg=window.*seg;
    
    if mod(i,2)==0
        seq_even=[seq_even wind_seg];
    else
        seq_odd=[seq_odd wind_seg];
    end
        
    %{  
    [a,g]=lpc(wind_seg,14);
    a_list=[a_list;a];
    g_list=[g_list; g]; 
    
    
    if i==0
        [resyn,zi]=filter(g,a,white_noise,zeros(length(a)-1,1));
    else    
        [resyn,zi]=filter(g,a,white_noise,zi);
    end
    %}
end
%% Resynthesize with white noise exitation
for j=1:64
    sequence=seq_even(:,j);
    [a,g]=lpc(sequence,P);
    
    white_noise=wgn(N,1,0);
    phi=zeros(P+1,1);
    for m=0:P
        phi(m+1)=dot(sequence(1:N-m),sequence(1+m:N)); 
    end
    
    G=sqrt(a*phi);
    if j==1
        [resyn,zi]=filter(G,a,white_noise,zeros(length(a)-1,1));
    else    
        [resyn,zi]=filter(G,a,white_noise,zi);
    end
    resyn_even=[resyn_even; resyn];
    
    
end

for k=1:63
    sequence=seq_odd(:,k);
    [a,g]=lpc(sequence,P);
    
    white_noise=wgn(N,1,0);
    phi=zeros(P+1,1);
    for m=0:P
        phi(m+1)=dot(sequence(1:N-m),sequence(1+m:N)); 
    end
    
    G=sqrt(a*phi);
    if k==1
        [resyn,zi]=filter(G,a,white_noise,zeros(length(a)-1,1));
    else    
        [resyn,zi]=filter(G,a,white_noise,zi);
    end
    resyn_odd=[resyn_odd; resyn];
    
    
end

resyn_white_noise=[resyn_even(1:160); 0.5.*(resyn_even(161:end-160)+resyn_odd); resyn_even(end-160+1: end)];
resyn_white_noise=resyn_white_noise./max(abs(resyn_white_noise));
audiowrite("resyn_white_noise.wav",resyn_white_noise,16000);
%% Resynthesize with periodic pulse train exitation    
f=100;
fs=16000;
impulse_train=zeros(N,1);
impulse_train(1:1/f*fs:end)=1;
resyn_even_impulse=[];
resyn_odd_impulse=[];
for j=1:64
    sequence=seq_even(:,j);
    [a,g]=lpc(sequence,P);  
    
    phi=zeros(P+1,1);
    for m=0:P
        phi(m+1)=dot(sequence(1:N-m),sequence(1+m:N)); 
    end
    
    G=sqrt(a*phi);
    if j==1
        [resyn,zi]=filter(G,a,impulse_train,zeros(length(a)-1,1));
    else    
        [resyn,zi]=filter(G,a,impulse_train,zi);
    end
    resyn_even_impulse=[resyn_even_impulse; resyn];
    
    
end

for k=1:63
    sequence=seq_odd(:,k);
    [a,g]=lpc(sequence,P);    
    phi=zeros(P+1,1);
    for m=0:P
        phi(m+1)=dot(sequence(1:N-m),sequence(1+m:N)); 
    end
    
    G=sqrt(a*phi);
    if k==1
        [resyn,zi]=filter(G,a,impulse_train,zeros(length(a)-1,1));
    else    
        [resyn,zi]=filter(G,a,impulse_train,zi);
    end
    resyn_odd_impulse=[resyn_odd_impulse; resyn];
    
    
end

resyn_impulse=[resyn_even_impulse(1:160); 0.5.*(resyn_even_impulse(161:end-160)+resyn_odd_impulse); resyn_even_impulse(end-160+1: end)];
resyn_impulse=resyn_impulse./max(abs(resyn_impulse));
audiowrite("resyn_impulse.wav",resyn_impulse,16000);
figure();
colormap jet
spectrogram(welc,hamming(320),160,512,16000,"yaxis");
title("Original signal");
figure();
colormap jet
spectrogram(resyn_impulse,hamming(320),160,512,16000,"yaxis");
title("Resynthesized signal with impulse train excitation");
figure();
colormap jet
spectrogram(resyn_white_noise,hamming(320),160,512,16000,"yaxis");
title("Resynthesized signal with white noise excitation");



