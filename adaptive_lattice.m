% Author: Yunyang Zeng
close all;
clear;
%% Read sound files
[Carspeech_pri,fs]=audioread("Carspeech_pri.wav");
[Carspeech_ref,fs]=audioread("Carspeech_ref.wav");
%% Solve Wiener Hopf Equation
Carspeech_pri=Carspeech_pri(85:end); % by checking the crosscorelation sequence, the peak appears at the 85th sample, so delay the pri channel signal by 85 samples to align with the ref.
Carspeech_ref=Carspeech_ref(1:end-84); % make the pri and ref the same length.
n=20; % number of adaptive filter coefficients.
autocor=xcorr(Carspeech_ref','unbiased',150); % get the autocorrelation sequence of the ref channel samples.
autocor=autocor(151:end); % only need half of the symetric autocorrelation sequence.
R=toeplitz(autocor(1:n)); % make autocorrelation sequence a toeplitz matrix to solve the Wiener Hopf equation
[V,D]=eig(R); % find eigenvalue of the autocorrelation matrix in order to calulate the range of mu
crosscor=xcorr(Carspeech_ref',Carspeech_pri','unbiased',150);% get the crosscorrelation sequence between the pri and ref channel samples.
crosscor=flip(crosscor(1:151));% P=E(dk*Xk)=E(pri[k]*REF[k]). The crosscorelation sequence is not symmetric, we need the left half part reversed.
P=crosscor(1:n)';
Ws=inv(R)*P;% Solve the Wiener Hopf Equation
%% Algorithm
mu=0.01*1/sum(D,"all");%set mu value
e_spread=max(D,[],"all")/(min(D(D>0),[],"all")+1e-16);%calculate eigenvalue spread
fprintf("eigen value spread:")
e_spread

Wk=ones(n,1);%initialize adaptive filter coefficients
k=ones(n-1,1)*0.1;%initialize lattice filter coefficients
f=ones(n,2);%initialize lattice filter f values
b=ones(n,2);%initialize lattice filter b values
LMSe=[];
mse=[];
tic
for i=n:length(Carspeech_ref)*2
    ind=mod(i,length(Carspeech_ref));
    if ind>=1
        xk=Carspeech_ref(ind);
        
    %Implementation of lattice filter
    if ind==1
        f(:,1)=f(:,1)*0;
        b(:,1)=b(:,1)*0;
        f(:,2)=f(:,2)*xk;
        b(1,2)=xk;
        b(2:end,2)=k*xk;
    else
        f_k=zeros(n,1);
        b_k=zeros(n,1);
        for j=1:n
            if j==1
                f_k(j)=xk;
                b_k(j)=xk;
            else
                f_k(j)=f_k(j-1)+b(j-1,2)*k(j-1);
                b_k(j)=k(j-1)*f_k(j-1)+b(j-1,2);
            end
        end
        f(:,1)=f(:,2);
        f(:,2)=f_k;
        b(:,1)=b(:,2);
        b(:,2)=b_k;
    end
    
    %adaptive filter
    dk=Carspeech_pri(ind);
    ek=dk-b_k'*Wk;
    Wkp1=Wk+2*mu*b_k*ek;
    Wk=Wkp1;
    
    % update lattice filter coefficients ks
    for g=1:(n-1)
        k_(g)=k(g)-2*mu*f_k(g+1)*b(g,1);
    end
    k=k_;
    
    % obtain the enhanced signal in the second iteration
    if (i/length(Carspeech_ref))>1
        LMSe(ind)=ek;
    end
        
    end 
    
end
toc
%% Compare Implementation with matlab dsp.AdaptiveLatticeFilter() function 
alf = dsp.AdaptiveLatticeFilter(20);
[y,err]=alf(Carspeech_ref,Carspeech_pri);
w=alf.Coefficients;
w=w';

fprintf("adaptive lattice filter coefficients obtained using my implentation");
Wk
fprintf("adaptive lattice filter coefficients obtained using matlab dsp.AdaptiveLatticeFilter() function ");
w

%% Play and save the resynthesized speech and plot
audiowrite("yunyangz_C10_4res.wav",LMSe,fs);
figure();
plot(LMSe);
title("waveform after noise cancellation");
figure();
plot(Carspeech_pri);
title("waveform before noise cancellation");

