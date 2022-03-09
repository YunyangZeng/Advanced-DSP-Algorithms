% Author: Yunyang Zeng
clear all;
close all;
%% Read audio
[welc,Fs]=audioread("welcome16k.wav");
welc_seg=welc(17000:17319);
%% a) Autocorrelation method
window=hamming(length(welc_seg));
windowed_welc=window.*welc_seg;
[autocorrelation,lg] = xcorr(windowed_welc,'biased');
autocorrelation(lg<0) = [];
autocorrelation=autocorrelation(1:15);
P=10;
[a,e,K]=levinson792(autocorrelation,P);
G=e(end)^0.5;
figure();
freqz(G,[1 -1*a],100);
title("Autocorrelation method, number of poles = "+P);
figure();
zplane(G,[1 -1*a]);
title("Z plane, Autocorrelation method");
fprintf("Autocorrelation method, a=")
a
fprintf("Autocorrelation method, k=")
K(end)
fprintf("Autocorrelation method, G=")
G
%% b) Covariance method
ac=autocorrelation;
N=length(welc_seg);
P=10;
phi=zeros(P+1,P+1);
y=0;
for i=0:P
    for k=0:P        
        for m=-i:N-i-1
            s_m=welc(17000+m);
            s_m_i_k=welc(17000+m+i-k);
            phi(i+1,k+1)=phi(i+1,k+1)+s_m*s_m_i_k;
            
        end
        y=y+1;
    end
end
R=phi(2:end,2:end);
P_=phi(2:end,1);

alpha=inv(R)*P_;
K_=alpha(end);
G_=(ac(1)-alpha(1:end)'*ac(2:P+1))^0.5;
figure();
freqz(G_,[1 -1*alpha'],100);
title("Covariance method, number of poles = "+P);
figure();
zplane(G_,[1 -1*alpha']);
title("Z plane, Covariance method");
fprintf("Covariance method, a=")
alpha'
fprintf("Covariance method, k=")
K_
fprintf("Covariance method, G=")
G_
fprintf("The Covariance method provides less mean-squared prediction error but is less computaionally efficient than the autocorrelation method")
