function [ periodogram_ ] = welch_method( noise,N,K )
%Uses welch method for periodogram estimation
 noise_=buffer(noise,N/K,ceil(N/(K*2)));
 noise_=noise_(:,2:end);
 periodogram_=0;
 
    for n=1:K*2-1
        
        periodogram=abs(fft(noise_(:,n),(N/K))).^2/(N/K);
        periodogram(1)=periodogram(2);
        periodogram_=periodogram+periodogram_;
    end
  periodogram_=periodogram_/(K*2-1);
end

