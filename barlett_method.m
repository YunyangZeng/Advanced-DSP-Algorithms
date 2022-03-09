function [ periodogram_ ] = barlett_method( noise,N,K )
% Uses barlett method for periodogram estimation
    noise_=buffer(noise,N/K);
    periodogram_=0;
    for n=1:K
        %sample_mean=mean(noise_(n,:));
        
        periodogram=abs(fft(noise_(:,n),(N/K))).^2/(N/K);
        periodogram(1)=periodogram(2);
        periodogram_=periodogram+periodogram_;
    end
    periodogram_=periodogram_/K;

end

