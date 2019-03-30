clear all; close all; clc

load handel
v = y'/2;

t = (1:length(v))/Fs;
L = 9;
n=length(v);
k=(2*pi/L)*[0:n/2-1 -n/2:-1];
ks=fftshift(k);
%%
v = v(1:73112);
t = t(1:73112);

width=[50 0.5];
tslide = 0:1:L;

for i = 1:length(width)
    Sgt_spec = [];
    Sst_spec = [];
    Smt_spec = [];
    figure(2*i - 1)
    for j=1:length(tslide)
        
        g = exp(-width(i)*(t - tslide(j)).^2); %gaussian filter
        s = (abs(t-tslide(j)) <= 1/width(i)); %shannon filter
        m = 2.*(1 - ((t-tslide(j))/width(i).^-1).^2)...
            .*exp(-((t-tslide(j)).^2)/...
            (2.*width(i).^-2))/(sqrt(3.*width(i).^-1)...
            .*pi^(1/4));%Ricker hat filter 
        
        Sg = g.*v; %filtered with gaussian
        Ss = s.*v; %filtered with shannon
        Sm = m.*v; %mexican hat filter
        
        Sgt = fft(Sg); %fft gaussian
        Sst = fft(Ss); %fft shannon
        Smt = fft(Sm); %fft mexican
        
        subplot(3,3,1:3), 
        plot(t,v,'k',t,g,'r',t,s,'g',t,m,'b')
        xlabel('Time [sec]');
        ylabel('Amplitude');
        title('Signal amd filters');
        
        subplot(3,3,4), 
        plot(t,Sg,'k')
        xlabel('Time [sec]');
        ylabel('Amplitude');
        title('Gaussian Filtered signal');
        
        subplot(3,3,5), 
        plot(t,Ss,'k')
        xlabel('Time [sec]');
        ylabel('Amplitude');
        title('Shannon Filtered signal');
        
        subplot(3,3,6), 
        plot(t,Sm,'k')
        xlabel('Time [sec]');
        ylabel('Amplitude');
        title('Ricker Filtered signal');
        
        subplot(3,3,7), 
        plot(ks,abs(fftshift(Sgt))/max(abs(Sgt)))
        xlabel('Frequency');
        ylabel('Amplitude');
        title('Gaussian Fourier signal');
        
        subplot(3,3,8),
        plot(ks,abs(fftshift(Sst))/max(abs(Sst)))
        xlabel('Frequency');
        ylabel('Amplitude');
        title('Shannon Fourier signal');
        
        subplot(3,3,9),
        plot(ks,abs(fftshift(Smt))/max(abs(Smt)))
        xlabel('Frequency');
        ylabel('Amplitude');
        title('Ricker Fourier signal');
        
        drawnow
        
        Sgt_spec = [Sgt_spec; abs(fftshift(Sgt))];
        Sst_spec = [Sst_spec; abs(fftshift(Sst))];
        Smt_spec = [Smt_spec; abs(fftshift(Smt))];
    end
    
    figure(2*i)
    subplot(3,1,1)
    pcolor(tslide,ks,Sgt_spec.'), ...
        shading interp, colormap(hot)
    str = sprintf('Gaussian Spectrogram with width parameter %f',...
        width(i));
    title(str)
    xlabel('Time Slide')
    ylabel('Frequency')
    
    subplot(3,1,2)
    pcolor(tslide,ks,Sst_spec.'), shading interp,...
        colormap(hot)
    str = sprintf('Shannon Spectrogram with width parameter %f',...
        width(i));
    title(str)
    xlabel('Time Slide')
    ylabel('Frequency')
    
    subplot(3,1,3)
    pcolor(tslide,ks,Smt_spec.'), shading interp,...
        colormap(hot)
    str = sprintf('Ricker Spectrogram with width parameter %f',...
        width(i));
    title(str)
    xlabel('Time Slide')
    ylabel('Frequency')
end