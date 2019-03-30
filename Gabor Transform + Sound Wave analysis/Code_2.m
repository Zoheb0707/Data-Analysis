clear all; close all; clc

L=16; % record time in seconds
y=audioread('music1.wav'); 
Fs=length(y)/L;
v = y'/2;
t = (1:length(v))/Fs;
n=length(v);
k=(2*pi/L)*[0:n/2-1 -n/2:-1]; 
ks=fftshift(k);

width=[100];
tslide = 0.5:0.53:L-1;

for i = 1:length(width)
    figure(2*i - 1)
    title('Piano')
    
    Sgt_spec = [];
    frequencies = [];
    for j=1:length(tslide)
        g = exp(-width(i)*(t - tslide(j)).^2);
        Sg = g.*v; %filtered with gaussian
        Sgt = fft(Sg); %fft gaussian
        
        %signal and filter plot
        subplot(3,1,1)
        plot(t,v,'k',t,g,'r')
        xlabel('Time [sec]');
        ylabel('Amplitude');
        title('Signal and filters (Piano)');
        
        subplot(3,1,2), 
        plot(t,Sg,'k')
        xlabel('Time [sec]');
        ylabel('Amplitude');
        title('Gaussian Filtered signal');
        
        subplot(3,1,3), 
        plot(ks,abs(fftshift(Sgt))/max(abs(Sgt)))
        xlabel('Frequency');
        ylabel('Amplitude');
        title('Gaussian Fourier signal');
        
        drawnow
        
        [val, index] = max(abs(fftshift(Sgt)));
        frequency = ks(index)/(2*pi);
        frequencies = [frequencies, frequency];
        
        Sgt_spec = [Sgt_spec; abs(fftshift(Sgt))];
    end
    frequencies = abs(frequencies);
    %disp(frequencies);
    figure(2*i)
    
    subplot(2,1,1)
    pcolor(tslide,ks,Sgt_spec.'), ...
        shading interp, colormap(hot)
    %axis([0 length(tslide) 300000 400000])
    str = sprintf('Gaussian Spectrogram with width parameter %f',...
        width(i));
    title(str)
    xlabel('Time slide')
    ylabel('Frequency')
    
    subplot(2,1,2)
    scatter(tslide,frequencies);
    xlabel('Time');
    ylabel('Frequency (Hz)');
    title('Frequency of notes played wrt time (Piano)');
    axis([0 L 100 400])
    drawnow   
end

disp(frequencies(1:length(tslide)));

%%
%clear all; close all; clc

L=14; % record time in seconds
y=audioread('music2.wav'); 
Fs=length(y)/L;
v = y'/2;
t = (1:length(v))/Fs;
n=length(v);
k=(2*pi/L)*[0:n/2-1 -n/2:-1]; 
ks=fftshift(k);

width=[100];
tslide = 0.2:0.53:L;

for i = 1:length(width)
    figure(4*i - 1)
    title('Recorder')
    
    Sgt_spec = [];
    frequencies = [];
    for j=1:length(tslide)
        g = exp(-width(i)*(t - tslide(j)).^2);
        Sg = g.*v; %filtered with gaussian
        Sgt = fft(Sg); %fft gaussian
        
        %signal and filter plot
        subplot(3,1,1)
        plot(t,v,'k',t,g,'r')
        xlabel('Time [sec]');
        ylabel('Amplitude');
        title('Signal and filters (Recorder)');
        
        subplot(3,1,2), 
        plot(t,Sg,'k')
        xlabel('Time [sec]');
        ylabel('Amplitude');
        title('Gaussian Filtered signal');
        
        subplot(3,1,3), 
        plot(ks,abs(fftshift(Sgt))/max(abs(Sgt)))
        xlabel('Frequency');
        ylabel('Amplitude');
        title('Gaussian Fourier signal');
        
        drawnow
        
        [val, index] = max(abs(fftshift(Sgt)));
        frequency = ks(index)/(2*pi);
        frequencies = [frequencies, frequency];
        
        Sgt_spec = [Sgt_spec; abs(fftshift(Sgt))];
    end
    frequencies = abs(frequencies);
    %disp(frequencies);
    figure(4*i)
    
    subplot(2,1,1)
    pcolor(tslide,ks,Sgt_spec.'), ...
        shading interp, colormap(hot)
    %axis([0 length(tslide) 300000 400000])
    str = sprintf('Gaussian Spectrogram with width parameter %f',...
        width(i));
    title(str)
    xlabel('Time slide')
    ylabel('Frequency')
    
    subplot(2,1,2)
    scatter(tslide,frequencies);
    xlabel('Time');
    ylabel('Frequency (Hz)');
    title('Frequency of notes played wrt time (Recorder)');
    axis([0 L 700 1200])
    drawnow   
end

disp(frequencies(2:length(tslide)));