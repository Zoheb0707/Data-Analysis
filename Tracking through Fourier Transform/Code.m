clear all;close all; clc;

load Testdata

L=15; % spatial domain
n=64; % Fourier modes
t2=linspace(-L,L,n+1); 
t=t2(1:n);

x2=linspace(-L,L,n+1); 
x=x2(1:n); 
y=x; 
z=x;

k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; 
ks=fftshift(k);

[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

%%
%1)
Utn_avg = zeros(1,n);
for j=1:20
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    Utn(:,:,:) = fftn(Un(:,:,:));
    
    Utn_avg = Utn_avg + Utn(:,:,:);
end

figure(1)
close all, isosurface(X,Y,Z,abs(Un),0.4)
axis([-20 20 -20 20 -20 20]), grid on, drawnow
title('Noisy Data')

Utn_avg_p=fftshift(Utn_avg)/20;
to_plot = abs(Utn_avg_p)/max(max(max(abs(Utn_avg_p))));
figure(2)
isosurface(Kx,Ky,Kz,to_plot,0.6)
axis([-7 7 -7 7 -7 7]), grid on, drawnow
title('Avg Frequency in the Fourier Domain')

[mxv,idx] = max(to_plot(:));
[r,c,p] = ind2sub(size(to_plot),idx); %3d indices
Frequency_signal = [Kx(r,c,p), Ky(r,c,p), Kz(r,c,p)];
disp(Frequency_signal);

%%
%2)
filter = exp(-0.2*( (Kx - Kx(r,c,p)).^2 + ...
(Ky - Ky(r,c,p)).^2+ (Kz - Kz(r,c,p)).^2));

figure(3)
isosurface(Kx,Ky,Kz,abs(filter),0.6)
axis([-20 20 -20 20 -20 20]), grid on, drawnow
title('Gaussian Filter')

coords = zeros(3,20);
for j=1:20
    Un(:,:,:)=reshape(Undata(j,:),n,n,n);
    Utn(:,:,:) = fftn(Un(:,:,:));
    Utn(:,:,:) = fftshift(Utn(:,:,:));
    Utfn(:,:,:) = filter.*Utn(:,:,:);
    ufn(:,:,:) = ifftshift(Utfn(:,:,:));
    ufn(:,:,:) = ifftn(ufn(:,:,:));
    abs_ufn = abs(ufn);
    [mxv,idx] = max(abs_ufn(:));
    [r,c,p] = ind2sub(size(abs_ufn),idx);
    coords(:,j) = [X(r,c,p),Y(r,c,p),Z(r,c,p)];
end
figure(4)
plot3(coords(1,:),coords(2,:), coords(3,:));
axis([-12 12 -12 12 -12 12]), grid on, drawnow
title('Path of the marble through time')
xlabel('X')
ylabel('Y')
zlabel('Z')
%%
%3)
disp(coords(:,20));