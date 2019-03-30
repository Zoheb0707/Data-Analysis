clear all; close all; clc
%% Loading video from reader
v = VideoReader('Data/test_1.mp4');
min_frames = v.NumberOFFrames;
v = VideoReader('Data/test_1.mp4');
frames_1 = read(v);
%% Converting to gray scale
frames_rgb = zeros(size(frames_1, 1)*...
    size(frames_1, 2),min_frames);
for i = 1:min_frames
    frames_rgb(:,i) = reshape(rgb2gray(frames_1(:,:,:,i)), ...
        size(frames_1, 1)*size(frames_1, 2), 1);
end
%% Creating DMD matrices
X = frames_rgb;
dt = 1;
t = 0:dt:min_frames;

X1 = X(:,1:end-1);
X2 = X(:,2:end);
%% SVD
[U,S,V] = svd(X1, 'econ');
%% Finding a good r to use
figure()
plot(diag(S), 'o');
title('Trend in singular values')
xlabel('Singular Value Index')
ylabel('Singular Values')
%% choosing r value based on the plot
r = 10;
%% Doing rank reduction
Ur = U(:,1:r);
Sr = S(1:r,1:r);
Vr = V(:,1:r);
%% Finding A_tilde and DMD modes
Atilde = Ur'*X2*Vr/Sr;
[W,D] = eig(Atilde);
Phi = X2*Vr/Sr*W;
%% DMD spectra
lambda = diag(D);
omega = log(lambda)/dt;

omega_to_use = min(abs(omega)); 
%% %% find DMD solution
b = Phi\X(:,1);
time_dynamics = zeros(r,min_frames);
for iter = 1:min_frames
    time_dynamics(:,iter) = (b.*exp(omega_to_use*t(iter)));
end
X_dmd = Phi*time_dynamics;
%% Finding sparse with residual
X_sparse = X - abs(X_dmd);
%% Finding residual
Residual = X_sparse .* (X_sparse < 0);

%% removing residual
X_dmd = abs(X_dmd);% + Residual;
X_sparse = X_sparse - Residual;
%% plot X_dmd
figure()
for i = 40:41
    subplot(1,3,1)
    to_show = X(:,i);
    to_show = reshape(to_show,[size(frames_1, 1), size(frames_1, 2)]);
    pcolor(flipud(abs(to_show))), shading interp, colormap gray;
    title('Original Video')
    
    subplot(1,3,2)
    to_show = X_dmd(:,i);
    to_show = reshape(to_show,[size(frames_1, 1), size(frames_1, 2)]);
    pcolor(flipud(abs(to_show))), shading interp, colormap gray;
    title('Background')
    
    subplot(1,3,3)
    to_show = X_sparse(:,i);
    to_show = reshape(to_show,[size(frames_1, 1), size(frames_1, 2)]);
    pcolor(flipud(abs(to_show))), shading interp, colormap gray;
    title('Foreground')
    drawnow
end