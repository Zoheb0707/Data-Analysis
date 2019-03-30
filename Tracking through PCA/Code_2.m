clear all; close all; clc;

%loading data

frames_1 = load('Data/cam1_2.mat');
frames_1 = frames_1.('vidFrames1_2');

frames_2 = load('Data/cam2_2.mat');
frames_2 = frames_2.('vidFrames2_2');

frames_3 = load('Data/cam3_2.mat');
frames_3 = frames_3.('vidFrames3_2');

%%
%making frames uniform
min_frames = 314;

frames_1 = frames_1(1:480,1:640, 1:3, 1:min_frames);
frames_2 = frames_2(1:480,1:640, 1:3, 1:min_frames);
frames_3 = frames_3(1:480,1:640, 1:3, 1:min_frames);

%%
%converting to rgb
%extracting the location of the can
M = zeros(6,min_frames);

for i = 1:min_frames
    
    frame_1_to_add = imresize(rgb2gray(frames_1(1:480,1:640, 1:3, i)...
        ), [120,160]);
    frame_2_to_add = imresize(rgb2gray(frames_2(1:480,1:640, 1:3, i)...
        ), [120,160]);
    frame_3_to_add = imresize(rgb2gray(frames_3(1:480,1:640, 1:3, i)...
        ), [120,160]);
    
    frame_1_to_add = double(frame_1_to_add);
    frame_2_to_add = double(frame_2_to_add);
    frame_3_to_add = double(frame_3_to_add);
    
    x_vals_1 = 0;
    y_vals_1 = 0;
    count_1 = 0;
    
    x_vals_2 = 0;
    y_vals_2 = 0;
    count_2 = 0;
    
    x_vals_3 = 0;
    y_vals_3 = 0;
    count_3 = 0;
    
    for j = 1:120
        for k = 1:160
            if frame_1_to_add(j,k) >= 240 && ...
                    k >= 70 && k <= 100
                x_vals_1 = x_vals_1 + j;
                y_vals_1 = y_vals_1 + k;
                count_1 = count_1 + 1;
            end
            if frame_2_to_add(j,k) >= 240 && ...
                    k >= 60 && k <= 90
                x_vals_2 = x_vals_2 + j;
                y_vals_2 = y_vals_2 + k;
                count_2 = count_2 + 1;
            end
            if frame_3_to_add(j,k) >= 240 && ...
                    j >= 60 && j <= 80 && k >= 60 && k <= 131
                x_vals_3 = x_vals_3 + j;
                y_vals_3 = y_vals_3 + k;
                count_3 = count_3 + 1;
            end
        end
    end
    
    x_vals_1 = x_vals_1/count_1;
    y_vals_1 = y_vals_1/count_1;
    
    x_vals_2 = x_vals_2/count_2;
    y_vals_2 = y_vals_2/count_2;
    
    x_vals_3 = x_vals_3/count_3;
    y_vals_3 = y_vals_3/count_3;
    
    M(1,i) = x_vals_1;
    M(2,i) = y_vals_1;
    
    M(3,i) = x_vals_2;
    M(4,i) = y_vals_2;
    
    M(5,i) = x_vals_3;
    M(6,i) = y_vals_3;
end
%%
%De meaning data
M(isnan(M))=0;
mean_data = zeros(6);
for i = 1:6
    mean_data(i) = mean(M(i,:));
    M(i,:) = M(i,:) - mean(M(i,:));
    
end
%%
%SVD
[U,S,V] = svd(M,'econ');
disp(diag(S).');
%%
%projection
proj = U.'*M;
%%
%plotting projection
figure(1)
plot(proj(1,1:300))
title('Projection using PCA')
xlabel('Frames')
ylabel('Displacement')
hold on
plot(proj(2,1:300))
hold on
plot(proj(3,1:300))
hold on
plot(proj(4,1:300))
legend('1','2','3','4');