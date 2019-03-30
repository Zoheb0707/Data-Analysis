clear all; close all; clc;
%% Loading songs

%Fs can be found by doing [song, Fs] = audioread();
Fs = 44100;

TS_1 = audioread("Data/Music_Part_1\TS1.mp3");
TS_2 = audioread("Data/Music_Part_1\TS2.mp3");
TS_3 = audioread("Data/Music_Part_1\TS3.mp3");
KL_1 = audioread("Data/Music_Part_1\KL1.mp3");
KL_2 = audioread("Data/Music_Part_1\KL2.mp3");
KL_3 = audioread("Data/Music_Part_1\KL3.mp3");
CP_1 = audioread("Data/Music_Part_1\CP1.mp3");
CP_2 = audioread("Data/Music_Part_1\CP2.mp3");
CP_3 = audioread("Data/Music_Part_1\CP3.mp3");

%% Taking avg of left and right stereo and transposing
TS_1 = ((TS_1(:,1) + TS_1(:,2))/2).';
TS_2 = ((TS_2(:,1) + TS_2(:,2))/2).';
TS_3 = ((TS_3(:,1) + TS_3(:,2))/2).';
KL_1 = ((KL_1(:,1) + KL_1(:,2))/2).';
KL_2 = ((KL_2(:,1) + KL_2(:,2))/2).';
KL_3 = ((KL_3(:,1) + KL_3(:,2))/2).';
CP_1 = ((CP_1(:,1) + CP_1(:,2))/2).';
CP_2 = ((CP_2(:,1) + CP_2(:,2))/2).';
CP_3 = ((CP_3(:,1) + CP_3(:,2))/2).';

%% Finding smallest length to make all songs uniform
min_length_1 = min(size(TS_1, 2), size(TS_2, 2));
min_length_2 = min(size(TS_3, 2), size(KL_1, 2));
min_length_3 = min(size(KL_2, 2), size(KL_3, 2));
min_length_4 = min(size(CP_1, 2), size(CP_2, 2));
min_length_5 = min(min_length_1, size(CP_3, 2));
min_length_6 = min(min_length_2, min_length_3);
min_length_7 = min(min_length_1, min_length_4);

min_length = min(min_length_6, min_length_7);
%% Creating song matrix
%removing 5 seconds from start and end of each song 
%as it is mostly silence
%we remove a bit more than that to make segments uniform
songs = [TS_1(1, 178657:min_length - 178657 - 7)
    TS_2(1, 178657:min_length - 178657 - 7)
    TS_3(1, 178657:min_length - 178657 - 7)
    KL_1(1, 178657:min_length - 178657 - 7)
    KL_2(1, 178657:min_length - 178657 - 7)
    KL_3(1, 178657:min_length - 178657 - 7)
    CP_1(1, 178657:min_length - 178657 - 7)
    CP_2(1, 178657:min_length - 178657 - 7)
    CP_3(1, 178657:min_length - 178657 - 7)];
%% Dividing the songs 40 into snippets of 4 seconds each
%for each song put segment 5 to 34 in train
% Putting 9*30 into train
%for each song put segment 35 to 37 in test
% and 9*3 into test
%thus, we have 270 train and 27 test points
train_snippets = zeros(30*9,178656);
test_snippets = zeros(27,178656);
for i = 1:9
    for j = 5:34
        train_snippets((i-1)*30 + j-4,:) = ...
            songs(i,(j-1)*178657 + 1:j*178657 - 1);
    end
end
for j = 35:37
    for k = 1:3
        for i = 1:9
                test_snippets(i + (k-1)*9,:) = ...
                    songs(i,(j-1)*178657 + 1:j*178657 - 1);
        end
    end
end
%% Creating spectogram shit
t = (1:178656)/Fs;
n=178656;
%each snipped is 4 seconds
L = 4;
k=(2*pi/L)*[0:n/2-1 -n/2:-1];
ks=fftshift(k);
%larg-ish t-slide cause else it matrix becomes to large to process
tslide = 0:0.5:L;

%% Getting Gaussian Spectogram matrices for train
Sgt_spec_train_matrix = zeros(270,length(tslide)*n);

width = 25;
for i = 1:270
    Sgt_spec_train = [];
    for j=1:length(tslide)
        g = exp(-width*(t - tslide(j)).^2); %gaussian filter
        Sg = g.*train_snippets(i,:); %filtered with gaussian
        Sgt = fft(Sg); %fft gaussian
        Sgt_spec_train = [Sgt_spec_train;abs(fftshift(Sgt))];
    end
    Sgt_spec_train_matrix(i,:) = ...
        reshape(Sgt_spec_train, 1, length(tslide)*n);
end
%% Getting Gaussian Spectogram matrices for test
Sgt_spec_test_matrix = zeros(27,length(tslide)*n);
for i = 1:27
    Sgt_spec_test = [];
    for j=1:length(tslide)
        g = exp(-width*(t - tslide(j)).^2); %gaussian filter
        Sg = g.*test_snippets(i,:); %filtered with gaussian
        Sgt = fft(Sg); %fft gaussian
        Sgt_spec_test = [Sgt_spec_test;abs(fftshift(Sgt))];
    end
    Sgt_spec_test_matrix(i,:) = ...
        reshape(Sgt_spec_test, 1, length(tslide)*n);
end
%% SVD
[U, S, V] = svd(Sgt_spec_train_matrix', 'econ');
%% Creating the projection matrix
projection_matrix = U.'*Sgt_spec_train_matrix.';
%% Creating test_projection vector
test_projection = U.'*Sgt_spec_test_matrix.';
%% Finding 3 avg projections
avg_projections = zeros(270,3);
for i = 1:3
    for j = 1:90
        avg_projections(:,i) = avg_projections(:,i) + ....
            projection_matrix(:,(i-1)*90 + j);
    end
    avg_projections(:,i) = avg_projections(:,i)/90;
end
%% Classifying test data
%for each row the expected out is [1,1,1,2,2,2,3,3,3]
indices = zeros(3,9);

for k = 1:3

    for i = 1:9

        min_index = 0;
        min_dist = Inf;

        for j = 1:3
            dist = sum(abs(avg_projections(:,j) - ...
                test_projection(:,i + (k-1)*9)));
            if dist < min_dist
                min_dist = dist;
                min_index = j;
            end
        end
        indices(k,i) = min_index;
    end

end
disp(indices);

%% saving spectogram train and test matrices as .csv
csvwrite("Spectogram_train_1.csv",Sgt_spec_train_matrix);
csvwrite("Spectogram_test_1.csv",Sgt_spec_test_matrix);