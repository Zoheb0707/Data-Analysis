clear all; close all; clc;
%% constructing avg face matrix

directory_Path = 'Data/Cropped/CroppedYale/';
directory = dir(directory_Path);

faces = [];
averageFaces = [];

for i = 3:length(directory)
    %finds each folder in CroppedYale
    current_directory = directory(i).name();
    %list of all file names in current folder
    file_list = dir(strcat(directory_Path, current_directory));
    %store sum of all faces of 1 type
    sum_faces_type_1 = zeros(1,32256);
    number_of_faces_type_1 = 0;
    
    for j = 3 : length(file_list)
        current_folder_path = strcat(directory_Path, current_directory);
        %find each file for each face
        currentFilePath = strcat(current_folder_path,...
            '/',file_list(j).name);
        
        %convert face from uint8 to double
        face = double(imread(currentFilePath));
        %reshape face to make it linear
        face_to_add = reshape(face, ...
            [1, size(face, 1)*size(face, 2)]);
        %add faces to matrix
        faces = [faces; face_to_add];
        
        %add face to sub matrix
        sum_faces_type_1 = sum_faces_type_1+face_to_add;
        number_of_faces_type_1 = number_of_faces_type_1 + 1;
    end
    %add average face to avg face matrix
    averageFaces = [averageFaces; ...
        sum_faces_type_1/number_of_faces_type_1];
end

%% De mean data
mn = mean(faces, 2);
for i=1:size(faces,1)
    faces(i,:) = faces(i,:) - mn(i,1);
end
%% SVD
[U, S, V] = svd(faces'/sqrt(32256-1), 'econ');
%% Drawing Eigen Faces
figure(1)
for i = 1:9
    subplot(3,3,i)
    eigenface = reshape(U(:,i),size(face, 1), size(face, 2));
    pcolor(flipud(eigenface)), shading INTERP, colormap(gray)
    title(sprintf('Eigen Face %d/2432', i));
end
figure(2)
diags = diag(S);
plot(log(diags(1:50)), 'ko', 'Linewidth', [2])
xlabel('Singular Value number')
ylabel('Log(Magnitude)')
title('Trend in Log(singular values) (1st 50)')
%% Projecting Avg faces
%plotting projections on bar graph
figure(3)
for i = 1:9
    projection = U(:,1:500).'*averageFaces(i, :).';
    subplot(3, 3, i)
    bar(projection), set(gca, 'Xlim', [0 100], 'Ylim', [-4000 4000]);
    title(sprintf('Projection of Avg Face %d on U_{short}', i));

end
%% Drawing Avg faces
figure(4)
for i = 1:9
    subplot(3, 3, i)
    imshow(uint8(reshape(averageFaces(i, :),...
        [size(face, 1), size(face, 2)])));
    title(sprintf('Face %d', i));
end
%% Reconstructing faces
V_T = V.';
V_T = V_T(1:500,:);
reconstructed = (U(:,1:500)*S(1:500,1:500)*V_T).'*sqrt(32256-1);
figure(5)
%plotting 9 reconstructed average faces
for i = 1:9
    subplot(3, 3, i)
    avg_face_reconstructed = zeros(1,32256);
    for j = 1:64
        avg_face_reconstructed = avg_face_reconstructed + ...
            reconstructed(j + (i-1)*64,:);
    end
    pcolor(flipud(reshape(avg_face_reconstructed, ...
        [size(face, 1), size(face, 2)]))), ...
        shading INTERP, colormap(gray);
    title(sprintf('Reconstructed Face %d', i));
end