clear all; close all; clc;
%% constructing avg face matrix

directory_Path = 'Data/UnCropped/yalefaces/';
directory = dir(directory_Path);

faces = [];
averageFaces = [];

for i = 1:15
    
    sum_faces_type_1 = zeros(1,77760);
    number_of_faces_type_1 = 0;
    
    for j = 1:11
        %find each file for each face
        currentFilePath = strcat(directory_Path, '/',...
            directory((i-1)*11 + j + 2).name);
        
        %convert face from uint8 to double
        face = double(imread(currentFilePath));
        %disp(face);
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
[U, S, V] = svd(faces'/sqrt(77760-1), 'econ');
%% Drawing Eigen Faces
figure(1)
for i = 1:9
    subplot(3,3,i)
    eigenface = reshape(U(:,i),size(face, 1), size(face, 2));
    pcolor(flipud(eigenface)), shading INTERP, colormap(gray)
    title(sprintf('Eigen Face %d/165', i));
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
    projection = U(:,1:60).'*averageFaces(i, :).';
    subplot(3, 3, i)
    bar(projection), set(gca, 'Xlim', [0 100], 'Ylim', [-4000 4000]);
    title(sprintf('Projection of Avg Face %d on U_{short}', i));

end
%% Drawing Avg faces
figure(4)
for i = 1:9
    subplot(3, 3, i)
    imshow(uint8(reshape(averageFaces(i, :), ...
        [size(face, 1), size(face, 2)])));
    title(sprintf('Face %d', i));
end
%% Reconstructing faces
V_T = V.';
V_T = V_T(1:60,:);
reconstructed = (U(:,1:60)*S(1:60,1:60)*V_T).'*sqrt(77760-1);
figure(5)
%plotting 9 reconstructed average faces
for i = 1:9
    subplot(3, 3, i)
    avg_face_reconstructed = zeros(1,77760);
    for j = 1:11
        avg_face_reconstructed = avg_face_reconstructed + ...
            reconstructed(j + (i-1)*11,:);
    end
    pcolor(flipud(reshape(avg_face_reconstructed, ...
        [size(face, 1), size(face, 2)]))), ...
        shading INTERP, colormap(gray);
    title(sprintf('Reconstructed Face %d', i));
end