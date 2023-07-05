clc, clear all, close all;
viewFlag = 0;

% Parameters
structSize = 5;
holeFillSize = 26;
connSize = 26;
immaxHeight = 5;
smallBitSize = 100; % size of small bits to remove
resizeRatio = 0.5;
imageSlice = 1; % slice to view
minParticleVol = 200; % minimum volume of particles to keep
lengthPerPixel = (5.55e-6)/resizeRatio; % length per pixel in image
dilationPixels = 2;

% Image Directory
inFile = 'sample_34.tif'; % image path
outFile = 'sample_34_segmented.tif';


%% Load tif image
Itmp = imread(inFile); % read in first frame
info = imfinfo(inFile);
numStack = length(info);

% sample_34.tif has 1220 frames
% size(X, DIMENSION), 1=rows, 2=cols

% construct an empty 3D matrix to represent image
img = zeros(size(Itmp,1), size(Itmp,2), numStack);

% read each frame of image and store in 3D matrix
for ii = 1:numStack
    % imread(file name, frame index, ?, ?);
    currentImage = imread(inFile,ii,'Info',info);
    img(:,:,ii) = currentImage;
end

% scale image (lossy compression)
imgScaled = imresize3(img, resizeRatio);
disp(['Finished loading image ', inFile]);


%% Threshold segmentation
% perform Otsu 1 level threshold: multithresh(A,N), max N=20
thresh = multithresh(imgScaled,1)+400;
disp(['Otsu threshold is ', num2str(thresh)]);

% create another matrix originally of all zeros
% =1 if solid, =0 if void
Ibinarized = imgScaled > thresh;

% clean the image and remove extraneous bits
Iopen = bwareaopen(Ibinarized, smallBitSize, connSize);
disp('Removed small extraneous bits')

% fill holes
Ifill = imfill(Iopen, connSize, 'holes');
disp('Filled holes')

% get rid of small extraneous bits
% close with a spherical structuring element of size structSize
tmpImage = ones(size(Ifill))-Ifill;

% close this image to clean it
Iclose = -(imclose(tmpImage,strel('sphere', structSize)));
Iclose = ones(size(Iclose))+Iclose;

% dilate by about 2 pixels to restore pixels lost due to the slow transition in intensity from void to solid 
Idilate = imdilate(Iclose, strel('sphere', dilationPixels));

if viewFlag == 1
    printFrame(Idilate, numStack*resizeRatio, "Closed and dilated image");
end


%% Get distance map from grain surface
distmap = bwdist(~Iclose);

% invert distance map
distmap = -distmap;
distmap(~Iclose) = Inf;

% set a maximum distance
distmap2 = imhmax(-distmap, immaxHeight, connSize);
distmap2(~Iclose) = Inf;

disp('Computed distance map');


%% Perform a watershed segmentation
% technique for separating overlapping objects
wshed = watershed(distmap2, connSize);

disp('Watershed complete');


%% Multiply watershed with the cleaned image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cutImage = uint16(wshed).*uint16(Iclose);


%% Get statistics of connected components (i.e., of particles)
disp('Computing region properties...')
stats = regionprops(cutImage,'Area','PixelList');
vols = [stats.Area];

% cutImage dimensions
s1 = size(cutImage,1);
s2 = size(cutImage,2);
s3 = size(cutImage,3);

% create array to keep particles only over a certain size
particlesLarger = zeros(s1,s2,s3);
ind2keep = find(vols>minParticleVol);
stats = stats(ind2keep);
grainStats = zeros(numel(ind2keep),2);

for ii=1:numel(ind2keep)
    % indices of voxels for this grain
    a=stats(ii).PixelList(:,1);
    b=stats(ii).PixelList(:,2);
    c=stats(ii).PixelList(:,3);
    
    % assign to arrays
    voltmp = numel(a);
    grainStats(ii,1) = mean(a)/resizeRatio;
    grainStats(ii,2) = mean(b)/resizeRatio;
    grainStats(ii,3) = mean(c)/resizeRatio;
    grainStats(ii,4) = (voltmp/resizeRatio^3*3/4/pi)^(1/3);
    grainStats(ii,5) = voltmp/resizeRatio^3;
    
    % save to new array that we can output and view in ImageJ/Fiji
    for jj=1:numel(a)
        particlesLarger(b(jj),a(jj),c(jj))=ii;
    end
end

%% Save the data and stack

save('GrainStats_2.mat','grainStats');
fid=fopen('GrainStats_2.txt','w+');
for jj=1:size(grainStats,1)
    fprintf(fid,'%f %f %f %f %i\n',grainStats(jj,3)/5,grainStats(jj,2)/5,...
        grainStats(jj,1)/5,grainStats(jj,4)/5,grainStats(jj,5));
end
% the factor of 1/5 above accounts for scaling of our radiography images,
fclose(fid);

% save the stack
% particlesLarger = uint16(particlesLarger);
imwrite(particlesLarger(:,:,1),outFile);

for ii=2:size(particlesLarger,3)
    imwrite(particlesLarger(:,:,ii),outFile,'WriteMode','append');
end
disp('Finished writing stack');

% Show the stack
data_kept_rgb = label2rgb(particlesLarger(:,:,imageSlice),'jet','k','shuffle');
disp([num2str(numel(ind2keep)),' grains.']);
