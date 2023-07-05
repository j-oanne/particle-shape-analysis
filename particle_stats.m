
%% Load the particle tif image

filename = '77.tif'; % change to file path

Itmp = imread(filename);
info = imfinfo(filename);
numStack = length(info);

% construct an empty 3D matrix to represent particle
img = zeros(size(Itmp,1), size(Itmp,2), numStack);

for ii=1:numStack
    % read in each frame
    img(:,:,ii) = imread(filename,ii,'Info',info);
end

disp(['File loaded: ', filename]);


%% Get particle form information
s = regionprops3(img,'all');

% principal dimensions
PA = sort(s.PrincipalAxisLength,'descend')';

a = PA(1); % long axis
b = PA(2); % intermediate axis
c = PA(3); % short axis

% aspect ratios
EI = b/a; % elongation index
FI = c/b; % flatness index
AR = (EI+FI)/2; % characteristic aspect ratio

disp(['Elongation Index: ', num2str(EI)]);
disp(['Flatness Index: ', num2str(FI)]);
disp(['Aspect Ratio: ', num2str(AR)]);


%% Particle roundness
% quantified with the local curvature distribution at a particular surface

N=1280; % number of triangular elements


% Marching Cubes with Gaussian Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
img(img>0)=1; % binarize matrix
img_org=img;

img=smooth3(img,'gaussian',5); % reduce noise in data set

isovalue=0.5;
[F,V] = extractIsosurface(img,isovalue);


% patch triangular mesh for particle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
view(3);
p = patch('Faces',F,'Vertices',V);
reducepatch(p, N);
% set(p, 'Visible', 'off');

% set figure presets
set(p,FaceColor=[0.5 0.5 0.5]);
set(p,FaceAlpha=0.4);
% set(p,EdgeColor="none");
camlight;
axis equal;
grid on;


% get maximum inscribed sphere
D = bwdist(~img_org, "euclidean");
[radius,index] = max(D, [], 'all');
[y,x,z] = ind2sub(size(D), index); % row, col, stack
center = [x,y,z];

% create sphere
[vv,ff]=icosphere(4);
vv = vv*radius+center.*ones(size(vv));

sph = patch('Faces',ff,'Vertices',vv);
reducepatch(sph, N);


% Gauss curvature of a sphere is K=1/r^2
% Mean curvature of a sphere is H=1/r
[Cmean,Cgaussian,~,~,~,~] = patchcurvature(p, true);
H = 1/radius;

% loop backwards to prevent concurrent modification errors
for ii=length(Cmean):-1:1
    % remove curvatures less than sphere because they are not corners
    if Cmean(ii) < H
        Cmean(ii) = [];
    end
end

R_M = mean(Cmean); % mean curvature roundness index
disp(['Mean curvature: ', num2str(R_M)]);


%% Compactness: sphericity and convexity

% volume and surface area should be returned as a scalar
% why are these returning matrices

vol = s.Volume(numel(s.Volume));
SA = s.SurfaceArea(numel(s.SurfaceArea));

S = ((36*pi*vol*vol)^(1/3))/SA;
Cx = vol/s.ConvexVolume(numel(s.ConvexVolume));

disp(['Sphericity: ', num2str(S)]);
disp(['Convexity: ', num2str(Cx)]);

