function [P, J] = anomaly_detection(filepath,  autoSeed, thresVal, seedThresh, minArea,maxAreaRatio, imgscale, maxDist, tfMean, tfFillHoles, tfSimplify)

% Syntax:  
%   P = anomaly_detection(filepath);
%   P = anomaly_detection(filepath, ridgetag, autoSeed)
%   P = anomaly_detection(..., thresVal, seedThresh, minArea, maxAreaRatio, maxDist, tfMean, tfFillHoles, tfSimpl)
%   [P, J] = anomaly_detection(...);
%
% Inputs:
%     filepath: Path of the intensity image
%     autoSeed: Boolean to indicate whether automate the seeding process {true}
%     thresVal: Absolute threshold level to be included     {5% of max-min}
%   seedThresh: [low, high] low and high bounds used to select seeds for 
%               low density and high density areas
%      minArea: Absolute threshold to rule out small regions
% maxAreaRatio: Relative threshold to rule out background/irrelevant
% regions
%     imgscale: downsampling for quick computation               {1.0}
%      maxDist: Maximum distance to the initial position in [px]      {Inf}
%       tfMean: Updates the initial value to the region mean (slow) {false}
%  tfFillHoles: Fills enclosed holes in the binary mask              {true}
%   tfSimplify: Reduces the number of vertices {true, if dpsimplify exists}
%
% Outputs:
%   P: VxN array (with V number of vertices, N number of dimensions)
%      P is the enclosing polygon for all associated pixel/voxel
%   J: Binary mask (with the same size as the input image) indicating
%      1 (true) for associated pixel/voxel and 0 (false) for outside
%   
% Requirements:
%   TheMathWorks Image Processing Toolbox for bwboundaries() and axes2pix()
%   Optional: Line Simplification by Wolfgang Schwanghart to reduce the 
%             number of polygon vertices (see the MATLAB FileExchange)
%
% Remarks:
%   The queue is not preallocated and the region mean computation is slow.
%   I haven't implemented a preallocation nor a queue counter yet for the
%   sake of clarity, however this would be of course more efficient.
% 
% Adapted from Region growing algorithm for intensity images
% contact: kai.wang.magic@gmail.com


% input arguments checking and default values
if nargin > 11
    error('Wrong number of input arguments!')
end

if ~exist('filepath', 'var')
    error('Please specify image file path!')
end

if ~exist('autoSeed', 'var')
    autoSeed= true;
end

if ~exist('seedThresh', 'var')
    seedThresh= [20, 200];
end

if ~exist('minArea', 'var')
    minArea= 10;
end

if ~exist('maxAreaRatio', 'var')
    maxAreaRatio= 0.15;
end

if ~exist('maxDist', 'var') || isempty(maxDist)
    maxDist = Inf;
end

if ~exist('tfMean', 'var') || isempty(tfMean)
    tfMean = false;
end

if ~exist('tfFillHoles', 'var')
    tfFillHoles = true;
end

if ~exist('imgscale', 'var')
    imgscale = 1.0;
end

I= imread(filepath);
I = imresize(I(:,:,1),imgscale); %downsize the image if you want a quick result


%% reinforce the edge and smooth the image
Itarget = imadjust(I);
Itarget = ordfilt2(Itarget,1,ones(3,3)); %3-by-3 minimum filter
Itarget = ordfilt2(Itarget,9,ones(3,3)); %3-by-3 maximum filter
Itarget = imclose(Itarget, ones(5));
cIM= Itarget;

nseed = 3;              % number of interative seed points
figure, imshow(I, [0 255]), hold all


%% initialize seed points
himage = findobj('Type', 'image');
if isempty(himage)
    himage = imshow(I, []);
end

if autoSeed == false
    % graphical user input for the initial position
    p = ginput(nseed);

    % get the pixel position concerning to the current axes coordinates
    for i=1:nseed
          initPos(i, 1) = round(axes2pix(size(cIM, 2), get(himage, 'XData'), p(i,2)));
          initPos(i, 2) = round(axes2pix(size(cIM, 1), get(himage, 'YData'), p(i,1)));
    end
else
    %automatic seeding
    seedMat= cIM> seedThresh(2);
    seedMat= seedMat | (cIM< seedThresh(1)& cIM >0);

    [initPos(:,1), initPos(:,2)]= find(cIM>seedThresh(2));
    [t1, t2]=find(cIM<seedThresh(1));
    initPos(end+1:end+size(t1,1),:)=[t1, t2];
    nseed= size(initPos,1);
end

if ~exist('thresVal', 'var') || isempty(thresVal)
    thresVal = double((max(cIM(:)) - min(cIM(:)))) * 0.05;
end


if isequal(ndims(cIM), 2)
    initPos(:,3) = 1;
elseif isequal(ndims(cIM),1) || ndims(cIM) > 3
    error('There are only 2D images and 3D image sets allowed!')
end

[nRow, nCol, nSli] = size(cIM);
imSize= nRow*nCol;
if any(initPos(:) < 1) || any(initPos(:,1) > nRow) || any(initPos(:,2) > nCol)
    error('Initial position out of bounds, please try again!')
end

if thresVal < 0 || maxDist < 0
    error('Threshold and maximum distance values must be positive!')
end

if ~isempty(which('dpsimplify.m'))
    if ~exist('tfSimplify', 'var')
        tfSimplify = true;
    end
    simplifyTolerance = 1;
else
    tfSimplify = false;
end

for iter= 1:nseed
    % check if this seed already was examined
    if seedMat(initPos(iter,1), initPos(iter,2))==0
        continue;
    end
    % initial pixel value
    regVal = double(cIM(initPos(iter,1), initPos(iter,2), initPos(iter,3)));

    % text output with initial parameters
    disp(['RegionGrowing Opening: Initial position (' num2str(initPos(iter,1))...
          '|' num2str(initPos(iter,2)) '|' num2str(initPos(iter,3)) ') with '...
          num2str(regVal) ' as initial pixel value!'])

    % preallocate array
    J = false(nRow, nCol, nSli);

    % add the initial pixel to the queue
    queue = [initPos(iter,1), initPos(iter,2), initPos(iter,3)];


    %%% START OF REGION GROWING ALGORITHM
    while size(queue, 1)
      % the first queue position determines the new values
      xv = queue(1,1);
      yv = queue(1,2);
      zv = queue(1,3);

      % .. and delete the first queue position
      queue(1,:) = [];

      % check the neighbors for the current position
      for i = -1:1
        for j = -1:1
          for k = -1:1

            if xv+i > 0  &&  xv+i <= nRow &&...          % within the x-bounds?
               yv+j > 0  &&  yv+j <= nCol &&...          % within the y-bounds?          
               zv+k > 0  &&  zv+k <= nSli &&...          % within the z-bounds?
               any([i, j, k])       &&...      % i/j/k of (0/0/0) is redundant!
               ~J(xv+i, yv+j, zv+k) &&...          % pixelposition already set?
               sqrt( (xv+i-initPos(iter,1))^2 +...
                     (yv+j-initPos(iter,2))^2 +...
                     (zv+k-initPos(iter,3))^2 ) < maxDist &&...   % within distance?
               cIM(xv+i, yv+j, zv+k) <= (regVal + thresVal) &&...% within range
               cIM(xv+i, yv+j, zv+k) >= (regVal - thresVal) % of the threshold?

               % current pixel is true, if all properties are fullfilled
               J(xv+i, yv+j, zv+k) = true; 
               seedMat(xv+i, yv+j, zv+k)= 0;
               % add the current pixel to the computation queue (recursive)
               queue(end+1,:) = [xv+i, yv+j, zv+k];

               if tfMean
                   regVal = mean(mean(cIM(J > 0))); % --> slow!
               end

            end        
          end
        end  
      end
    end
    %%% END OF REGION GROWING ALGORITHM


    % loop through each slice, fill holes and extract the polygon vertices
    P = [];
    for cSli = 1:nSli
        if ~any(J(:,:,cSli))
            continue
        end

        % use bwboundaries() to extract the enclosing polygon
        if tfFillHoles
            % fill the holes inside the mask
            J(:,:,cSli) = imfill(J(:,:,cSli), 'holes');    
            B = bwboundaries(J(:,:,cSli), 8, 'noholes');
        else
            B = bwboundaries(J(:,:,cSli));
        end

        newVertices = [B{1}(:,2), B{1}(:,1)];

        % simplify the polygon via Line Simplification
        if tfSimplify
            newVertices = dpsimplify(newVertices, simplifyTolerance);        
        end

        % number of new vertices to be added
        nNew = size(newVertices, 1);

        % append the new vertices to the existing polygon matrix
        if isequal(nSli, 1) % 2D
            P(end+1:end+nNew, :) = newVertices;
        else                % 3D
            P(end+1:end+nNew, :) = [newVertices, repmat(cSli, nNew, 1)];
        end
    end

    % text output with final number of vertices
    disp(['RegionGrowing Ending: Found ' num2str(length(find(J)))...
          ' pixels within the threshold range (' num2str(size(P, 1))...
          ' polygon vertices)!']);
 
    % rule out small regions and background or normal parts
    if length(find(J))< imSize*maxAreaRatio && length(find(J))> minArea
        plot(P(:,1), P(:,2), 'LineWidth', 2)
    end
end

