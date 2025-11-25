%% River Flats Delineation using Minima of Hyposmetric Curvature

clear
close all
clc
%% Set parameter values
% Load in the DEM as a GRIDobj. Variable 'name' was used to help directory
% sorting for different DEMs and printing/saving outputs with specific 
% names. This variable may not be useful for other users depending on how
% one labels and stores their DEMs. The varibale is not necessary for the
% calculations. 
name = 'Sagavanirktok';

% the ac variable is the drainage area threshold. It will define what is 
% and is not channel. May need to be adjusted by user. 
ac = 1000; % Drainage area threshold in km^2 

DEM = GRIDobj([name, '_dem.tif']);
% DEM = GRIDobj([name, '/', name, '_dem.tif']);

% imagesc(DEM)
% look at the DEM and take the average latitude. Enter it in the variable
% "lat" and proceed.
lat = 69.5; % ENTER LATITUDE HERE

%% Load and Process Maps
% This algorithm needs to put the DEM in resolution of km, so if DEM is 
% in arcsec, find the latitude of DEM for the km conversion. If DEM is 
% already in km/m do not do the arcsec conversion (but do change from m to 
% km). 

% If DEM is in arcsec, use this: ##########################################
res = cos(pi/180*lat)*(1852/60)*3; % Pixel resolution in meters at given latitude
DEM.cellsize = res/1000; % units now km
FD = FLOWobj(DEM,'preprocess','carve');
A = flowacc(FD).*(FD.cellsize)*(30.86666667*3/1000); %Flow accumulation grid (3 arcsec)
% A = flowacc(FD).*(FD.cellsize)*(30.86666667/1000); %Flow accumulation grid (1 arcsec)
% #########################################################################

% If DEM is already in km or m, use this: #################################
% DEM.cellsize = DEM.cellsize./1000; % meters to km
% FD = FLOWobj(DEM,'preprocess','carve'); 
% A = flowacc(FD); % Flow Accumulation grid
% A = flowacc(FD)*FD.cellsize; % Flow Accumulation grid
% #########################################################################

S = STREAMobj(FD, A>=ac); % define stream based on drainage area threshold
A = A.Z; % flow accumulation grid as a matrix (not GRIDobj) 
max_ix = find(A == max(A(:))); % linear index vector of the outlet of the river
L = drainagebasins(FD, max_ix); % define drainage basin using outlet
L = L>0; %Binary Drainage Basin: 1 = in the drainage basin, 0 = not 
L = L.Z; % Matrix form of the Binary Drainage Basin (not GRIDobj)
zdiff = vertdistance2stream(FD,S,DEM); % Height above nearest drainage
zdiff=zdiff.Z; % Height above nearest drainage matrix form

%% Find maximum radius to pad matrices 
% padarray to avoid out of bounds errors. Find the largest possible radius
% from the outlet river pixel using hydraulic geometry. This will be
% largest due to the power law relation between drainage area and valley
% width. 

c = 0.3; % power law taken from literature
k = 0.84; % adjust calculation to pixels
W_v = k*A(max_ix)^c; % diameter in terms of km
W_v = W_v*1000; % diameter in terms of m
W_v = W_v/res;  % diameter in terms of pixels
W_v = W_v/2;    % radius in terms of pixels
rad_max = double(int16(W_v)); 

%% Preprocess matrices for kernel algorithm 
% padarray to avoid out of bounds errors.
L = padarray(L, [rad_max, rad_max],0, 'both');
A = padarray(A, [rad_max, rad_max],nan, 'both');
zdiff = padarray(zdiff, [rad_max, rad_max],nan, 'both');

A(~L) = 0; % take only river pixels within drainage basin
A2 = A>=ac; %Binary river pixel grid: 1=river pixel, 0=not

A2(~L) = 0; % Set non-basin points to 0

% make a HAND grid that does not include channel points, we do not want 
% them in the kenrel algorithm. the hypsometric analysis will only use
% integers and no river pixels. Set non-basin point to NaN (0 would imply
% part of the river and skew the hyposmetric analysis.
zdiff(~L) = NaN; 
zdiff_int = zdiff;
zdiff_int(A2) = NaN;
zdiff_int = int16(zdiff_int); 

sz = size(A2); 
inKernel = zeros(sz, 'int16');  % # times seen within a kernel
plains = zeros(sz, 'int16');    % # times appeared in kernel and identified as river flat
subKernel = zeros(sz, 'logical');% kernels that do not have enough elevation diversity

%% Kernel algorithm

idx = find(A2); % vector index of all indices where there's a river pixel

% Iterate over all channel points
for i = idx'
    
    % Create a structuring element: F defines the kernel
    % find the radius of this river pixel based on drainage area
    W_v = k*A(i)^c; % diameter in terms of km
    W_v = W_v*1000; % diameter in terms of m
    W_v = W_v/res;  % diameter in terms of pixels
    W_v = W_v/2;    % radius in terms of pixels
    rad = double(int16(W_v)); 
    sz2 = [rad*2+1, rad*2+1];
    SE = strel('disk', rad, 0);
    F = int16(SE.Neighborhood); %kernel
    
    % Define the kernel centered at the current river pixel:
    [ptx, pty] = ind2sub(sz, i);
    
    F_z = zdiff_int(ptx-rad:ptx+rad, pty-rad:pty+rad); %the kernal circle
    F_z(F_z<0) = 0;
    F_z(~F) = nan;
    
    zvals=F_z(F & ~isnan(F_z) & F_z~=-9999); %Elevation values within kernel
    
    % add condition: if (# of unique values or values larger than 0) <= 5
    % mark channel pt and do not include it in the algorithm 
    uvals = unique(zvals); 
    
    if length(uvals) <= 5 % 5 was chosen because of how the hypsometric analysis works
        % the number 5 may need to increase if half_diff in normDiff_thresh
        % increases (due to the numerical differentiation method) 
        subKernel(ptx, pty) = 1; % not enough unique elevation, hypsometeric analysis fails
    else        
        % Find zthresh, the maximum height where a river flat occurs in this area
        zthresh = normDiff_thresh(zvals, 1);
        
        plains1 = zeros(sz2, 'int16');
        plains1(F_z<=zthresh)=true;
        plains1(~F) = false;

        %Update the global values with this kernel
        plains(ptx-rad:ptx+rad, pty-rad:pty+rad) = plains(ptx-rad:ptx+rad, pty-rad:pty+rad)+plains1;
        inKernel(ptx-rad:ptx+rad, pty-rad:pty+rad) = inKernel(ptx-rad:ptx+rad, pty-rad:pty+rad)+F;
    end
end     

%% Compute the likelihood of a pixel being flooded by kernel
likelihoodGrid = single(plains)./single(inKernel); % # times identified/#total times seen
likelihoodGrid(isnan(zdiff_int))=nan; % nan means pixel is not in the basin

%% get elevations with 40-60% likelihood of flooding for each river pixel

elevations = NaN(size(A));
% with current center and tolerance value, likelihood range of 40-60% will
% be found. These values can be adjusted to find different geologic river
% flats. 
tol = 0.1; % set tolerance for percentiles
center = 0.5; % center of percentile range tolerance will work around

idx = find(A2 & ~subKernel);

for i = idx'
    
    % Create a structuring element: F defines the kernel
    % recompute the radius for each river pixel to assess the 40-60%
    % likelihood elevation 
    W_v = k*A(i)^c; % diameter in terms of km
    W_v = W_v*1000; % diameter in terms of m
    W_v = W_v/res;  % diameter in terms of pixels
    W_v = W_v/2;    % radius in terms of pixels
    rad = double(int16(W_v));
    SE = strel('disk', rad, 0);
    F = int16(SE.Neighborhood); %kernel
    
    % define the kernel centered at point pt
    [ptx, pty] = ind2sub(sz, i);
    
    % likelihoodGrid within kernel
    F_likelihoodGrid = likelihoodGrid(ptx-rad:ptx+rad, pty-rad:pty+rad); %the kernal circle
    F_likelihoodGrid(~F) = NaN;

    % HAND elevation within kernel
    F_z = zdiff(ptx-rad:ptx+rad, pty-rad:pty+rad); %the kernal circle
    F_z(F_z<0) = 0;
    F_z(~F) = nan;
    
    test = find(F_likelihoodGrid>=(center-tol) & F_likelihoodGrid<=(center+tol) & F);
    
    % if there is no pixel that falls in the target range, make an exception
    if isempty(test)
        % find all pixels greater than the upper bound of the range chosen
        greaterThanUpper = find(F_likelihoodGrid>center+tol & F);
        
        % if nothing is greater than the upper bound, find highest 
        % likelihood below the lower bound
        if isempty(greaterThanUpper)
            elevations(i) = max(F_z(F_likelihoodGrid == max(F_likelihoodGrid(:))));
            
        % otherwise use lowest max elevation of those above the upper bound
        else
            elevations(i) = max(F_z(F & F_likelihoodGrid == min(F_likelihoodGrid(greaterThanUpper))));
        end
        
    % this is the normal condition we expect    
    else
        elevations(i) = median(F_z(F & abs(F_likelihoodGrid - center)<=tol));
    end
end

% note we find the median value in the requested range because we want an
% average of that range. If there is nothing in that range, we find the
% lowest liklihood value ABOVE that range (more likely flooded). If there 
% is no likelihood value in that range or above it, we find the highest
% liklihood value BELOW the range. We find the single value above/below
% instead of a median because it is closer to representing the requested 
% range value than a median values outside the range (e.g., 61% vs median 
% of 60-100%).

%% un-pad grids that were altered for kernel algoirthm 
% when ND (nearest drainage) was created, it used the original DEM.size to keep track of
% linear indexes. We padded the maps to avoid an out of bounds error during
% the kernel agorithm. Now we need to un-pad them to get the original
% size

% un-pad matrices (or grids) 
elevations = elevations(rad_max+1:end-rad_max, rad_max+1:end-rad_max);
A = A(rad_max+1:end-rad_max, rad_max+1:end-rad_max);
A2 = A2(rad_max+1:end-rad_max, rad_max+1:end-rad_max);
zdiff = zdiff(rad_max+1:end-rad_max, rad_max+1:end-rad_max);
likelihoodGrid = likelihoodGrid(rad_max+1:end-rad_max, rad_max+1:end-rad_max);
L = L(rad_max+1:end-rad_max, rad_max+1:end-rad_max);

% Create nearest drainage grid (output is matrix form, not GRIDobj). For
% each pixel in the basin, finds the river pixel (defined by our drainage
% area threshold) which the basin pixel drains into. 
ND = nearestdrainage(FD, S, DEM);
ND(~L) = NaN;

%% Flood Basin Based on Nearest Drainage

DB_idx = find(L); % linear index of each pixel in the drainage basin
rf_marked = zeros(size(A)); % River Flats map. Value will be 0 if it is not 
% part of the river flat. If a pixel is in the river flat, the value will
% be the likely flooding height. 

for i = DB_idx'
    
    cp_idx = ND(i); % find the channel point the pixel drains to
    cp_Z = elevations(cp_idx); % use the elevations grid to find flooding height of river pixel
    
    if cp_Z >= zdiff(i) % if river flooding height > elevation of basin pixel, 
        rf_marked(i) = cp_Z; % mark basin pixel with flooding height 
    end
    
end

% make binary grid (logical matrix, 1 = river flat, 0 = not river flat
rf_binary = rf_marked>0;
% set non-river flat pixels to NaN value (not a number) 
rf_marked(rf_marked<=0) = NaN;

%% Show map
% Final River flats product 
close all
g=figure;
imagesc(rf_marked, 'AlphaData', rf_marked>0)
cb=colorbar;
cb.Label.String = '[m]';
title(name)

%% Save file as Geotiff

% File and folder names may need to be changed for user.

% [~,R_geo] = geotiffread([name, '/', name, '_dem.tif']);
% info = geotiffinfo([name, '/', name,'_dem.tif']);
% geotiffwrite([name, '/', name,'_river_pxls.tif'], A2, R_geo, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
% geotiffwrite([name, '/', name,'_fp_floodingElevations.tif'], rf_marked, R_geo, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
% geotiffwrite([name, '/', name,'_fp_binary_200rad.tif'], rf_binary, R_geo, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
