%{ 
Hypsometric Curvature Minimum Calculator 

To find flooding height within a kernel, we want to find where the
hypsometric curve transitions from high slope (corresponding to river flats
region) to low slope (corresponding to steeping, bounding terrain). This
change in slope should be seen as a minimum of the curvature of the
hypsometric curve. This function derives the hypsometric curvature
(technically it is the normalized difference of the slope, but we found
this lead to more prominant minima at river flats elevations), finds the
index of that minimum, and returns the elevation corresponding to that
index. 
Please note: you will lose [half_ddiff] * 2 nodes on each side of the vector
when the curvature is calcuated. We suggest using a [half_diff] value of 1 
as it loses less data (when using integers as recommended, 0 m and 1 m will
likely be lost and 2 m will be the lowest flat elevation). 

INPUT VARIABLES
zvals: Nx1 vector of all the elevations found in the kernel.
half_diff: the spread of local values we use to find the slope at each
elevation value, must be integer (recommend 1). 

OUTPUT
thresh: the flooding elevation based on the minimum of the curvature of the
elevations withint the kernel

%}

% zvals should be fed in as integer values in the main script. If not,
% make sure they are integers here. However, some functions listed
% below need single types as input.

function thresh = normDiff_thresh(zvals, half_diff)

    if isempty(zvals) 
        thresh=0;
        return;
    end

    zvals = single(zvals); % convert to single floats 
    [uvals, ~, ic] = uniquetol(zvals); % unique set of elevations (needs singles)
    
    if (length(uvals)- 4*half_diff) < 10
        % otherwise there is a small number of 2nd derivatives, 4* because 2 at each side of vector
        thresh=0;
        return
    end
    cdf = cumsum(accumarray(ic,1)); % cummulative distribution function 
    
    % get values of cdf and elevations spaced [smooth] away from each other
    smooth=2*half_diff+1;
    area_high = cdf(smooth:end);
    area_low = cdf(1:end-smooth+1);
    % must be converted to double in order to calculate da_dz below
    z_high = double(uvals(smooth:end)); 
    z_low = double(uvals(1:end-smooth+1));

    % compute slope and normalized difference of each point
    da_dz = (area_high-area_low)./(z_high-z_low); % area / elevation in meters

    % repeat process with hypsometric slope values (NOT topographic slope)
    da_dz_high = double(da_dz(smooth:end));
    da_dz_low = double(da_dz(1:end-smooth+1));
       
   
    % calculate normalized difference of the slope
    normDiff0 = (da_dz_high-da_dz_low)./(da_dz_high+da_dz_low);
    

    % fit the two vectors back to the same size as the elevations vector.
    % This is needed for indexing at the end of the function and is helpful
    % if you want to make plots of what is happening inside the function
   
    normDiff = [ones(2*half_diff,1)*max(normDiff0); normDiff0; ones(2*half_diff,1)*max(normDiff0)];
    
    % find the minimums of the normDiff
    min_list = find(islocalmin(normDiff));
    idx = [];
    
    % go through each minimum and check to see if it is a "prominant"
    % minimum (i.e., the two normDiff values on either side of the minimum
    % are greater than their other neighbor)
    for i = 1:length(min_list)
        if normDiff(min_list(i)-1)<=normDiff(min_list(i)-2) && normDiff(min_list(i)+1)<=normDiff(min_list(i)+2)
            idx = min_list(i); 
            break % when you find one, break out of for loop
        end
    end
    
    % if no minimums are "prominant" minimums, then use the last minimum
    if isempty(idx)
        idx = find(min_list(end));%ES xx xhange
    end
    
    % if there are no minimums, use the last elevation value
    if isempty(idx)
        thresh = uvals(length(uvals));
        return
    else   
        % this is the usual case: there will be a prominant minimum that
        % was found in the first for loop, and is being used now to find
        % the elevation value associated with it
        thresh = uvals(idx);
        return
    end
end 