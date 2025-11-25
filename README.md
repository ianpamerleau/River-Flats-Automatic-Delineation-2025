# River-Flats-Automatic-Delineation-2025
Code and Data for manuscript "Delineation of river flats with local hypsometry"

‘delineate_flats.m’ is the main script that outputs the river flats map. In addition to the two functions described below and provided, one will need the suite of funcitons from [Topotoolbox](https://www.mathworks.com/matlabcentral/fileexchange/50124-topotoolbox/). The only things needed as input for the code is a DEM file, the average latitude of the DEM, the resoluiton of the DEM, and the drainage area threshold.  
‘normDiff_thresh.m’ is the function that finds the ‘kernel flats elevation’ using the hypsometric analysis.  
‘nearestdrainage.m’ finds the index of the river pixel each watershed pixel drains into (all pixels not defined as part of the river are watershed pixels).  
