function ND = nearestdrainage(FD,S,DEM)

%NEARESTDRAINAGE linear index to the nearest stream
%
% Syntax
%
%     ND = nearestdrainage(FD,S,DEM)
%
% Description
%
%     nearestdrainage finds the linear index to the nearest stream for each
%     point in the basin. This is just an altered verstion of
%     vertdistance2stream, which finds the HAND elevation of the nearest
%     drainage point. Here, I just took the index of that drainage point.
%
% Input arguments
%
%     FD    instance of FLOWobj
%     S     stream network (class: STREAMobj)
%     DEM   digital elevation model (class: GRIDobj)
%
% Output arguments
%
%     ND    nearest drainage point as a matrix

narginchk(3,3)

validatealignment(S,DEM);
validatealignment(FD,DEM);

ND = DEM;
ND.Z = NaN(DEM.size);
ND.Z(S.IXgrid) = S.IXgrid;

ix = FD.ix;
ixc = FD.ixc;
for r = numel(ix):-1:1
    if isnan(ND.Z(ix(r)))
        ND.Z(ix(r)) = ND.Z(ixc(r));
    end
end

ND = ND.Z;
