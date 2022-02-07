function [binmatrix] = ccp_setup_grid(minlat,maxlat,minlon,maxlon,bspace,bsize)
%CCP_SETUP_GRID Setup geographical binning matrix
%   CCP_SETUP_GRID (minlat,maxlat,minlon,maxlon,bspace,bsize) sets up the
%   geographical binning space.  Specify the lat/lon range of the
%   geographical area, and set the spacing between each nodes for the x
%   and y directions.  The spacing is given in kilometers.  The size of
%   each bin is given by a radius for which the node is the center.
%
%   INPUT:
%          minlat = minimum geographic area latitude
%          maxlat = maximum geographic area latitude
%          minlon = minimum geographic area longitude
%          maxlon = maximum geographic area longitude
%          bspace = node spacing (in kilometers)
%          bsize = size of bins (in kilometers)
%
%   NOTES:
%   1. Since it is hard to estimate a range and spacing that gives you whole
%   number nodes, we take the range divided by the spacing and round up.
%   We then take that integer and multiply by the spacing.  This gives us a
%   new range, which we center in the original input geographic coordinates.
%   2. Because the spacing is given in kilometers, the number of nodes in
%   each row of the grid changes as you move in latitude due to the
%   conversion between kilometers and degrees around the small circles.
%   That way, the grid doesn't become finer as you move towards the poles,
%   which would not make a lot of sense in data analysis.
%   3. Using the circle option only outputs the bins and the distance
%   (radius) constraints of each bin.  Therefore, it requires a little bit
%   more programming setup at the next step in order to use these bins
%   (i.e. calculate the distance to each node and assign the point to those bins
%   for which the distance calculated is less than or equal to the radius).
%   4. Bin numbering starts from southwest and goes east.
%          11   12   13    14   15
%          6     7    8     9   10
%          1     2    3     4    5
%
%   EXAMPLE:
%   binmatrix = ccp_setup_grid (35,48,-126,-110,100,200);
%   Output binmatrix (304x5) columns:
%   Bin # | node lat | node lon | radius (km) | "blank"
%
%   node lat | node lon | radius (km) | "blank cell" | "blank cell"

%   Author: Kevin C. Eagar
%   Date Created: 02/20/2007
%   Last Updated: 05/18/2010
% Feb. 12, 2017, Yunfeng Chen, force ynodes and xnodes to be integer,
% otherwise may cause round off problem (i.e., missing nodes at center 
% latitude)
% calculate y nodes first
%--------------------------------------------------------------------------
deltay = km2deg(bspace);
oldrangey = abs(minlat-maxlat);
yintervals = ceil(oldrangey/deltay);
newrangey = yintervals*deltay;
diffyrange = (newrangey-oldrangey)/2;
maxlat = maxlat+diffyrange;
minlat = minlat-diffyrange;
ynodes = yintervals+1;

% set up binning matrix and calculate x nodes using km2deg conversion at
% each latitude as we go 
%--------------------------------------------------------------------------
bincount = 0;
for ny = 1:ynodes
    nodelat = minlat+((ny-1)*deltay);

    % calculate the x nodes
    %--------------------------------------------------------------
    kmperdeg = 2*pi*6371*sind(90-nodelat)/360;
    deltax = bspace/kmperdeg;
    oldrangex = abs(minlon-maxlon);
    xintervals = oldrangex/deltax;
    newrangex = ceil(xintervals)*deltax;
    diffxrange = (newrangex-oldrangex)/2;
    ymaxlon = maxlon+diffxrange;
    yminlon = minlon-diffxrange;
    xnodes = round(newrangex/deltax);
    for nx = 1:xnodes+1
        nodelon = yminlon+((nx-1)*deltax);
        bincount = bincount+1;
        binmatrix{bincount,1} = nodelat;
        binmatrix{bincount,2} = nodelon;
        binmatrix{bincount,3} = bsize;
        binmatrix{bincount,4} = {};
        binmatrix{bincount,5} = {};
    end
end