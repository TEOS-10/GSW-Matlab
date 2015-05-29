function [delta_SA_ref, in_ocean] = gsw_delta_SA_ref(p,long,lat)

% gsw_delta_SA_ref                Absolute Salinity anomaly reference value
%                                                (excluding the Baltic Sea)
%==========================================================================
%
% USAGE:  
%  [delta_SA_ref, in_ocean] = gsw_delta_SA_ref(p,long,lat)
%
% DESCRIPTION:
%  Calculates the Absolute Salinity anomaly reference value, SA - SR, in 
%  the open ocean by spatially interpolating the global reference data set  
%  of deltaSA_ref to the location of the seawater sample.  
% 
%  The Absolute Salinity anomaly reference value in the Baltic Sea is 
%  evaluated separately, since it is a function of Practical Salinity, not  
%  of space.  The present function returns a delta_SA_ref of zero for data  
%  in the Baltic Sea.  The correct way of calculating Absolute Salinity in  
%  the Baltic Sea is by calling gsw_SA_from_SP.  
%
% INPUT:
%  p     =  sea pressure                                           [ dbar ] 
%          ( i.e. absolute pressure - 10.1325 dbar )
%  long  =  Longitude in decimal degrees                     [ 0 ... +360 ]
%                                                      or [ -180 ... +180 ]
%  lat   =  Latitude in decimal degrees north               [ -90 ... +90 ]
%
%  p, long & lat need to be vectors and have the same dimensions.
%
% OUTPUT:
%  delta_SA_ref  =  Absolute Salinity anomaly reference value      [ g/kg ]
%  in_ocean      =  0, if long and lat are a long way from the ocean 
%                =  1, if long and lat are in the ocean
%  Note. This flag is only set when the observation is well and truly on
%    dry land; often the warning flag is not set until one is several 
%    hundred kilometres inland from the coast. 
%
% AUTHOR: 
%  David Jackett                                      [ help_gsw@csiro.au ]
%
% MODIFIED:
%  Paul Barker and Trevor McDougall
%
% VERSION NUMBER: 3.0 (23rd March, 2011)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm 
%   for estimating Absolute Salinity in the global ocean.  Submitted to 
%   Ocean Science. A preliminary version is available at Ocean Sci. Discuss.,
%   6, 215-242.  
%   http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_delta_SA_ref:  Requires three inputs')
end %if

[mp,np] = size(p);
[mla,nla] = size(lat);
[mlo,nlo] = size(long);

if (mp ~= mla) | (mp ~=mlo) | (np ~= nla) | (np ~= nlo)
    error('gsw_delta_SA_ref: Inputs need be of the same size')
end %if

%--------------------------------------------------------------------------
% Start of the calculation (extracting from a look up table)
%--------------------------------------------------------------------------

gsw_data = 'gsw_data_v3_0.mat';

gsw_data_file = which(gsw_data);

load (gsw_data_file,'deltaSA_ref','lats_ref','longs_ref','p_ref',...
    'ndepth_ref');

nx = length(longs_ref); 
ny = length(lats_ref); 
nz = length(p_ref); 
n0 = length(p);

dlongs_ref = longs_ref(2) - longs_ref(1); 
dlats_ref = lats_ref(2) - lats_ref(1);

indsx0 = floor(1 + (nx-1)*(long - longs_ref(1))./(longs_ref(nx) - longs_ref(1)));
indsx0 = indsx0(:); 
inds = find(indsx0 == nx); 
indsx0(inds) = nx - 1;
              
indsy0 = floor(1 + (ny-1)*(lat - lats_ref(1))./(lats_ref(ny) - lats_ref(1)));
indsy0 = indsy0(:); 
inds = find(indsy0 == ny); 
indsy0(inds) = ny - 1;

indsz0 = sum(ones(nz,1)*p(:)' >= p_ref(:)*ones(1,n0));
indsz0 = indsz0(:);                             % adjust in the vertical                                            
                                            
indsn1 = sub2ind([ny,nx],indsy0,indsx0);        % casts containing data
indsn2 = sub2ind([ny,nx],indsy0,indsx0+1);
indsn3 = sub2ind([ny,nx],indsy0+1,indsx0+1);
indsn4 = sub2ind([ny,nx],indsy0+1,indsx0);

nmax = max([ndepth_ref(indsn1)';ndepth_ref(indsn2)';ndepth_ref(indsn3)';ndepth_ref(indsn4)']);

inds1 = find(indsz0(:)' > nmax);                % casts deeper than GK maximum

p(inds1) = p_ref(nmax(inds1));                  % have reset p here so have to reset indsz0

indsz0 = sum(ones(nz,1)*p(:)' >= p_ref(:)*ones(1,n0));
indsz0 = indsz0(:); 
inds = find(indsz0 == nz); 
indsz0(inds) = nz - 1;

inds0 = sub2ind([nz,ny,nx],indsz0,indsy0,indsx0);
   
data_indices = [indsx0,indsy0,indsz0,inds0]; 
data_inds = data_indices(:,3); 
    
r1 = (long(:) - longs_ref(indsx0))./(longs_ref(indsx0+1) - longs_ref(indsx0));
s1 = (lat(:) - lats_ref(indsy0))./(lats_ref(indsy0+1) - lats_ref(indsy0));
t1 = (p(:) - p_ref(indsz0))./(p_ref(indsz0+1) - p_ref(indsz0));
    
nksum = 0;
no_levels_missing = 0;

sa_upper = nan(size(data_inds)); 
sa_lower = nan(size(data_inds));
delta_SA_ref = nan(size(data_inds));
in_ocean = ones(size(delta_SA_ref));

for k = 1:nz-1
    
    inds_k = find(indsz0 == k);
    nk = length(inds_k);
    
    if nk>0
        nksum = nksum+nk;
        indsx = indsx0(inds_k);
        indsy = indsy0(inds_k);
        indsz = k*ones(size(indsx));
        inds_di = find(data_inds == k);             % level k interpolation
        dsa = nan(4,n0);
        inds1 = sub2ind([nz,ny,nx], indsz, indsy, indsx);
        dsa(1,inds_k) = deltaSA_ref(inds1);
        inds2 = sub2ind([nz,ny,nx], indsz, indsy, indsx+1);
        dsa(2,inds_k) = deltaSA_ref(inds2);                % inds0 + ny*nz
        inds3 = sub2ind([nz,ny,nx], indsz, indsy+1, indsx+1);
        dsa(3,inds_k) = deltaSA_ref(inds3);           % inds0 + ny*nz + nz
        inds4 = sub2ind([nz ny,nx], indsz, indsy+1, indsx);
        dsa(4,inds_k) = deltaSA_ref(inds4);                   % inds0 + nz
                       
        inds = find(260<=long(:) & long(:)<=295.217 & ...
            0<=lat(:) & lat(:)<=19.55 & indsz0(:)==k);
        if ~isempty(inds)
            dsa(:,inds) = gsw_dsa_add_barrier(dsa(:,inds),long(inds), ...
                lat(inds),longs_ref(indsx0(inds)),lats_ref(indsy0(inds)),dlongs_ref,dlats_ref);
        end
        
        inds = find(isnan(sum(dsa))' & indsz0==k);
        if ~isempty(inds)
            dsa(:,inds) = gsw_dsa_add_mean(dsa(:,inds));
        end
        
        sa_upper(inds_di) = (1-s1(inds_di)).*(dsa(1,inds_k)' + ...
            r1(inds_di).*(dsa(2,inds_k)'-dsa(1,inds_k)')) + ...
            s1(inds_di).*(dsa(4,inds_k)' + ...
            r1(inds_di).*(dsa(3,inds_k)'-dsa(4,inds_k)'));  % level k+1 interpolation
                
        dsa = nan(4,n0);
        inds1 = sub2ind([nz,ny,nx], indsz+1, indsy, indsx);
        dsa(1,inds_k) = deltaSA_ref(inds1);
        inds2 = sub2ind([nz,ny,nx], indsz+1, indsy, indsx+1);
        dsa(2,inds_k) = deltaSA_ref(inds2);                % inds1 + ny*nz
        inds3 = sub2ind([nz,ny,nx], indsz+1, indsy+1, indsx+1);
        dsa(3,inds_k) = deltaSA_ref(inds3);           % inds1 + ny*nz + nz
        inds4 = sub2ind([nz ny,nx], indsz+1, indsy+1, indsx);
        dsa(4,inds_k) = deltaSA_ref(inds4);                   % inds1 + nz
                
        inds = find(260<=long(:) & long(:)<=295.217 & ...
            0<=lat(:) & lat(:)<=19.55 & indsz0(:)==k);
        if ~isempty(inds)
            dsa(:,inds) = gsw_dsa_add_barrier(dsa(:,inds),long(inds), ...
                lat(inds),longs_ref(indsx0(inds)),lats_ref(indsy0(inds)),dlongs_ref,dlats_ref);
        end
        
        inds = find(isnan(sum(dsa))' & indsz0==k);
        
        if ~isempty(inds)
            dsa(:,inds) = gsw_dsa_add_mean(dsa(:,inds));
        end
        
        sa_lower(inds_di) = (1-s1(inds_di)).*(dsa(1,inds_k)' + ...
            r1(inds_di).*(dsa(2,inds_k)'-dsa(1,inds_k)')) + ...
            s1(inds_di).*(dsa(4,inds_k)' + ...
            r1(inds_di).*(dsa(3,inds_k)'-dsa(4,inds_k)'));
        
        inds_different = find(isfinite(sa_upper(inds_di)) & isnan(sa_lower(inds_di)));
        
        if ~isempty(inds_different)
            sa_lower(inds_di(inds_different)) = sa_upper(inds_di(inds_different));
        end
        
        delta_SA_ref(inds_di) = sa_upper(inds_di) + t1(inds_di).*(sa_lower(inds_di) - sa_upper(inds_di));
        
    else
        no_levels_missing = no_levels_missing + 1;
    end
end

inds = find(~isfinite(delta_SA_ref)); 
delta_SA_ref(inds) = 0;

in_ocean(inds) = 0;

end

%##########################################################################

function delta_SA_ref = gsw_dsa_add_mean(dsa)

% gsw_dsa_add_mean
%==========================================================================
%
% USAGE:
%  delta_SA_ref = gsw_dsa_add_mean(dsa)
%
% DESCRIPTION:
%  Replaces NaN's with nanmean of the 4 adjacent neighbours
%
% INPUT:
%  dsa  =  Absolute Salinity anomaly of the 4 adjacent neighbours  [ g/kg ]
%
% OUTPUT:
%  delta_SA_ref  =  nanmean of the 4 adjacent neighbours           [ g/kg ]
%
% AUTHOR: 
%  David Jackett
%
% VERSION NUMBER: 3.0
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T. J., D. R. Jackett and F. J. Millero, 2010: An algorithm 
%   for estimating Absolute Salinity in the global ocean.  Submitted to 
%   Ocean Science, a preliminary version is available at Ocean Sci. Discuss.,
%   6, 215-242.  
%   http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf
%   and the computer software is available from http://www.TEOS-10.org
%
%==========================================================================

dsa_mean = mean(dsa); 
inds_nan = find(isnan(dsa_mean)); 
no_nan = length(inds_nan);

for kk = 1:no_nan
    col = inds_nan(kk);
    inds_kk = find(isnan(dsa(:,col)));
    [Inn] = find(~isnan(dsa(:,col)));
    if ~isempty(Inn)
        dsa(inds_kk,col) = mean(dsa(Inn,col));
    end
end

delta_SA_ref = dsa;

end

%##########################################################################

function delta_SA_ref = gsw_dsa_add_barrier(dsa,long,lat,longs_ref,lats_ref,dlongs_ref,dlats_ref)

% gsw_dsa_add_barrier
%==========================================================================
%
% USAGE:
%  delta_SA_ref = gsw_dsa_add_barrier(dsa,long,lat,longs_ref,lats_ref,dlongs_ref,dlats_ref)
%
% DESCRIPTION:
%  Adds a barrier through Central America (Panama) and then averages
%  over the appropriate side of the barrier
%
% INPUT:
%  dsa         =  Absolute Salinity anomaly                                [ g/kg ]
%  long        =  Longitudes of data in decimal degrees east               [ 0 ... +360 ]
%  lat         =  Latitudes of data in decimal degrees north               [ -90 ... +90 ]
%  longs_ref   =  Longitudes of regular grid in decimal degrees east       [ 0 ... +360 ]
%  lats_ref    =  Latitudes of regular grid in decimal degrees north       [ -90 ... +90 ]
%  dlongs_ref  =  Longitude difference of regular grid in decimal degrees  [ deg longitude ]
%  dlats_ref   =  Latitude difference of regular grid in decimal degrees   [ deg latitude ]
%
% OUTPUT:
%  delta_SA_ref   =  Reference Absolute Salinity anomaly                   [ g/kg ]
%
% AUTHOR: 
%  David Jackett
%
% VERSION NUMBER: 3.0
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T. J., D. R. Jackett and F. J. Millero, 2010: An algorithm 
%   for estimating Absolute Salinity in the global ocean.  Submitted to 
%   Ocean Science, a preliminary version is available at Ocean Sci. Discuss.,
%   6, 215-242.  
%   http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf
%   and the computer software is available from http://www.TEOS-10.org
%
%==========================================================================

longs_pan = [260.0000 272.5900 276.5000 278.6500 280.7300 295.2170];

lats_pan  = [ 19.5500  13.9700   9.6000   8.1000   9.3300   0];

lats_lines0 = interp1(longs_pan,lats_pan,long);

lats_lines1 = interp1(longs_pan,lats_pan,longs_ref);
lats_lines2 = interp1(longs_pan,lats_pan,(longs_ref+dlongs_ref));

for k0 = 1:length(long)
    if lats_lines0(k0) <= lat(k0)
        above_line0 = 1;
    else
        above_line0 = 0;
    end
    if lats_lines1(k0) <= lats_ref(k0)
        above_line(1) = 1;
    else
        above_line(1) = 0;
    end
    if lats_lines1(k0) <= (lats_ref(k0) + dlats_ref)
        above_line(4) = 1;
    else
        above_line(4) = 0;
    end
    if lats_lines2(k0) <= lats_ref(k0)
        above_line(2) = 1;
    else
        above_line(2) = 0;
    end
    if lats_lines2(k0) <= (lats_ref(k0) + dlats_ref)
        above_line(3) = 1;
    else
        above_line(3) = 0;
    end
    inds = find(above_line ~= above_line0);      % indices of different sides of CA line
    dsa(inds,k0) = nan;
end

dsa_mean = mean(dsa); 
inds_nan = find(isnan(dsa_mean)); 
no_nan = length(inds_nan);

for kk = 1:no_nan
    col = inds_nan(kk);
    inds_kk = find(isnan(dsa(:,col)));
    [Inn] = find(~isnan(dsa(:,col)));
    if ~isempty(Inn)
     dsa(inds_kk,col) = mean(dsa(Inn,col));
    end
end

delta_SA_ref = dsa;

end
