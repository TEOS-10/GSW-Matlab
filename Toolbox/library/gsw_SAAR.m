function [SAAR, in_ocean] = gsw_SAAR(p,long,lat)

% gsw_SAAR       Absolute Salinity Anomaly Ratio (excluding the Baltic Sea)
%==========================================================================
%
% USAGE:  
%  [SAAR, in_ocean] = gsw_SAAR(p,long,lat)
%
% DESCRIPTION:
%  Calculates the Absolute Salinity Anomaly Ratio, SAAR, in the open ocean
%  by spatially interpolating the global reference data set of SAAR to the
%  location of the seawater sample.  
% 
%  This function uses version 3.0 of the SAAR look up table. 
%
%  The Absolute Salinity Anomaly Ratio in the Baltic Sea is evaluated 
%  separately, since it is a function of Practical Salinity, not of space. 
%  The present function returns a SAAR of zero for data in the Baltic Sea. 
%  The correct way of calculating Absolute Salinity in the Baltic Sea is by 
%  calling gsw_SA_from_SP.  
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
%  SAAR      =  Absolute Salinity Anomaly Ratio                [ unitless ]
%  in_ocean  =  0, if long and lat are a long way from the ocean 
%            =  1, if long and lat are in the ocean
%  Note. This flag is only set when the observation is well and truly on
%    dry land; often the warning flag is not set until one is several 
%    hundred kilometres inland from the coast. 
%
% AUTHOR: 
%  David Jackett                                       [ help@teos-10.org ]
%
% MODIFIED:
%  Paul Barker and Trevor McDougall
%
% VERSION NUMBER: 3.01 (31st May, 2011)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., D.R. Jackett and F.J. Millero, 2010: An algorithm 
%   for estimating Absolute Salinity in the global ocean.  Submitted to 
%   Ocean Science. A preliminary version is available at Ocean Sci. 
%   Discuss., 6, 215-242.  
%   http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_SAAR:  Requires three inputs')
end %if

[mp,np] = size(p);
[mla,nla] = size(lat);
[mlo,nlo] = size(long);

if (mp ~= mla) | (mp ~=mlo) | (np ~= nla) | (np ~= nlo)
    error('gsw_SAAR: Inputs need be of the same size')
end %if

if ~isempty(find(p < -1.5))
    error('gsw_SAAR: pressure needs to be positive')
end

[Ipn] = find(p < 0 & p > -1.5);
if ~isempty(Ipn)
    p(Ipn) = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation (extracting from a look up table)
%--------------------------------------------------------------------------

gsw_data = 'gsw_data_v3_0.mat';

gsw_data_file = which(gsw_data);

load (gsw_data_file,'SAAR_ref','lats_ref','longs_ref','p_ref', ...
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
SAAR = nan(size(data_inds));
in_ocean = ones(size(SAAR));

for k = 1:nz-1
    
    inds_k = find(indsz0 == k);
    nk = length(inds_k);
    
    if nk>0
        nksum = nksum+nk;
        indsx = indsx0(inds_k);
        indsy = indsy0(inds_k);
        indsz = k*ones(size(indsx));
        inds_di = find(data_inds == k);             % level k interpolation
        saar = nan(4,n0);
        inds1 = sub2ind([nz,ny,nx], indsz, indsy, indsx);
        saar(1,inds_k) = SAAR_ref(inds1);
        inds2 = sub2ind([nz,ny,nx], indsz, indsy, indsx+1);
        saar(2,inds_k) = SAAR_ref(inds2);                   % inds0 + ny*nz
        inds3 = sub2ind([nz,ny,nx], indsz, indsy+1, indsx+1);
        saar(3,inds_k) = SAAR_ref(inds3);              % inds0 + ny*nz + nz
        inds4 = sub2ind([nz ny,nx], indsz, indsy+1, indsx);
        saar(4,inds_k) = SAAR_ref(inds4);                      % inds0 + nz
                       
        inds = find(260<=long(:) & long(:)<=295.217 & ...
            0<=lat(:) & lat(:)<=19.55 & indsz0(:)==k);
        if ~isempty(inds)
            saar(:,inds) = gsw_saar_add_barrier(saar(:,inds),long(inds), ...
                lat(inds),longs_ref(indsx0(inds)),lats_ref(indsy0(inds)),dlongs_ref,dlats_ref);
        end
        
        inds = find(isnan(sum(saar))' & indsz0==k);
        if ~isempty(inds)
            saar(:,inds) = gsw_saar_add_mean(saar(:,inds));
        end
        
        sa_upper(inds_di) = (1-s1(inds_di)).*(saar(1,inds_k)' + ...
            r1(inds_di).*(saar(2,inds_k)'-saar(1,inds_k)')) + ...
            s1(inds_di).*(saar(4,inds_k)' + ...
            r1(inds_di).*(saar(3,inds_k)'-saar(4,inds_k)'));  % level k+1 interpolation
                
        saar = nan(4,n0);
        inds1 = sub2ind([nz,ny,nx], indsz+1, indsy, indsx);
        saar(1,inds_k) = SAAR_ref(inds1);
        inds2 = sub2ind([nz,ny,nx], indsz+1, indsy, indsx+1);
        saar(2,inds_k) = SAAR_ref(inds2);                   % inds1 + ny*nz
        inds3 = sub2ind([nz,ny,nx], indsz+1, indsy+1, indsx+1);
        saar(3,inds_k) = SAAR_ref(inds3);              % inds1 + ny*nz + nz
        inds4 = sub2ind([nz ny,nx], indsz+1, indsy+1, indsx);
        saar(4,inds_k) = SAAR_ref(inds4);                      % inds1 + nz
                
        inds = find(260<=long(:) & long(:)<=295.217 & ...
            0<=lat(:) & lat(:)<=19.55 & indsz0(:)==k);
        if ~isempty(inds)
            saar(:,inds) = gsw_saar_add_barrier(saar(:,inds),long(inds), ...
                lat(inds),longs_ref(indsx0(inds)),lats_ref(indsy0(inds)),dlongs_ref,dlats_ref);
        end
        
        inds = find(isnan(sum(saar))' & indsz0==k);
        
        if ~isempty(inds)
            saar(:,inds) = gsw_saar_add_mean(saar(:,inds));
        end
        
        sa_lower(inds_di) = (1-s1(inds_di)).*(saar(1,inds_k)' + ...
            r1(inds_di).*(saar(2,inds_k)'-saar(1,inds_k)')) + ...
            s1(inds_di).*(saar(4,inds_k)' + ...
            r1(inds_di).*(saar(3,inds_k)'-saar(4,inds_k)'));
        
        inds_different = find(isfinite(sa_upper(inds_di)) & isnan(sa_lower(inds_di)));
        
        if ~isempty(inds_different)
            sa_lower(inds_di(inds_different)) = sa_upper(inds_di(inds_different));
        end
        
        SAAR(inds_di) = sa_upper(inds_di) + t1(inds_di).*(sa_lower(inds_di) - sa_upper(inds_di));
        
    else
        no_levels_missing = no_levels_missing + 1;
    end
end

inds = find(~isfinite(SAAR)); 
SAAR(inds) = 0;

in_ocean(inds) = 0;

end

%##########################################################################

function SAAR = gsw_saar_add_mean(saar)

% gsw_saar_add_mean
%==========================================================================
%
% USAGE:
%  SAAR = gsw_saar_add_mean(saar)
%
% DESCRIPTION:
%  Replaces NaN's with nanmean of the 4 adjacent neighbours
%
% INPUT:
%  saar  =  Absolute Salinity Anomaly Ratio of the 4 adjacent neighbours  
%                                                              [ unitless ]
%
% OUTPUT:
%  SAAR  =  nanmean of the 4 adjacent neighbours               [ unitless ]
%
% AUTHOR: 
%  David Jackett
%
% MODIFIED:
%  Paul Barker and Trevor McDougall
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
%   Ocean Science, a preliminary version is available at Ocean Sci. 
%   Discuss., 6, 215-242.  
%   http://www.ocean-sci-discuss.net/6/215/2009/osd-6-215-2009-print.pdf
%   and the computer software is available from http://www.TEOS-10.org
%
%==========================================================================

saar_mean = mean(saar); 
inds_nan = find(isnan(saar_mean)); 
no_nan = length(inds_nan);

for kk = 1:no_nan
    col = inds_nan(kk);
    inds_kk = find(isnan(saar(:,col)));
    [Inn] = find(~isnan(saar(:,col)));
    if ~isempty(Inn)
        saar(inds_kk,col) = mean(saar(Inn,col));
    end
end

SAAR = saar;

end

%##########################################################################

function SAAR = gsw_saar_add_barrier(saar,long,lat,longs_ref,lats_ref,dlongs_ref,dlats_ref)

% gsw_saar_add_barrier
%==========================================================================
%
% USAGE:
%  SAAR = gsw_saar_add_barrier(saar,long,lat,longs_ref,lats_ref,dlongs_ref,dlats_ref)
%
% DESCRIPTION:
%  Adds a barrier through Central America (Panama) and then averages
%  over the appropriate side of the barrier
%
% INPUT:
%  saar        =  Absolute Salinity Anomaly Ratio                          [ unitless ]
%  long        =  Longitudes of data in decimal degrees east               [ 0 ... +360 ]
%  lat         =  Latitudes of data in decimal degrees north               [ -90 ... +90 ]
%  longs_ref   =  Longitudes of regular grid in decimal degrees east       [ 0 ... +360 ]
%  lats_ref    =  Latitudes of regular grid in decimal degrees north       [ -90 ... +90 ]
%  dlongs_ref  =  Longitude difference of regular grid in decimal degrees  [ deg longitude ]
%  dlats_ref   =  Latitude difference of regular grid in decimal degrees   [ deg latitude ]
%
% OUTPUT:
%  SAAR        =  Absolute Salinity Anomaly Ratio                          [ unitless ]
%
% AUTHOR: 
%  David Jackett
%
% MODIFIED:
%  Paul Barker and Trevor McDougall
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
%   Ocean Science, a preliminary version is available at Ocean Sci. 
%   Discuss., 6, 215-242.  
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
    saar(inds,k0) = nan;
end

saar_mean = mean(saar); 
inds_nan = find(isnan(saar_mean)); 
no_nan = length(inds_nan);

for kk = 1:no_nan
    col = inds_nan(kk);
    inds_kk = find(isnan(saar(:,col)));
    [Inn] = find(~isnan(saar(:,col)));
    if ~isempty(Inn)
     saar(inds_kk,col) = mean(saar(Inn,col));
    end
end

SAAR = saar;

end
