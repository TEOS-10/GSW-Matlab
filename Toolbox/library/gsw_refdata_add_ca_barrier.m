function [ref_data_ca_barrier, weights_ca_barrier] = gsw_refdata_add_ca_barrier(ref_data,long,lat,longs_ref,lats_ref,dlongs_ref,dlats_ref,weights)

% gsw_refdata_add_ca_barrier
%==========================================================================
%
% USAGE:
%  [ref_data_ca_barrier, weights_ca_barrier] = gsw_refdata_add_ca_barrier...
%               (ref_data,long,lat,longs_ref,lats_ref,dlongs_ref,dlats_ref)
%
% DESCRIPTION:
%  Adds a barrier through Central America (Panama) and then averages
%  over the appropriate side of the barrier
%
% INPUT:
%  ref_data    =  reference data                                           [ data units ]
%  long        =  Longitudes of data in decimal degrees east               [ 0 ... +360 ]
%  lat         =  Latitudes of data in decimal degrees north               [ -90 ... +90 ]
%  longs_ref   =  Longitudes of regular grid in decimal degrees east       [ 0 ... +360 ]
%  lats_ref    =  Latitudes of regular grid in decimal degrees north       [ -90 ... +90 ]
%  dlongs_ref  =  Longitude difference of regular grid in decimal degrees  [ deg longitude ]
%  dlats_ref   =  Latitude difference of regular grid in decimal degrees   [ deg latitude ]
%
% Optional
%  weights     =  weights                                                  [ 0 ... 1 ]
%
% OUTPUT:
%  ref_data_ca_barrier  =  reference data averaged on the appropiate side 
%                          of Central America                              [ data units ]
%  weights_ca_barrier  =  weights averaged on the appropiate side 
%                          of Central America                              [ 0 ... 1 ]
%
% AUTHOR: 
%  David Jackett
%
% MODIFIED:
%  Paul Barker and Trevor McDougall
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
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
    ref_data(above_line ~= above_line0,k0) = nan;     % indices of different sides of CA line
    weights(above_line ~= above_line0,k0) = nan;
end

if any(isnan(ref_data(:)))
    [ref_data_ca_barrier,weights_ca_barrier] = gsw_refdata_replace_nan(ref_data,weights);
else
    ref_data_ca_barrier = ref_data;
    weights_ca_barrier = weights;
end

end
