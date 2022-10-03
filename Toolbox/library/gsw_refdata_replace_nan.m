function [ref_data_no_nan,weights_no_nan] = gsw_refdata_replace_nan(ref_data,weights)

% gsw_refdata_replace_nan
%==========================================================================
%
% USAGE:
%  ref_data_no_nan = gsw_refdata_replace_nan(ref_data,weights)
%
% DESCRIPTION:
%  Replaces NaN's with nanmean of the 4 adjacent neighbours
%
% INPUT:
%  ref_data  =  reference data of the 4 adjacent neighbours  [ data units ]
%
% Optional
%  weights   =  weights
%
% OUTPUT:
%  ref_data_no_nan  =  nanmean of the 4 adjacent neighbours  [ data units ]
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

[ref_data_nanmean,weights_nan_mean] = gsw_weighted_nanmean(ref_data,weights);

ref_data_nanmean(2,:) = ref_data_nanmean;
ref_data_nanmean(3:4,:) = ref_data_nanmean;
nan_data = isnan(ref_data);
[Inan] = find(isnan(ref_data));
ref_data_mean_nan = nan_data(Inan).*ref_data_nanmean(Inan);
ref_data(Inan) = ref_data_mean_nan;

if nargin > 1
    weights_nan_mean(2,:) = weights_nan_mean;
    weights_nan_mean(3:4,:) = weights_nan_mean;
    weights_mean_nan = nan_data(Inan).*weights_nan_mean(Inan);
    weights(Inan) = weights_mean_nan;
    
    % Apply weights
    weighted_ref_data = ref_data.*(weights) + ref_data_nanmean.*(1-weights);
    weights_no_nan = weights;
else
    weighted_ref_data = ref_data;
end

ref_data_no_nan = weighted_ref_data;

end