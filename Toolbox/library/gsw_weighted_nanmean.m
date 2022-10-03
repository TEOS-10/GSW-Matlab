function [weighted_mean, weighted_weights] = gsw_weighted_nanmean(data,weights)

% gsw_weighted_nanmean                      weighted mean of non-NaN values
%==========================================================================
%
% USAGE:
%  [weighted_mean, weighted_weights] = gsw_weighted_nanmean(data,{weights})
%
% DESCRIPTION:
%  This function returns the weighted mean of the data, where the weights
%  are given by weights and a coresponding weighted weight for the weighted
%  mean.  If the weights are not supplied then the weights are considerd to
%  be all equal. This programme treats NaN's as missing values.   
%
% INPUT:
%  data    = data values
%
% Optional
%  weights = weights of the data
%
% OUTPUT:
%  weighted_mean     =  weighted mean of the data
%  weighted_weights  =  weighted weights of the weights
%
% AUTHOR: 
%  Paul Barker and Trevor McDougall, based on nanmean from Matlab.
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

if nargin == 1
    weights = ones(size(data));
end

weighted_data = weights.*data;

Inan = isnan(weighted_data); 
weighted_data(Inan) = 0;
weights(Inan) = 0;

sum_weights = sum(weights);     
sum_weights(sum_weights == 0) = NaN;  

weighted_mean = sum(weighted_data)./sum_weights; 

weights_nn = sum(~Inan);  % non-NaN weights
weights_nn(weights_nn == 0) = NaN;
weighted_weights = sum_weights./weights_nn;

end
