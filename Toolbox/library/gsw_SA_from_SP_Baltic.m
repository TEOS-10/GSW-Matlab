function SA_baltic = gsw_SA_from_SP_Baltic(SP,long,lat)

% gsw_SA_from_SP_Baltic      Calculates Absolute Salinity in the Baltic Sea 
%==========================================================================
%
% USAGE:  
%  SA_baltic = gsw_SA_from_SP_Baltic(SP,long,lat)
%
% DESCRIPTION:
%  Calculates Absolute Salinity in the Baltic Sea, from Practical Salinity.
%  Since SP is non-negative by definition, this function changes any 
%  negative input values of SP to be zero.  
%  Note. This programme will only produce Absolute Salinity values for the
%  Baltic Sea.
%
% INPUT:
%   SP    =  Practical Salinity  (PSS-78)                      [ unitless ]
%   long  =  Longitude in decimal degrees east               [ 0 ... +360 ]    
%   lat   =  Latitude in decimal degrees north              [ -90 ... +90 ]  
%
% OUTPUT:
%   SA_baltic    =  Absolute Salinity in the Baltic Sea            [ g/kg ]
%
% AUTHOR: 
%  David Jackett, Trevor McDougall & Paul Barker [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (23rd July, 2010)
%
% REFERENCES:
%  Feistel, R., S. Weinreben, H. Wolf, S. Seitz, P. Spitzer, B. Adel, 
%   G. Nausch, B. Schneider and D. G. Wright, 2010: Density and Absolute 
%   Salinity of the Baltic Sea 2006-2009.  Ocean Science, 6, 3-24.
%   http://www.ocean-sci.net/6/3/2010/os-6-3-2010.pdf 
%
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

if ~(nargin == 3)
   error('gsw_SA_from_SP_Baltic.m:  Requires 3 inputs')
end %if

% These few lines ensure that SP is non-negative.
[I_neg_SP] = find(SP < 0);
if ~isempty(I_neg_SP)
    SP(I_neg_SP) = 0;
end

xb1 = 12.6; 
xb2 = 7; 
xb3 = 26; 
xb1a = 45; 
xb3a = 26;

yb1 = 50; 
yb2 = 59; 
yb3 = 69;

inds_baltic = find(xb2<long & long<xb1a & yb1<lat & lat<yb3);

SA_baltic = nan(size(SP));

if ~isempty(inds_baltic)
    xx_left = interp1([yb1,yb2,yb3],[xb1,xb2,xb3],lat(inds_baltic));
    xx_right = interp1([yb1,yb3],[xb1a,xb3a],lat(inds_baltic));
    inds_baltic1 = find(xx_left<=long(inds_baltic) & long(inds_baltic)<=xx_right);
    SA_baltic(inds_baltic(inds_baltic1)) = ((35.16504 - 0.087)/35)*SP(inds_baltic(inds_baltic1)) + 0.087; 
    SA_baltic = reshape(SA_baltic,size(long));
end

end
