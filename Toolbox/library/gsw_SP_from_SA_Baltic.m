function SP_baltic = gsw_SP_from_SA_Baltic(SA,long,lat)

% gsw_SP_from_SA_Baltic    Calculates Practical Salinity for the Baltic Sea 
%==========================================================================
%
% USAGE:  
%  SP_baltic = gsw_SP_from_SA_Baltic(SA,long,lat)
%
% DESCRIPTION:
%  Calculates Practical Salinity for the Baltic Sea, from a value computed
%  analytically from Absolute Salinity.
%  Note. This programme will only produce Practical Salinty values for the
%    Baltic Sea.
%
% INPUT:
%  SA    =  Absolute Salinity in the Baltic Sea                 [ g kg^-1 ]
%  long  =  Longitude in decimal degress east                [ 0 ... +360 ]    
%  lat   =  Latitude in decimal degress north               [ -90 ... +90 ]  
%
% OUTPUT:
%  SP_baltic  =  Practical Salinity                            [ unitless ]
%
% AUTHOR: 
%  David Jackett, Trevor McDougall & Paul Barker       [ help@teos-10.org ]
%
% VERSION NUMBER: 3.01 (29th March, 2011) 
%
% REFERENCES:
%  Feistel, R., S. Weinreben, H. Wolf, S. Seitz, P. Spitzer, B. Adel, 
%   G. Nausch, B. Schneider and D. G. Wright, 2010c: Density and Absolute 
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
   error('gsw_SP_from_SA_Baltic.m:  Requires 3 inputs')
end %if

xb1 = 12.6; 
xb2 = 7; 
xb3 = 26; 
xb1a = 45; 
xb3a = 26;

yb1 = 50; 
yb2 = 59; 
yb3 = 69;

inds_baltic = find(xb2<long & long<xb1a & yb1<lat & lat<yb3);

SP_baltic = nan(size(SA));

if ~isempty(inds_baltic)
    xx_left = interp1([yb1,yb2,yb3],[xb1,xb2,xb3],lat(inds_baltic));
    xx_right = interp1([yb1,yb3],[xb1a,xb3a],lat(inds_baltic));
    inds_baltic1 = find(xx_left<=long(inds_baltic) & long(inds_baltic)<=xx_right);
    if ~isempty(inds_baltic1)
        SP_baltic(inds_baltic(inds_baltic1)) = (35/(35.16504 - 0.087))*(SA(inds_baltic(inds_baltic1)) - 0.087);
    end
%    SP_baltic = reshape(SP_baltic,size(long));
end

end

