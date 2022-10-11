function entropy = gsw_entropy_t_exact(SA,t,p)

% gsw_entropy_from_t                           specific entropy of seawater
%==========================================================================
%
% This function has changed its name, it is now called 
% gsw_entropy_from_t.
%
% USAGE:
%  entropy  =  gsw_entropy_from_t(SA,t,p)
%
% DESCRIPTION:
%  Calculates specific entropy of seawater. 
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  entropy  =  specific entropy                                [ J/(kg*K) ]
%
% AUTHOR: 
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
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

warning('The function "entropy_t_exact" has changed to "gsw_entropy_from_t"');

entropy = gsw_entropy_from_t(SA,t,p);

end