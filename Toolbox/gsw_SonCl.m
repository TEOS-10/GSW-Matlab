function SonCl = gsw_SonCl

% gsw_SonCl                                          SP to Chlorinity ratio
%==========================================================================
%
% USAGE:
%  SonCl = gsw_SonCl
%
% DESCRIPTION:
%  The ratio of Practical Salinity, SP, to Chlorinity, 1.80655 kg/g for
%  Reference Seawater (Millero et al., 2008).  This is the ratio that was 
%  used by the JPOTS committee in their construction of the 1978 Practical
%  Salinity Scale (PSS-78) to convert between the laboratory measurements 
%  of seawater samples (which were measured in Chlorinity) to Practical 
%  Salinity.  
%
% OUTPUT:
%  SonCl  =  SP to Chlorinity ratio                           [ (g/kg)^-1 ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  Millero, F. J., R. Feistel, D. G. Wright, and T. J. McDougall, 2008: 
%   The composition of Standard Seawater and the definition of the 
%   Reference-Composition Salinity Scale, Deep-Sea Res. I, 55, 50-72. 
%    See section 5 below Eqn. (5.5) of this paper.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

SonCl = 1.80655;

end
