function t90 = gsw_t90_from_t48(t48)

% gsw_t90_from_t48              ITS-90 temperature from IPTS-48 temperature
%==========================================================================
%
% USAGE:
%  t90 = gsw_t90_from_t48(t48)
%
% DESCRIPTION:
%  Calculates IPTS-48 temperature to International Temperature Scale 1990 
%  (ITS-90) temperature.  This conversion should be applied to all in-situ
%  data collected prior to 31/12/1967.
%  
% INPUT:
%  t48  =  in-situ temperature  (IPTS-48)                         [ deg C ]
%
% OUTPUT:
%  t90  =  in-situ temperature  (ITS-90)                          [ deg C ]
%
% AUTHOR: 
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.02 (16th November, 2012)
%  This function is unchanged from version 2.0 (24th September, 2010).
%
% REFERENCES:
%  International Temperature Scales of 1948, 1968 and 1990, an ICES
%   note, available from http://www.ices.dk/ocean/procedures/its.htm
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables
%--------------------------------------------------------------------------

if ~(nargin == 1)
   error('gsw_t90_from_t48:  Requires only one input argument')
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

t90 = (t48 - (4.4e-6)*t48.*(100 - t48))./1.00024;

end