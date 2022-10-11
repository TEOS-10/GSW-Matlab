function t90 = gsw_t90_from_t68(t68)

% gsw_t90_from_t68              ITS-90 temperature from IPTS-68 temperature
%==========================================================================
%
% USAGE:
%  t90 = gsw_t90_from_t68(t68)
%
% DESCRIPTION:
%  Converts IPTS-68 temperature to International Temperature Scale 1990 
%  (ITS-90) temperature.  This conversion should be applied to all in-situ
%  data collected between 1/1/1968 and 31/12/1989.
%
% INPUT:
%  t68  =  in-situ temperature  (IPTS-68)                         [ deg C ]
%
% OUTPUT:
%  t90  =  in-situ temperature  (ITS-90)                          [ deg C ]
%
% AUTHOR: 
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  International Temperature Scales of 1948, 1968 and 1990, an ICES
%   note, available from http://www.ices.dk/ocean/procedures/its.htm
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 2.1 and appendix A.1 of this TEOS-10 manual. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables
%--------------------------------------------------------------------------

if ~(nargin == 1)
   error('gsw_t90_from_t68:  Requires only one input argument')
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

t90 = t68.*0.999760057586179;                         % t90 = t68./1.00024;

end
