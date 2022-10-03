function t_ice = gsw_t_from_pt0_ice(pt0_ice,p)

% gsw_t_from_pt0_ice                             in-situ temperature of ice
% =========================================================================
%
% USAGE:
%  t_ice = gsw_t_from_pt0_ice(pt0_ice,p)
%
% DESCRIPTION:
%  Calculates in-situ temperature from the potential temperature of ice Ih 
%  with reference pressure, p_ref, of 0 dbar (the surface), and the 
%  in-situ pressure.
%
% INPUT:
%  pt0_ice  =  potential temperature of ice Ih with reference pressure of 
%              zero dbar (ITS-90)                                 [ deg C ]
%  p        =  sea pressure                                        [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  pt0 & p must have the same dimensions.
%
% OUTPUT:
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix I of this TEOS-10 Manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
    error('gsw_t_from_pt0_ice:  Requires two inputs')
end 

[mt,nt] = size(pt0_ice);
[mp,np] = size(p);

if (mt ~= mp | nt ~= np )
    error('gsw_t_from_pt0_ice: Input arguments do not have the same dimensions')
end %if

if mt == 1
    pt0_ice = pt0_ice.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

p0 = zeros(size(pt0_ice));
t_ice = gsw_pt_from_t_ice(pt0_ice,p0,p);

if transposed
    t_ice = t_ice.';
end

end
