function t_i = gsw_t_interp(t,p,p_i)

% gsw_t_interp                         interpolation of in-situ temperature
%                                                          to p_i on a cast
%==========================================================================
%
% USAGE:
%  t_i = gsw_t_interp(t,p,p_i)
%
% DESCRIPTION:
%  This function interpolates vertical casts of values of in situ
%  temperature, t, to the arbitrary pressures p_i. 
%
%  The interpolation method uses sixteen PCHIPs (Piecewise Cubic Hermite 
%  Interpolating Polynomials), one for each of sixteen different linear 
%  combinations of Conservative Temperature and a scaled version of the 
%  independent variable.  Each of these sixteen PCHIPs use a scaled version
%  of the "bottle number" as the independent variable.  A final seventeenth 
%  PCHIP is used to relate the interpolated data back to pressure space 
%  (rather than "botttle number" space).  The interpolation method is 
%  described as the MR-PCHIP method in Barker and McDougall (2020). 
%
%  The scaling factor that we use for the independent variable is 0.33
%  times the maximum magnitude (over all data pairs) of the slope on the 
%  [CT - bottle_number] diagram.  With the oceanographic data that we have 
%  tested we have found that using multiplying factors between 0.1 and 1 
%  (instead of 0.33) also gave good interpolants, so the method is not 
%  sensitive to this choice of "scale_factor". 
%
%  The code uses Conservative Temperature (of TEOS-10) as the temperature 
%  variable in the interpolation procedure.  The output of the code is 
%  converted back to in-situ temperature.  The conversions between in-situ 
%  and Conservative Temperatures are done using the constant value of 
%  Absolute Salinity equal to the Standard Ocean Reference Salinity 
%  (35.16504 g/kg) which is found by calling gsw_SSO.  We have found that
%  the use of this salinity rather than the estimate salinity from 
%  MRST-PCHIP method in Barker and McDougall (2020) results in an 
%  temperature error less than 1 mK.
%
%  Any interpolated bottles that have pressures shallower than the 
%  shallowest observed bottle are set equal to the shallowest observed 
%  bottle.
%
%  Note that this interpolation scheme requires at least four observed
%  bottles on the cast.
%
% INPUT:
%  t   = in situ temperature (ITS-90)                             [ deg C ]
%  p   = sea pressure                                              [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  p_i = specific query points at which the interpolated SA_i 
%          and CT_i are required                                   [ dbar ]
%
%  p may have dimensions Mx1 or 1xN or MxN, where t is MxN.
%  p_i needs to be either a vector or a matrix and have dimensions M_ix1
%  or M_ixN.
%
% OUTPUT:
%  t_i = interpolated in situ temperature values at pressures p_i [ deg C ]
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (15th July, 2020)
%
% REFERENCES:
%  Barker, P.M., and T.J. McDougall, 2020: Two interpolation methods using 
%   multiply-rotated piecewise cubic hermite interpolating polynomials. 
%   J. Atmosph. Ocean. Tech., 37, pp. 605-619. 
%   doi: 10.1175/JTECH-D-19-0211.1. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

[pl,number_of_profiles] = size(t);
[mp,np] = size(p);
[interp_profile_length,np_i] = size(p_i);

if (pl < 4) & (np == 1)
    error('gsw_t_interp:  There must be at least 4 bottles') 
elseif (number_of_profiles == np) & (mp == 1)
    p = p(ones(1,pl), :);
elseif (pl == mp) & (np == 1)
    p = p(:,ones(1,number_of_profiles));
elseif (number_of_profiles == mp) & (np == 1)
    p = p.';
    p = p(ones(1,pl), :);
elseif (pl == np) & (mp == 1)
    p = p.';
    p = p(:,ones(1,number_of_profiles));
elseif (pl == np) & (number_of_profiles == mp)
    p = p.';
elseif (pl == mp) & (number_of_profiles == np)
    % ok
else
    error('gsw_t_interp: Inputs array dimensions arguments do not agree')
end

if interp_profile_length == 1 & np_i > 1
    p_i = p_i.';
    dp_i = diff(p_i);
    if any(dp_i) < 0
        warning('gsw_t_interp: interpolating pressure must be monotonic')
        return
    end
    [interp_profile_length,np_i] = size(p_i);
elseif interp_profile_length == number_of_profiles & np_i~= number_of_profiles & all(diff(p_i,1,2)) >= 0
    p_i = p_i.';
    [interp_profile_length,np_i] = size(p_i);
elseif any(diff(p_i,1,1)) < 0
    warning('gsw_t_interp: interpolating pressure must be monotonic')
    return
else
    % Data shape and interval are ok.
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

scale_factor = 0.33;

t_i = NaN(interp_profile_length, number_of_profiles);

for Iprofile = 1:number_of_profiles
    
    data_bottles = t(:,Iprofile) + p(:,Iprofile);
    [Inn] = find(~isnan(data_bottles));
    if length(Inn) < 2
        t_i(:,Iprofile) = NaN;
        continue
    end
    
    t_obs = t(Inn,Iprofile);
    p_obs = p(Inn,Iprofile);
    
    if np_i ~= 1
        p_i_tmp = p_i(:,Iprofile);
    else
        p_i_tmp = p_i;
    end
     
    CT_obs = gsw_CT_from_t(gsw_SSO*ones(size(t_obs)),t_obs,p_obs);
    
    CT_i = gsw_data_interp(CT_obs,p_obs,p_i_tmp,scale_factor);
    
    t_i(:,Iprofile) = gsw_t_from_CT(gsw_SSO*ones(size(CT_i)),CT_i,p_i_tmp);

end

end
