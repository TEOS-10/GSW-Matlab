function tracer_i = gsw_tracer_interp(tracer,p,p_i,scale_factor)

% gsw_tracer_interp              interpolation of a tracer to p_i on a cast
%==========================================================================
%
% USAGE:
%  tracer_i = gsw_tracer_interp(tracer,p,p_i,{scale_factor})
%
% DESCRIPTION:
%  This function interpolates vertical casts of values of tracer data to
%  the arbitrary pressures p_i. 
%
%  The interpolation method uses sixteen PCHIPs (Piecewise Cubic Hermite 
%  Interpolating Polynomials), one for each of sixteen different linear 
%  combinations of the tracer data and a scaled version of the independent 
%  variable.  Each of these sixteen PCHIPs use a scaled version of the 
%  "bottle number" as the independent variable.  A final seventeenth 
%  PCHIP is used to relate the interpolated data back to pressure space 
%  (rather than "botttle number" space).  The interpolation method is 
%  described as the MR-PCHIP method in Barker and McDougall (2020). 
%
%  When the tracer is in situ temperature we have found a suitable value 
%  for the scale_factor is 0.33, so that the final scaling factor is 0.33
%  times the maximum magnitude (over all data pairs) of the slope on the 
%  [tracer - bottle_number] diagram.  We expect 0.33 will be a suitable 
%  scale_factor for other tracer data. 
%
%  When values of Conservative Temperature are also available along with
%  each value of tracer data, a better interpolation code than the present
%  one is gsw_tracer_CT_interp(tracer,CT,p,p_i,factor).  
%
%  Any interpolated bottles that have pressures shallower than the 
%  shallowest observed bottle are set equal to the shallowest observed 
%  bottle.
%
%  Note that this interpolation scheme requires at least four observed
%  bottles on the cast.
%
% INPUT:
%  t    =  tracer                                                     [ ? ]
%  p    =  sea pressure                                            [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  p_i  =  specific query points at which the interpolated tracer_i 
%          are required                                            [ dbar ]
%
% OPTIONAL:
%  scale_factor = scaling factor of the maximum magnitude of the slope on 
%                 the [tracer - bottle_number] diagram.  The default value
%                 is 0.33.                                     [ unitless ]
%
%  p may have dimensions Mx1 or 1xN or MxN, where tracer is MxN.
%  p_i needs to be either a vector or a matrix and have dimensions M_ix1
%  or M_ixN.
%  scale_factor, if provided, must be a scalar.
%
% OUTPUT:
%  tracer_i = interpolated tracer values at pressures p_i             [ ? ]
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th June, 2020)
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

[pl,number_of_profiles] = size(tracer);
[mp,np] = size(p);
[interp_profile_length,np_i] = size(p_i);

if (pl < 4) & (np == 1)
    error('gsw_tracer_interp:  There must be at least 4 bottles')
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
    error('gsw_tracer_interp: Inputs array dimensions arguments do not agree')
end

if interp_profile_length == 1 & np_i > 1
    p_i = p_i.';
    dp_i = diff(p_i);
    if any(dp_i) < 0
        warning('gsw_tracer_interp: interpolating pressure must be monotonic')
        return
    end
    [interp_profile_length,np_i] = size(p_i);
elseif interp_profile_length == number_of_profiles & np_i~= number_of_profiles & all(diff(p_i,1,2)) >= 0
    p_i = p_i.';
    [interp_profile_length,np_i] = size(p_i);
elseif any(diff(p_i,1,1)) < 0
    warning('gsw_tracer_interp: interpolating pressure must be monotonic')
    return
else
    % Data shape and interval are ok.
end

if ~exist('scale_factor','var')
    scale_factor = 0.33;
else
    if ~isscalar(scale_factor)
        error('gsw_tracer_interp: scale_factor must be a scalar')
    end
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

tracer_i = NaN(interp_profile_length, number_of_profiles);

for Iprofile = 1:number_of_profiles
    
    data_bottles = tracer(:,Iprofile) + p(:,Iprofile);
    [Inn] = find(~isnan(data_bottles));
    if length(Inn) < 2
        tracer_i(:,Iprofile) = NaN;
        continue
    end
    
    tracer_obs = tracer(Inn,Iprofile);
    p_obs = p(Inn,Iprofile);
    
    if np_i ~= 1
        p_i_tmp = p_i(:,Iprofile);
    else
        p_i_tmp = p_i;
    end
       
     tracer_i(:,Iprofile) = gsw_data_interp(tracer_obs,p_obs,p_i_tmp,scale_factor);
    
end

end
