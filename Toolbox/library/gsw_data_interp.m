function data_i = gsw_data_interp(data,p,p_i,scale_factor)

% gsw_data_interp                                 data interpolation to p_i 
%==========================================================================
%
% USAGE:
%  data_i = gsw_data_interp(data,p,p_i,scale_factor)
%
% DESCRIPTION:
%  This function interpolates vertical casts of values of data to the
%  arbitrary pressures p_i.  
%
%  The interpolation method uses sixteen PCHIPs (Piecewise Cubic Hermite 
%  Interpolating Polynomials), one for each of sixteen different linear 
%  combinations of the data and a scaled version of the independent variable.  
%  Each of these sixteen PCHIPs use a scaled version of the "bottle number"
%  as the independent variable.  A final seventeenth PCHIP is used to 
%  relate the interpolated data back to p space (rather than "botttle 
%  number" space).  The interpolation method is described as the MR-PCHIP 
%  method in Barker and McDougall (2020). 
%
%  The scaling factor that we use for the independent variable is 0.33
%  times the maximum magnitude (over all data pairs) of the slope on the 
%  [data - bottle_number] diagram.  With the oceanographic data that we 
%  have tested we have found that using multiplying factors between 0.1 and
%  1 (instead of 0.33) also gave good interpolants, so the method is not 
%  sensitive to this choice of "scale_factor". 
%
%  Any interpolated bottles that have pressures shallower than the 
%  shallowest observed bottle are set equal to the shallowest observed 
%  bottle.
%
%  Note that this interpolation scheme requires at least four observed
%  bottles on the cast.
%
% INPUT:
%  data  =  data                                                      [ ? ]
%  p     =  sea pressure                                           [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  p_i  =  specific query points at which the interpolated tracer_i 
%          are required                                            [ dbar ]
%
% Optional:
%  scale_factor = scaling factor of the data in the world ocean. 
%                  The default value is 0.33.                  [ unitless ]
%
%  data and p must have the same dimensions, MxN.
%  p_i needs to be either a vector or a matrix and have dimensions M_ix1
%  or M_ixN.
%  scale_factor, if specified, needs to ba a scalar.
%
% OUTPUT:
%  data_i = interpolated data values at pressures p_i                 [ ? ]
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% References
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

[dummy,number_of_profiles] = size(data);
[interp_profile_length,np_i] = size(p_i);

if ~exist('scale_factor','var')
    scale_factor = 0.33;
elseif ~isscalar(scale_factor)
    error('gsw_data_interp: scale_factor must be a scalar')
else
    % scale_factor is ok.
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

data_i = NaN(interp_profile_length, number_of_profiles);

sin_pi_on_16  = 1.950903220161283e-1;   %sin_pi_on_16 = sin(pi./16)
cos_pi_on_16  = 9.807852804032304e-1;   %cos_pi_on_16 = cos(pi./16)
sin_pi_on_8   = 3.826834323650898e-1;   %sin_pi_on_8 = sin(pi./8)
cos_pi_on_8   = 9.238795325112867e-1;   %cos_pi_on_8 = cos(pi./8)
sin_3pi_on_16 = 5.555702330196022e-1;   %sin_3pi_on_16 = sin(3pi./16)
cos_3pi_on_16 = 8.314696123025452e-1;   %cos_3pi_on_16 = cos(3pi./16)
sin_pi_on_4   = 7.071067811865475e-1;   %sin_pi_on_8 = sin(pi./4)
cos_pi_on_4   = 7.071067811865476e-1;   %cos_pi_on_8 = cos(pi./4)
sin_5pi_on_16 = 8.314696123025452e-1;   %sin_pi_on_8 = sin(5pi./16)
cos_5pi_on_16 = 5.555702330196023e-1;   %cos_pi_on_8 = cos(5pi./16)
sin_3pi_on_8  = 9.238795325112867e-1;   %sin_pi_on_8 = sin(3pi./8)
cos_3pi_on_8  = 3.826834323650898e-1;   %cos_pi_on_8 = cos(3pi./8)
sin_7pi_on_16 = 9.807852804032304e-1;   %sin_pi_on_8 = sin(7pi./16)
cos_7pi_on_16 = 1.950903220161283e-1;   %cos_pi_on_8 = cos(7pi./16)

for Iprofile = 1:number_of_profiles
    
    data_bottles = data(:,Iprofile) + p(:,Iprofile);
    [Inn] = find(~isnan(data_bottles));
    if length(Inn) < 2
        data_i(:,Iprofile) = NaN;
        continue
    end
        
    data_obs = data(Inn,Iprofile);
    p_obs = p(Inn,Iprofile);
    if np_i > 1
        p_i_tmp = p_i(:,Iprofile);
    else
        p_i_tmp = p_i(:);
    end
        
    pl = length(p_obs);
    Ishallow = 1:(pl-1);
    Ideep = 2:pl;
    dp_tmp = p_obs(Ideep) - p_obs(Ishallow);
    if any(dp_tmp <= 0)
        warning('gsw_data_interp: pressure must be monotonic')
        [p_sort,Ipsort] = sort(p_obs);
        data_sort = data_obs(Ipsort);
        [p_obs, dummy, Iunique] = unique(p_sort);
        data_obs = accumarray(Iunique,data_sort,[],@mean); % This averages observations at the same pressure 
        pl = length(p_obs);
        Ishallow = 1:(pl-1);
        Ideep = 2:pl;
    end
        
    independent_variable = [0:pl-1].';
    
    [Itointerp] = find(p_i_tmp >= min(p_obs) & p_i_tmp <= max(p_obs));
    independent_variable_i = gsw_pchip_interp(independent_variable,p_obs,p_i_tmp(Itointerp));
    
    ddata_tmp = data_obs(Ideep) - data_obs(Ishallow);
    factor = scale_factor.*max(abs(ddata_tmp));
    
    dummy = [1:pl].';
    tor = factor.*dummy;
    
    v1_tmp = data_obs;
    q1_tmp = tor;
    [v1_i,q1_i] = gsw_pchip_interp_SA_CT(v1_tmp,q1_tmp,independent_variable,independent_variable_i);
    
    v2_tmp = tor.*sin_pi_on_16 + data_obs.*cos_pi_on_16;
    q2_tmp = tor.*cos_pi_on_16 - data_obs.*sin_pi_on_16;
    [v2_i,q2_i] = gsw_pchip_interp_SA_CT(v2_tmp,q2_tmp,independent_variable,independent_variable_i);
        
    v3_tmp = tor.*sin_pi_on_8 + data_obs.*cos_pi_on_8;
    q3_tmp = tor.*cos_pi_on_8 - data_obs.*sin_pi_on_8;
    [v3_i,q3_i] = gsw_pchip_interp_SA_CT(v3_tmp,q3_tmp,independent_variable,independent_variable_i);
        
    v4_tmp = tor.*sin_3pi_on_16 + data_obs.*cos_3pi_on_16;
    q4_tmp = tor.*cos_3pi_on_16 - data_obs.*sin_3pi_on_16;
    [v4_i,q4_i] = gsw_pchip_interp_SA_CT(v4_tmp,q4_tmp,independent_variable,independent_variable_i);
    
    v5_tmp = tor.*sin_pi_on_4 + data_obs.*cos_pi_on_4;
    q5_tmp = tor.*cos_pi_on_4 - data_obs.*sin_pi_on_4;
    [v5_i,q5_i] = gsw_pchip_interp_SA_CT(v5_tmp,q5_tmp,independent_variable,independent_variable_i);
        
    v6_tmp = tor.*sin_5pi_on_16 + data_obs.*cos_5pi_on_16;
    q6_tmp = tor.*cos_5pi_on_16 - data_obs.*sin_5pi_on_16;
    [v6_i,q6_i] = gsw_pchip_interp_SA_CT(v6_tmp,q6_tmp,independent_variable,independent_variable_i);
        
    v7_tmp = tor.*sin_3pi_on_8 + data_obs.*cos_3pi_on_8;
    q7_tmp = tor.*cos_3pi_on_8 - data_obs.*sin_3pi_on_8;
    [v7_i,q7_i] = gsw_pchip_interp_SA_CT(v7_tmp,q7_tmp,independent_variable,independent_variable_i);
      
    v8_tmp = tor.*sin_7pi_on_16 + data_obs.*cos_7pi_on_16;
    q8_tmp = tor.*cos_7pi_on_16 - data_obs.*sin_7pi_on_16;
    [v8_i,q8_i] = gsw_pchip_interp_SA_CT(v8_tmp,q8_tmp,independent_variable,independent_variable_i);
    
    data_i_1 = v1_i;

    data_i_2 = -q2_i.*sin_pi_on_16 + v2_i.*cos_pi_on_16;
    
    data_i_3 = -q3_i.*sin_pi_on_8 + v3_i.*cos_pi_on_8;
    
    data_i_4 = -q4_i.*sin_3pi_on_16 + v4_i.*cos_3pi_on_16;
    
    data_i_5 = -q5_i.*sin_pi_on_4 + v5_i.*cos_pi_on_4;
    
    data_i_6 = -q6_i.*sin_5pi_on_16 + v6_i.*cos_5pi_on_16;
    
    data_i_7 = -q7_i.*sin_3pi_on_8 + v7_i.*cos_3pi_on_8;
    
    data_i_8 = -q8_i.*sin_7pi_on_16 + v8_i.*cos_7pi_on_16;
    
    data_i_sum = 0.125.*(data_i_1 + data_i_2 + data_i_3 + data_i_4 + data_i_5 + data_i_6 + data_i_7 + data_i_8);

    [min_p_obs, Imin_p_obs] = min(p_obs);
    if min_p_obs ~= 0       
        [Isurface] = find(p_i_tmp < min_p_obs);
        data_i(Isurface,Iprofile) = data_obs(Imin_p_obs);
        data_i(Itointerp,Iprofile) = data_i_sum;
    else 
        data_i(Itointerp,Iprofile) = data_i_sum;        
    end
      
end

end
