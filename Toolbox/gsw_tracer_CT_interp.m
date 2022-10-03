function [tracer_i, CT_i] = gsw_tracer_CT_interp(tracer,CT,p,p_i,factor)

% gsw_tracer_CT_interp                          tracer and CT interpolation
%                                                          to p_i on a cast
%==========================================================================
%
% USAGE:
%  [tracer_i, CT_i] = gsw_tracer_CT_interp(tracer,CT,p,p_i,factor)
%
% DESCRIPTION:
%  This function interpolates vertical casts of values of a tracer
%  and Conservative Temperature to the arbitrary pressures p_i. The
%  interpolation method is designed to respect the shape of the tracer-CT 
%  diagram.  That is, the interpolated tracer_i and CT_i values look 
%  realistic when plotted on the tracer-CT diagram.  The interpolation 
%  method uses sixteen PCHIPs (Piecewise Cubic Hermite Interpolating 
%  Polynomials), one for each of sixteen different linear combinations of 
%  the tracer and CT input data.  Each of these sixteen PCHIPs use the 
%  "bottle number" as the independent variable.  A final seventeenth PCHIP
%  is used to relate the interpolated data back to pressure space (rather 
%  than "botttle number" space).  The interpolation method is is the 
%  MRST-PCHIP method described in Barker and McDougall (2020), with the
%  tracer data being used in place of Absoluate Salinity data.  
%
%  This function requires scaling the tracer and temperature data so that
%  the tracer-CT diagram reflects the relative variation of the tracer and 
%  temperature in the world ocean.  Specifically, "factor" should be chosen
%  to be the ratio of the global range of CT to that of the tracer variable
%  in the world ocean.  A list of suitable values of "factor" for various
%  tracers is given here.  
% 
%      TRACER              UNITS             FACTOR
%    Absolute Salinity     g/kg                9
%    dissolved oxygen        ?                 ?
%    AOU                     ?                 ?
%    silicic acid            ?                 ?
%    nitrate                 ?                 ?
%    phosphate               ?                 ?
%    carbon 14               ?                 ?
%    tritium                 ?                 ?
%    eastward velocity      m/s               100
%    westward velocity      m/s               100
%    
%  If an input value of "factor" is not given in the function call, it is 
%  set equal to 9.
%
%  Any interpolated bottles that have pressures shallower than the 
%  shallowest observed bottle are set equal to the shallowest observed 
%  bottle.
%
%  Note that this interpolation scheme requires at least four observed
%  bottles on the cast.
%
% INPUT:
%  tracer =  tracer                                                   [ ? ]
%  CT     =  Conservative Temperature (ITS-90)                    [ deg C ]
%  p      =  sea pressure                                          [ dbar ]
%              ( i.e. absolute pressure - 10.1325 dbar )
%  p_i    =  specific query points at which the interpolated SA_i 
%            and CT_i are required                                 [ dbar ]
%
% OPTIONAL INPUT:
%  factor = ratio between the ranges of Conservative Temperature and tracer
%           in the world ocean. The default value is 9.               [ ? ]
%
%  tracer & CT need to have the same dimensions.
%  p may have dimensions Mx1 or 1xN or MxN, where tracer & CT are MxN.
%  p_i needs to be either a vector or a matrix and have dimensions M_ix1
%  or M_ixN.
%  factor, if supplied, needs to be a scalar.
%
% OUTPUT:
%  tracer_i  =  interpolated tracer values at pressures p_i           [ ? ]
%  CT_i      =  interpolated CT values at pressures p_i           [ deg C ]
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
[mt,nt] = size(CT);
[mp,np] = size(p);
[interp_profile_length,np_i] = size(p_i);

if (pl ~= mt) | (number_of_profiles ~= nt)
    error('gsw_tracer_CT_interp: tracer and CT need to have the same dimensions')
end

if (pl < 4) & (np == 1)
    error('gsw_tracer_CT_interp:  There must be at least 4 bottles')
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
    error('gsw_tracer_CT_interp: Inputs array dimensions arguments do not agree')
end

if interp_profile_length == 1 & np_i > 1
    p_i = p_i.';
    dp_i = diff(p_i);
    if any(dp_i) < 0
        warning('gsw_tracer_CT_interp: interpolating pressure must be monotonic')
        return
    end
    [interp_profile_length,np_i] = size(p_i);
elseif interp_profile_length == number_of_profiles & np_i~= number_of_profiles & all(diff(p_i,1,2)) >= 0
    p_i = p_i.';
    [interp_profile_length,np_i] = size(p_i);
elseif any(diff(p_i,1,1)) < 0
    warning('gsw_tracer_CT_interp: interpolating pressure must be monotonic')
    return
else
    % Data shape and interval are ok.
end

if ~exist('factor','var')
    factor = 9;
else
    if ~isscalar(factor)
        error('gsw_tracer_CT_interp: factor must be a scalar')
    end
end
rec_factor = 1/factor;

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

tracer_i = NaN(interp_profile_length, number_of_profiles);
CT_i = tracer_i;

sin_pi_on_16 = 1.950903220161283e-1;  %sin_pi_on_16 = sin(pi./16)
cos_pi_on_16 = 9.807852804032304e-1;  %cos_pi_on_16 = cos(pi./16)
sin_pi_on_8 = 3.826834323650898e-1;   %sin_pi_on_8 = sin(pi./8)
cos_pi_on_8 = 9.238795325112867e-1;   %cos_pi_on_8 = cos(pi./8)
sin_3pi_on_16 = 5.555702330196022e-1; %sin_3pi_on_16 = sin(3pi./16)
cos_3pi_on_16 = 8.314696123025452e-1; %cos_3pi_on_16 = cos(3pi./16)
sin_pi_on_4 = 7.071067811865475e-1;   %sin_pi_on_8 = sin(pi./4)
cos_pi_on_4 = 7.071067811865476e-1;   %cos_pi_on_8 = cos(pi./4)
sin_5pi_on_16 = 8.314696123025452e-1; %sin_pi_on_8 = sin(5pi./16)
cos_5pi_on_16 = 5.555702330196023e-1; %cos_pi_on_8 = cos(5pi./16)
sin_3pi_on_8 = 9.238795325112867e-1;  %sin_pi_on_8 = sin(3pi./8)
cos_3pi_on_8 = 3.826834323650898e-1;  %cos_pi_on_8 = cos(3pi./8)
sin_7pi_on_16 = 9.807852804032304e-1; %sin_pi_on_8 = sin(7pi./16)
cos_7pi_on_16 = 1.950903220161283e-1; %cos_pi_on_8 = cos(7pi./16)

for Iprofile = 1:number_of_profiles
    
    data_bottles = tracer(:,Iprofile) + CT(:,Iprofile) + p(:,Iprofile);
    [Inn] = find(~isnan(data_bottles));
    if length(Inn) < 2
        tracer_i(:,Iprofile) = NaN;
        CT_i(:,Iprofile) = NaN;
        continue
    end
    
    tracer_obs = tracer(Inn,Iprofile);
    CT_obs = CT(Inn,Iprofile);
    p_obs = p(Inn,Iprofile);
    dummy = 1e-3.*round(1e3.*p_obs);
    p_obs = dummy;
    if np_i > 1
        p_i_tmp = p_i(:,Iprofile);
    else
        p_i_tmp = p_i(:);
    end
    dummy = 1e-3.*round(1e3.*p_i_tmp);
    p_i_tmp = dummy;
    
    pl = length(p_obs);
    Ishallow = 1:(pl-1);
    Ideep = 2:pl;
    dp_tmp = p_obs(Ideep) - p_obs(Ishallow);
    if any(dp_tmp <= 0)
        warning('gsw_tracer_CT_interp: pressure must be monotonic')
        [p_sort,Ipsort] = sort(p_obs);
        tracer_sort = tracer_obs(Ipsort);
        CT_sort = CT_obs(Ipsort);        
        [p_obs, dummy, Iunique] = unique(p_sort); 
        tracer_obs = accumarray(Iunique,tracer_sort,[],@mean); % This averages observations at the same pressure
        CT_obs = accumarray(Iunique,CT_sort,[],@mean);        
        pl = length(p_obs);
    end
    
    p_all = unique(sort([p_obs; p_i_tmp]));
    [Iobs_plus_interp] = find(p_all >= min(p_obs) & p_all <= max(p_obs));
    [Isurf_and_obs_plus_interp] = find(p_all <= max(p_obs));
    [dummy, Iout, I1] = intersect(p_i_tmp,p_all(Isurf_and_obs_plus_interp));
    [dummy, I2, I3] = intersect(p_obs,p_all(Iobs_plus_interp));
    
    clear independent_variable
    independent_variable(1:pl,1) = [0:pl-1];
    independent_variable_obs_plus_interp = gsw_pchip_interp(independent_variable,p_obs,p_all(Iobs_plus_interp));
           
    scaled_tracer_obs = factor.*tracer_obs;
    
    v1_tmp = CT_obs;
    q1_tmp = scaled_tracer_obs;
    [v1_i,q1_i] = gsw_pchip_interp_SA_CT(v1_tmp,q1_tmp,independent_variable,independent_variable_obs_plus_interp);
    
    v2_tmp = scaled_tracer_obs.*sin_pi_on_16 + CT_obs.*cos_pi_on_16;
    q2_tmp = scaled_tracer_obs.*cos_pi_on_16 - CT_obs.*sin_pi_on_16;
    [v2_i,q2_i] = gsw_pchip_interp_SA_CT(v2_tmp,q2_tmp,independent_variable,independent_variable_obs_plus_interp);
        
    v3_tmp = scaled_tracer_obs.*sin_pi_on_8 + CT_obs.*cos_pi_on_8;
    q3_tmp = scaled_tracer_obs.*cos_pi_on_8 - CT_obs.*sin_pi_on_8;
    [v3_i,q3_i] = gsw_pchip_interp_SA_CT(v3_tmp,q3_tmp,independent_variable,independent_variable_obs_plus_interp);
        
    v4_tmp = scaled_tracer_obs.*sin_3pi_on_16 + CT_obs.*cos_3pi_on_16;
    q4_tmp = scaled_tracer_obs.*cos_3pi_on_16 - CT_obs.*sin_3pi_on_16;
    [v4_i,q4_i] = gsw_pchip_interp_SA_CT(v4_tmp,q4_tmp,independent_variable,independent_variable_obs_plus_interp);
    
    v5_tmp = scaled_tracer_obs.*sin_pi_on_4 + CT_obs.*cos_pi_on_4;
    q5_tmp = scaled_tracer_obs.*cos_pi_on_4 - CT_obs.*sin_pi_on_4;
    [v5_i,q5_i] = gsw_pchip_interp_SA_CT(v5_tmp,q5_tmp,independent_variable,independent_variable_obs_plus_interp);
        
    v6_tmp = scaled_tracer_obs.*sin_5pi_on_16 + CT_obs.*cos_5pi_on_16;
    q6_tmp = scaled_tracer_obs.*cos_5pi_on_16 - CT_obs.*sin_5pi_on_16;
    [v6_i,q6_i] = gsw_pchip_interp_SA_CT(v6_tmp,q6_tmp,independent_variable,independent_variable_obs_plus_interp);
        
    v7_tmp = scaled_tracer_obs.*sin_3pi_on_8 + CT_obs.*cos_3pi_on_8;
    q7_tmp = scaled_tracer_obs.*cos_3pi_on_8 - CT_obs.*sin_3pi_on_8;
    [v7_i,q7_i] = gsw_pchip_interp_SA_CT(v7_tmp,q7_tmp,independent_variable,independent_variable_obs_plus_interp);
      
    v8_tmp = scaled_tracer_obs.*sin_7pi_on_16 + CT_obs.*cos_7pi_on_16;
    q8_tmp = scaled_tracer_obs.*cos_7pi_on_16 - CT_obs.*sin_7pi_on_16;
    [v8_i,q8_i] = gsw_pchip_interp_SA_CT(v8_tmp,q8_tmp,independent_variable,independent_variable_obs_plus_interp);
    
    CT_i_1 = v1_i;
    tracer_i_1 = rec_factor.*q1_i;

    CT_i_2 = -q2_i.*sin_pi_on_16 + v2_i.*cos_pi_on_16;
    tracer_i_2 = rec_factor.*(q2_i.*cos_pi_on_16 + v2_i.*sin_pi_on_16);
    
    CT_i_3 = -q3_i.*sin_pi_on_8 + v3_i.*cos_pi_on_8;
    tracer_i_3 = rec_factor.*(q3_i.*cos_pi_on_8 + v3_i.*sin_pi_on_8);
    
    CT_i_4 = -q4_i.*sin_3pi_on_16 + v4_i.*cos_3pi_on_16;
    tracer_i_4 = rec_factor.*(q4_i.*cos_3pi_on_16 + v4_i.*sin_3pi_on_16);
    
    CT_i_5 = -q5_i.*sin_pi_on_4 + v5_i.*cos_pi_on_4;
    tracer_i_5 = rec_factor.*(q5_i.*cos_pi_on_4 + v5_i.*sin_pi_on_4);
    
    CT_i_6 = -q6_i.*sin_5pi_on_16 + v6_i.*cos_5pi_on_16;
    tracer_i_6 = rec_factor.*(q6_i.*cos_5pi_on_16 + v6_i.*sin_5pi_on_16);
    
    CT_i_7 = -q7_i.*sin_3pi_on_8 + v7_i.*cos_3pi_on_8;
    tracer_i_7 = rec_factor.*(q7_i.*cos_3pi_on_8 + v7_i.*sin_3pi_on_8);
    
    CT_i_8 = -q8_i.*sin_7pi_on_16 + v8_i.*cos_7pi_on_16;
    tracer_i_8 = rec_factor.*(q8_i.*cos_7pi_on_16 + v8_i.*sin_7pi_on_16);
    
    CT_i_obs_plus_interp = 0.125.*(CT_i_1 + CT_i_2 + CT_i_3 + CT_i_4 + CT_i_5 + CT_i_6 + CT_i_7 + CT_i_8);
    tracer_i_obs_plus_interp = 0.125.*(tracer_i_1 + tracer_i_2 + tracer_i_3 + tracer_i_4 + tracer_i_5 + tracer_i_6 + tracer_i_7 + tracer_i_8);
              
    [min_p_obs, Imin_p_obs] = min(p_obs);
    if min_p_obs ~= 0        
        [Isurface] = find(p_i_tmp < min_p_obs);        
        tracer_i_tooutput = NaN(length(Isurf_and_obs_plus_interp),1);
        CT_i_tooutput = tracer_i_tooutput;

        tracer_i_tooutput(Isurface) = tracer_i_obs_plus_interp(I3(Imin_p_obs));
        CT_i_tooutput(Isurface) = CT_i_obs_plus_interp(I3(Imin_p_obs));
        
        tracer_i_tooutput(Iobs_plus_interp) = tracer_i_obs_plus_interp;
        CT_i_tooutput(Iobs_plus_interp) = CT_i_obs_plus_interp;
    else
        tracer_i_tooutput = tracer_i_obs_plus_interp;
        CT_i_tooutput = CT_i_obs_plus_interp;
    end
    
    tracer_i(Iout,Iprofile) = tracer_i_tooutput(I1);
    CT_i(Iout,Iprofile) = CT_i_tooutput(I1);

end

end
