function [SA_i, CT_i] = gsw_SA_CT_interp(SA,CT,p,p_i)

% gsw_SA_CT_interp                                  SA and CT interpolation
%                                                          to p_i on a cast
%==========================================================================
%
% USAGE:
%  [SA_i, CT_i] = gsw_SA_CT_interp(SA,CT,p,p_i)
%
% DESCRIPTION:
%  Interpolate Absolute Salinity and Conservative Temperature values to
%  arbitrary pressures using the SA-CT diagram.  This programme requires
%  the observed profile to be stable.  Any interpolated bottles that have 
%  pressures shallower than the shallowest observed bottle are set equal to
%  the shallowest observed bottle.
%
%  Note that this interpolation scheme requires at least four observed
%  bottles on the cast.
%
%  If the user has a licence for either the Optimization toolbox or Tomlab
%  CPLEX, then the returned interpolated profile will be stable, However if
%  one of theres lices are not pressent then the interpolated profile may
%  contain unstabilities.
%  
%  This programme is based on the computationally-efficient expression
%  for specific volume in terms of SA, CT and p (Roquet et al., 2015).
%   
%  Note that this 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel". 
%
% INPUT:
%  SA   =  Absolute Salinity                                       [ g/kg ]
%  CT   =  Conservative Temperature (ITS-90)                      [ deg C ]
%  p    =  sea pressure                                            [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  p_i  =  pressures to interpolate to.
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions Mx1 or 1xN or MxN, where SA & CT are MxN.
%  p_i needs to be either a vector or a matrix and have dimensions M_ix1 
%  or M_ixN.
%
% OUTPUT:
%  SA_i = interpolated SA values at pressures p_i.
%  CT_i = interpolated CT values at pressures p_i.
%
% AUTHOR:
%  Paul Barker, Trevor McDougall and Simon Wotherspoon [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.1 (22nd June, 2017)
%
% References
%  Barker, P.M., T.J. McDougall and S.J. Wotherspoon, 2017: An 
%   interpolation method for oceanographic data. JOAT. (To be submitted).
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4)
    error('gsw_SA_CT_interp:  Requires four inputs')
end

[pl,number_of_profiles] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);
[interp_profile_length,np_i] = size(p_i);

if (pl ~= mt) | (number_of_profiles ~= nt)
    error('gsw_SA_CT_interp: SA and CT need to have the same dimensions')
end

if (mp == 1) & (np == 1)
    %error('gsw_SA_CT_interp:  There must be at least 4 bottles')
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
    error('gsw_SA_CT_interp: Inputs array dimensions arguments do not agree')
end 

if interp_profile_length == 1 & np_i > 1
    p_i = p_i.';
    dp_i = diff(p_i);
    if any(dp_i) < 0
        warning('gsw_SA_CT_interp: interpolating pressure must be monotonic')
        return
    end
    [interp_profile_length,np_i] = size(p_i);
elseif interp_profile_length == number_of_profiles & np_i~= number_of_profiles & all(diff(p_i,1,2)) >= 0
    p_i = p_i.';
    [interp_profile_length,np_i] = size(p_i);
elseif any(diff(p_i,1,1)) < 0
    warning('gsw_SA_CT_interp: interpolating pressure must be monotonic')
    return
else
    % Data shape and interval are ok.
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------
    
SA_i = NaN(interp_profile_length, number_of_profiles);
CT_i = SA_i;

for Iprofile = 1:number_of_profiles
    
    data_bottles = SA(:,Iprofile) + CT(:,Iprofile) + p(:,Iprofile);
    [Inn] = find(~isnan(data_bottles));
    if length(Inn) < 2
        SA_i(:,Iprofile) = NaN;
        CT_i(:,Iprofile) = NaN;
        continue
    end
            
    SA_tmp = SA(Inn,Iprofile);
    CT_tmp = CT(Inn,Iprofile);
    p_tmp = p(Inn,Iprofile);
    if np_i > 1
        p_i_tmp = p_i(:,Iprofile);
    else
        p_i_tmp = p_i(:);
    end
    
    CTf_tmp = gsw_CT_freezing_poly(SA_tmp,p_tmp);
    if any(CT_tmp < (CTf_tmp - 0.1))
        [Ifrozen] = find(CT_tmp < (CTf_tmp - 0.1));
        CT_tmp(Ifrozen) = CTf_tmp(Ifrozen);
    end
    
    pl = length(p_tmp);
    Ishallow = 1:(pl-1);
    Ideep = 2:pl;
    dp_tmp = p_tmp(Ideep) - p_tmp(Ishallow);
    if any(dp_tmp <= 0)
        warning('gsw_SA_CT_interp: pressure must be monotonic')
        [p_sort,Ipsort] = (sort(p_tmp));
        SA_sort = SA_tmp(Ipsort);
        CT_sort = CT_tmp(Ipsort);
        [p_tmp,Ipunique] = unique(p_sort);
        SA_tmp = SA_sort(Ipunique);
        CT_tmp = CT_sort(Ipunique);
        pl = length(p_tmp);
        Ishallow = 1:(pl-1);
        Ideep = 2:pl;
        dp_tmp = p_tmp(Ideep) - p_tmp(Ishallow);       
        %return
    end
        
    [N2, pmid] = gsw_Nsquared(SA_tmp,CT_tmp,p_tmp);

    if any(N2 <= 0)
        [matlab_version, matlab_release_date] = version();
        if exist('tomlabVersion') == 2 | ...
                (license('checkout', 'Optimization_Toolbox') == 1 & datenum(matlab_release_date) < 736574) | ...
                exist('cplexqp.p') == 6  
            if license('checkout', 'Optimization_Toolbox') == 1
                warning off  
            end
            [SA_tmp, CT_tmp] = gsw_stabilise_SA_CT(SA_tmp,CT_tmp,p_tmp);
            [N2, pmid] = gsw_Nsquared(SA_tmp,CT_tmp,p_tmp);
        end
    end

    if all(N2 > 0)
        
        intergral_N2 = zeros(pl,1);
        intergral_N2(2:end,1) = cumsum(N2.*dp_tmp);
        p_all = unique(sort([p_tmp; p_i_tmp]));  
        [dummy, Iout, I1] = intersect(p_i_tmp,p_all);
        [dummy, I2, I3] = intersect(p_tmp,p_all);
        
        intergral_N2_all = interp1q(p_tmp,intergral_N2,p_all);
        
        [Itointerp] = find(p_all >= min(p_tmp) & p_all <= max(p_tmp));
        
        SA_i_all = nan(size(p_all));
        CT_i_all = SA_i_all;
        SA_i_linear = SA_i_all;       
        CT_i_linear = SA_i_all;

        [SA_i_all(Itointerp), CT_i_all(Itointerp)] = gsw_pchip_interp_SA_CT(SA_tmp,CT_tmp,intergral_N2,intergral_N2_all(Itointerp));
        
        [SA_i_linear(Itointerp), CT_i_linear(Itointerp)] = gsw_spline_interp_SA_CT(SA_tmp,CT_tmp,p_tmp,p_all(Itointerp),0.8,10000);

        v_i_all = gsw_specvol(SA_i_all, CT_i_all, p_all);
        v_i_linear = gsw_specvol(SA_i_linear, CT_i_linear, p_all);
        
        [Ireplacenan] = find(isnan(v_i_all) & ~isnan(v_i_linear));
        SA_i_all(Ireplacenan) = SA_i_linear(Ireplacenan);
        CT_i_all(Ireplacenan) = CT_i_linear(Ireplacenan);
        v_i_all(Ireplacenan) = v_i_linear(Ireplacenan);
        
        SA_mid = 0.5*(SA_tmp(1:end-1) + SA_tmp(2:end));
        CT_mid = 0.5*(CT_tmp(1:end-1) + CT_tmp(2:end));
        p_mid = 0.5*(p_tmp(1:end-1) + p_tmp(2:end));
        
        v_shallower = gsw_specvol(SA_tmp(1:end-1),CT_tmp(1:end-1),p_mid);
        v_deeper = gsw_specvol(SA_tmp(2:end),CT_tmp(2:end),p_mid);
        v_mid = gsw_specvol(SA_mid,CT_mid,p_mid);
        delta_v_local = -v_mid + 0.5*(v_shallower + v_deeper);
        v_error = 2.*delta_v_local + 1e-7;
        v_error_all = gsw_linear_interp(v_error,p_mid,p_all);
                              
        [max_v_tmp, Imax_v_tmp]  = max(v_i_all(I3));
        max_v_data = max_v_tmp + abs(v_error_all(I3(Imax_v_tmp)));
        if any(max(v_i_all) > max_v_data)
            toolight = 1;
            while toolight == 1
                [Itoolight] = find(v_i_all > max_v_data);
                [Ishallower] = find((I3 - Itoolight(1)) <= 0);
                Iabove = I2(Ishallower(end));
                Iabove_i = I3(Ishallower(end));
                if (Iabove+1) > I3(end)
                    Ibelow_i = I3(end);
                else
                    Ibelow_i = I3(Iabove + 1);
                end
                SA_i_all(Iabove_i:Ibelow_i) = SA_i_linear(Iabove_i:Ibelow_i);
                CT_i_all(Iabove_i:Ibelow_i) = CT_i_linear(Iabove_i:Ibelow_i);
                v_i_all(Iabove_i:Ibelow_i) = v_i_linear(Iabove_i:Ibelow_i);
                if ~any(max(v_i_all) > max_v_data)
                    toolight = 0;
                end
            end
        end
        
        [min_v_tmp, Imin_v_tmp]  = min(v_i_all(I3));
        min_v_data = min_v_tmp - abs(v_error_all(I3(Imin_v_tmp)));
        if any(min(v_i_all) < min_v_data)
            tooheavy = 1;
            while tooheavy == 1
                [Itooheavy] = find(v_i_all < min_v_data);
                [Ishallower] = find((I3 - Itooheavy(1)) <= 0);
                Iabove = I2(Ishallower(end));
                Iabove_i = I3(Ishallower(end));
                if (Iabove+1) > I3(end)
                    Ibelow_i = I3(end);
                else
                    Ibelow_i = I3(Iabove + 1);
                end
                SA_i_all(Iabove_i:Ibelow_i) = SA_i_linear(Iabove_i:Ibelow_i);
                CT_i_all(Iabove_i:Ibelow_i) = CT_i_linear(Iabove_i:Ibelow_i);
                v_i_all(Iabove_i:Ibelow_i) = v_i_linear(Iabove_i:Ibelow_i);
                if ~any(min(v_i_all) < min_v_data)
                    tooheavy = 0;
                end
            end
        end
        
        CTf_i_all = gsw_CT_freezing_poly(SA_i_all,p_all);
        if any(CT_i_linear < (CTf_i_all - 0.1))
            [ICTf_i_obs] = find(CT_i_linear < (CTf_i_all - 0.1));
            CTf_i_all(ICTf_i_obs) = CT_i_linear(ICTf_i_obs);
        end
        if any(CT_i_all < (CTf_i_all - 0.1))
            frozen = 1;
            while frozen == 1
                [Ifrozen] = find(CT_i_all < (CTf_i_all - 0.1));
                [Ishallower] = find((I3 - Ifrozen(1)) <= 0);
                Iabove = I2(Ishallower(end));
                Iabove_i = I3(Ishallower(end));
                if (Iabove+1) > I3(end)
                    Ibelow_i = I3(end);
                else
                    Ibelow_i = I3(Iabove + 1);
                end
                SA_i_all(Iabove_i:Ibelow_i) = SA_i_linear(Iabove_i:Ibelow_i);
                CT_i_all(Iabove_i:Ibelow_i) = CT_i_linear(Iabove_i:Ibelow_i);
                CTf_i_all(Iabove_i:Ibelow_i) = gsw_CT_freezing_poly(SA_i_all(Iabove_i:Ibelow_i),p_all(Iabove_i:Ibelow_i));
                if any(CT_i_linear < (CTf_i_all - 0.1))
                    [ICTf_i_obs] = find(CT_i_linear < (CTf_i_all - 0.1));
                    CTf_i_all(ICTf_i_obs) = CT_i_linear(ICTf_i_obs);
                end
                if ~any(CT_i_all < (CTf_i_all - 0.1))
                    frozen = 0;
                end
            end
        end
        
        [min_p_tmp, Iminp] = min(p_tmp);
        if min_p_tmp ~= 0
            SA_i_all(p_i_tmp < min_p_tmp) = SA_i_all(I3(Iminp));
            CT_i_all(p_i_tmp < min_p_tmp) = CT_i_all(I3(Iminp));
        end
                
        if exist('tomlabVersion') == 2 | ...
                (license('checkout', 'Optimization_Toolbox') == 1 & length(I1) < 2000) | ...
                exist('cplexqp.p') == 6  
            [SA_out, CT_out] = gsw_stabilise_SA_CT(SA_i_all(I1),CT_i_all(I1),p_all(I1));
            SA_i(Iout,Iprofile) = SA_out(:);
            CT_i(Iout,Iprofile) = CT_out(:);
        else
            SA_i(Iout,Iprofile) = SA_i_all(I1);
            CT_i(Iout,Iprofile) = CT_i_all(I1);
        end
        
    else
         [SA_i(:,Iprofile), CT_i(:,Iprofile)] = gsw_spline_interp_SA_CT(SA_tmp,CT_tmp,p_tmp,p_i,0.8,10000);
%        [SA_i(:,Iprofile), CT_i(:,Iprofile)] = gsw_linear_interp_SA_CT(SA_tmp,CT_tmp,p_tmp,p_i);
     end
    
end

end


function data_i = gsw_linear_interp(data,p,p_i)

% gsw_linear_interp                   linear interpolation to p_i on a cast
%==========================================================================
% This function interpolates the cast with respect to the interpolating 
% variable p. This function finds the values of the inputed data at p_i on
% this cast.
%
% VERSION NUMBER: 3.05 (6th July 2016)
%
% This function was adapted from Matlab's interp1q.
%==========================================================================

p = p(:);
data = data(:);
p_i = p_i(:);
data_i = NaN(size(p_i));

[min_p, Imin_p] = min(p);

data_i(p_i <= min_p) = data(Imin_p);% Set equal to the shallowest bottle.

[max_p,Imax_p] = max(p);
data_i(p_i >= max_p) = data(Imax_p);% Set equal to the deepest bottle.

xi = p_i(p_i >= min_p & p_i <= max_p);

if ~isempty(xi)
    
    x = p;
    
    siz = size(xi);
    if ~isscalar(xi)
        [xxi, k] = sort(xi);
        [dummy, j] = sort([x;xxi]);
        r(j) = 1:length(j);
        r = r(length(x)+1:end) - (1:length(xxi));
        r(k) = r;
        r(xi==x(end)) = length(x)-1;
        ind = find((r>0) & (r<length(x)));
        ind = ind(:);
        var_ri = NaN(length(xxi),size(data,2),superiorfloat(x,data,xi));
        rind = r(ind);
        xrind = x(rind);
        u = (xi(ind)-xrind)./(x(rind+1)-xrind);
        datarind = data(rind,:);
        if exist('bsxfun','builtin') == 5
            var_ri(ind,:) = datarind + bsxfun(@times,data(rind+1,:)-datarind,u);
        else
            var_ri(ind,:) = datarind + (data(rind+1,:)-datarind).*u;
        end
    else
        % Special scalar xi case
        r = find(x <= xi,1,'last');
        r(xi==x(end)) = length(x)-1;
        if isempty(r) || r<=0 || r>=length(x)
            var_ri = NaN(1,size(data,2),superiorfloat(x,data,xi));
        else
            u = (xi-x(r))./(x(r+1)-x(r));
            varr = data(r,:);
            if exist('bsxfun','builtin') == 5
                var_ri = varr + bsxfun(@times,data(r+1,:)-varr,u);
            else
                var_ri = varr + (data(r+1,:)-varr).*u;
            end
        end
    end
    
    if min(size(var_ri)) == 1 && numel(xi) > 1
        var_ri = reshape(var_ri,siz);
    end
    
    data_i(p_i >= min_p & p_i <= max_p) = var_ri;
end
end

