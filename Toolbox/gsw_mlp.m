function mlp = gsw_mlp(SA,CT,p)

% gsw_mlp                                              mixed-layer pressure
%==========================================================================
% 
% USAGE:
%  mlp = gsw_mlp(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the mixed-layer pressure as described in de Boyer Montégut 
%  et al. (2004).  The mlp is always deeper than 20 dbar, if the initial
%  estimate of the mlp is less than 20 dbar, the temperature and salinity  
%  of the bottles in the top 5 dbar are set to that of the bottle closest 
%  to 5 dbar.  This removes the effect if a thin layer of fresh water, 
%  such as that from a river outflow or from rain.
%
%  Note that this 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  mlp  =  mixed-layer pressure                                    [ dbar ]
%  
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05.5 (3rd June, 2016)
%
% REFERENCES:
%  de Boyer Montégut, C., G. Madec, A.S. Fischer, A. Lazar and D. Iudicone 
%   2004: Mixed layer depth over the global ocean: An examination of 
%   profile data and a profile-based climatology, J. Geophys. Res., 109,
%   C12003, doi:10.1029/2004JC002378.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling, 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
    error('gsw_mlp: requires 3 input arguments')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_mlp: SA and CT must have same dimensions')
end

if (mp == 1) & (np == 1)
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)
    p = p(ones(1,ms), :);     
elseif (ms == mp) & (np == 1)  
    p = p(:,ones(1,ns));   
elseif (ns == mp) & (np == 1)  
    p = p.';    
    p = p(ones(1,ms),:); 
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_mlp: Inputs array dimensions arguments do not agree')
end

[ms,ns] = size(SA);
mlp = nan(ns,1);

rho0 = gsw_rho(SA,CT,zeros(ms,ns));

%--------------------------------------------------------------------------
% This function calculates density using the computationally-efficient 
% 75-term expression for specific volume in terms of SA, CT and p.  If one
% wanted to compute specific volume with the full TEOS-10 Gibbs function 
% expression, the following lines of code will enable this.
%
%  t = gsw_t_from_CT(SA,CT,p);
%  rho0 = gsw_rho_t_exact(SA,t,zeros(mp,np));
%
%---------------This is the end of the alternative code--------------------

for Iprofile = 1:ns
    
    [Inn] = find(~isnan(rho0(:,Iprofile)));
    
    if ~isempty(Inn)
        [dummy, I1] = (sort(p(Inn,Iprofile)));
        [dummy2, I2] = unique(dummy);
        
        p_tmp = p(Inn(I1(I2)),Iprofile);
        rho0_tmp = rho0(Inn(I1(I2)),Iprofile);
        
        min_rho0 = min(rho0_tmp);
        
        if isnan(min_rho0)  %  the profile is junk
            mlp(Iprofile) = NaN;
        elseif min(p_tmp) > 20  % the profile starts at pressure greater than 20 dbar
            mlp(Iprofile) = NaN;
        else
            diff_rho0 = (min_rho0 + 0.3) - rho0_tmp;
            I3 = find(diff_rho0 > 0);
            mlp(Iprofile) = p_tmp(I3(end));
        end
        
        dmlp = mlp(Iprofile) - min(p_tmp);
        
        % If the mlp is less than 20 dbar it is possible that this density 
        % difference is being effected by a thin fresh water layer at the 
        % surface, set the salinities and temperatures of the bottles in 
        % the top section of the cast to be equal to that of the bottle 
        % closest to 5 dbar.
        if dmlp < 20
            
            SA_tmp = SA(Inn(I1(I2)),Iprofile);
            CT_tmp = CT(Inn(I1(I2)),Iprofile);
            
            [dummy, I4] = min(abs(p_tmp - 5));
            SA_tmp(1:I4-1) = SA_tmp(I4);
            CT_tmp(1:I4-1) = CT_tmp(I4);
            
            rho0_tmp(1:I4-1) = gsw_rho(SA_tmp(1:I4-1),CT_tmp(1:I4-1),zeros(length(1:I4-1),1));
            
            min_rho0_tmp = min(rho0_tmp);
            
            diff_rho0_tmp = (min_rho0_tmp + 0.3) - rho0_tmp;
            
            I3 = find(diff_rho0_tmp > 0);
            mlp(Iprofile) = p_tmp(I3(end));        
            dmlp = mlp(Iprofile) - min(p_tmp);
            
            if dmlp < 20
                mlp(Iprofile) = NaN;
            end
            
        end
    else
        mlp(Iprofile) = NaN;
    end
end

end