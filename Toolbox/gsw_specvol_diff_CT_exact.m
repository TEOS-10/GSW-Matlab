function specvol_diff_CT_exact = gsw_specvol_diff_CT_exact(SA,CT,p_shallow,p_deep)

% specvol_diff_CT_exact    specific volume difference between two pressures
%==========================================================================
% 
% USAGE:  
%  specvol_diff_CT_exact = gsw_specvol_diff_CT_exact(SA,CT,p_shallow,p_deep)
%
% DESCRIPTION:
%  Calculates the difference of the specific volume of seawater between 
%  two different pressures, p_deep (the deeper pressure) and p_shallow
%  (the shallower pressure), at the same values of SA and CT. 
%
%  Note that this function uses the full Gibbs function.  There is an 
%  alternative to calling this function, namely gsw_specvol_diff(SA,CT,p),
%  which uses the computationally efficient 75-term expression for specific
%  volume in terms of SA, CT and p (Roquet et al., 2015).  
%
% INPUT:
%  SA         =  Absolute Salinity                                 [ g/kg ]
%  CT         =  Conservative Temperature (ITS-90)                [ deg C ]
%  p_shallow  =  upper sea pressure                                [ dbar ]
%                ( i.e. shallower absolute pressure - 10.1325 dbar ) 
%  p_deep     =  lower sea pressure                                [ dbar ]
%                ( i.e. deeper absolute pressure - 10.1325 dbar )
%
%  SA & CT  need to have the same dimensions.
%  p_shallow and p_deep may have dimensions Mx1 or 1xN or MxN, 
%  where SA and CT are MxN.
%
% OUTPUT:
%  specvol_diff_CT_exact  =  difference of specific volume       [ m^3/kg ]
%                             (deep minus shallow)                      
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.13 (1st June, 2021)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.7.2) of this TEOS-10 Manual. 
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4)
   error('gsw_specvol_diff_CT_exact:  Requires four inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT); 
[mpu,npu] = size(p_shallow);
[mpl,npl] = size(p_deep);

if (mt ~= ms | nt ~= ns)
    error('gsw_specvol_diff_CT_exact: SA and CT must have same dimensions')
end

if (mpu == 1) & (npu == 1)                     
    p_shallow = p_shallow*ones(size(SA));
elseif (ns == npu) & (mpu == 1)              
    p_shallow = p_shallow(ones(1,ms), :);         
elseif (ms == mpu) & (npu == 1)              
    p_shallow = p_shallow(:,ones(1,ns));        
elseif (ns == mpu) & (npu == 1)          
    p_shallow = p_shallow.';                           
    p_shallow = p_shallow(ones(1,ms), :);               
elseif (ms == mpu) & (ns == npu)
    % ok
end

if (mpl == 1) & (npl == 1)                      
    p_deep = p_deep*ones(size(SA));
elseif (ns == npl) & (mpl == 1)                    
    p_deep = p_deep(ones(1,ms), :);               
elseif (ms == mpl) & (npl == 1)                  
    p_deep = p_deep(:,ones(1,ns));                 
elseif (ns == mpl) & (npl == 1)         
    p_deep = p_deep.';                              
    p_deep = p_deep(ones(1,ms), :);            
elseif (ms == mpl) & (ns == npl)
    % ok
else
    error('gsw_specvol_diff_CT_exact: Inputs array dimensions arguments do not agree')
end

if ms == 1
    SA = SA.';
    CT = CT.';
    p_shallow = p_shallow.';
    p_deep = p_deep.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

t_shallow = gsw_t_from_CT(SA,CT,p_shallow);
t_deep = gsw_t_from_CT(SA,CT,p_deep);

specvol_diff_CT_exact = gsw_specvol_t_exact(SA,t_deep,p_deep) ...
    - gsw_specvol_t_exact(SA,t_shallow,p_shallow);

if transposed
    specvol_diff_CT_exact = specvol_CT_diff_exact.';
end

end
