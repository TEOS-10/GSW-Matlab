function entropy = gsw_entropy_from_pt(SA,pt)

% gsw_entropy_from_pt                     specific entropy of seawater as a 
%                                         function of potential temperature 
%==========================================================================
%
% USAGE:
%  entropy  =  gsw_entropy_from_pt(SA,pt)
%
% DESCRIPTION:
%  Calculates specific entropy of seawater. 
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  pt  =  potential temperature (ITS-90)                          [ deg C ]
%
%  SA & pt need to have the same dimensions.
%
% OUTPUT:
%  entropy  =  specific entropy                                [ J/(kg*K) ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%      
% VERSION NUMBER: 3.03 (29th April, 2013)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.10 of this TEOS-10 Manual. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_entropy_from_pt:  Requires two inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(pt);

if (mt ~= ms | nt ~= ns)
    error('gsw_entropy_from_pt: SA and pt must have same dimensions')
end

if ms == 1
    SA = SA.';
    pt = pt.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

n0 = 0; 
n1 = 1;

pr0 = zeros(size(SA)); 
entropy = -gsw_gibbs(n0,n1,n0,SA,pt,pr0);

if transposed
    entropy = entropy.';
end

end