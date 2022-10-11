function [SA_i, CT_i] = gsw_rr68_interp_SA_CT(SA,CT,p,p_i)

% gsw_rr68_interp_SA_CT              Reiniger and Ross (1968) interpolation
%                                                          to p_i on a cast
%==========================================================================
%
% USAGE:
%  [SA_i,CT_i] = gsw_rr68_interp_SA_CT(SA,CT,p,p_i)
%
% DESCRIPTION:
%  Interpolate Absolute Salinity and Conservative Temperature values to
%  arbitrary pressures using the Reiniger and Ross (1968) interpolation
%  scheme.
%  Note that this interpolation scheme requires at least four observed
%  bottles on the cast.
%  Any interpolated bottles that have pressures shallower than the 
%  shallowest observed bottle are set equal to the shallowest observed 
%  bottle.
%
% INPUT: 
%  SA   =  Absolute Salinity                                  [ g/kg ]
%  CT   =  Conservative Temperature (ITS-90)                 [ deg C ]
%  p    =  sea pressure                                       [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  p_i  =  pressures to interpolate to.
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  SA_i = interpolated SA values at pressures p_i.
%  CT_i = interpolated CT values at pressures p_i.
%
% AUTHOR:
%  Paul Barker and Trevor McDougall             [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.12 (25th May, 2020)
%
% References
%  Reiniger, R.F. and C.K. Ross, 1968: A method of interpolation with
%   application to oceanographic data. Deep-Sea Res., 15, 185-193.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4)
    error('gsw_rr68_interp_SA_CT:  Requires four inputs')
end

[Ibad] = find(isnan(SA) | isnan(CT) | isnan(p) | p<-5 | p>12000);
if ~isempty(Ibad)
    SA(Ibad) = [];
    CT(Ibad) = [];
    p(Ibad) = [];
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);
[mp_i,np_i] = size(p_i);

if (ms~=mt) | (ns~=nt) | (ms~=mp) | (ns~=np)
    error('gsw_rr68_interp_SA_CT: SA, CT and p need to have the same dimensions')
end

SA_i = NaN(mp_i,np_i);
CT_i = SA_i;

if isempty(SA)
    warning('gsw_rr68_interp_SA_CT: There is no data to interpolate')
    return
end

if mp < 4 % need at least four bottles to perform this interpolation
    return
end

Ishallow = 1:(mp-1);
Ideep = 2:mp;
d_p = (p(Ideep) - p(Ishallow));

if any(d_p <= 0)
    disp([SA,CT,p])
    warning('gsw_rr68_interp_SA_CT: pressure must be monotonic')
    return
end

SA = SA(:);
CT = CT(:);
p = p(:);
p_i = p_i(:);

[mp,np] = size(p);
[mp_i,np_i] = size(p_i);

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% m = the power variable
m = 1.7;

Ip = (1:mp);
Ip = Ip(:);
Ip_i = (1:mp_i);
Ip_i = Ip_i(:);

[p_at_p_i,Ip_at_p_i,Ip_i_at_p] = intersect(p,p_i);
Ip_ishallow = Ip_i(p_i >= p(1) & p_i <= p(2)); 
Ip_icentral = Ip_i(p_i >= p(2) & p_i <= p(mp-1));
Ip_ideep = Ip_i(p_i >= p(mp-1) & p_i <= p(mp));

if ~isempty(Ip_icentral) 
    % Calculate the 2 outer extrapolated values and the inner interpolated value:
    % eqn (3d)
    Ip_central = interp1q(p,Ip,p_i(Ip_icentral));
    
    Ip_2 = floor(Ip_central);
    Ip_1 = Ip_2 - 1;

    SA_12 = SA(Ip_1) + ((SA(Ip_2) - SA(Ip_1)).*(p_i(Ip_icentral) - p(Ip_1))./(p(Ip_2) - p(Ip_1)));
    CT_12 = CT(Ip_1) + ((CT(Ip_2) - CT(Ip_1)).*(p_i(Ip_icentral) - p(Ip_1))./(p(Ip_2) - p(Ip_1)));
        
    Ip_3 = ceil(Ip_central);
    Ip_4 = Ip_3 + 1;
    SA_34 = SA(Ip_3) + ((SA(Ip_4) - SA(Ip_3)).*(p_i(Ip_icentral) - p(Ip_3))./(p(Ip_4) - p(Ip_3)));
    CT_34 = CT(Ip_3) + ((CT(Ip_4) - CT(Ip_3)).*(p_i(Ip_icentral) - p(Ip_3))./(p(Ip_4) - p(Ip_3)));
        
    %SA_23 = SA(Ip_2) + ((SA(Ip_3) - SA(Ip_2)).*(p_i(Ip_icentral) - p(Ip_2))./(p(Ip_3) - p(Ip_2))));
    %CT_23 = CT(Ip_2) + ((CT(Ip_3) - CT(Ip_2)).*(p_i(Ip_icentral) - p(Ip_2))./(p(Ip_3) - p(Ip_2))));
    % The follwing line produces the same result as the two commented lines above.
    [SA_23, CT_23] = gsw_linear_interp_SA_CT(SA,CT,p,p_i(Ip_icentral));
    SA_23 = SA_23(:);
    CT_23 = CT_23(:);
    
    % Construct the Reiniger & Ross reference curve equation.
    % eqn (3a)
    SAref_num = (abs(SA_23 - SA_34).^m).*SA_12 + (abs(SA_12 - SA_23).^m).*SA_34 ;
    CTref_num = (abs(CT_23 - CT_34).^m).*CT_12 + (abs(CT_12 - CT_23).^m).*CT_34;
    
    SAref_denom = abs(SA_23 - SA_34).^m + abs(SA_12 - SA_23).^m ;
    CTref_denom = abs(CT_23 - CT_34).^m + abs(CT_12 - CT_23).^m ;
    
    [Idenom0] = find(SAref_denom == 0);
    if ~isempty(Idenom0)
        SA_23(Idenom0) = SA_23(Idenom0) + 1e-6;
        SAref_num(Idenom0) = (abs(SA_23(Idenom0) - SA_34(Idenom0)).^m).*SA_12(Idenom0) + (abs(SA_12(Idenom0) - SA_23(Idenom0)).^m).*SA_34(Idenom0);
        SAref_denom(Idenom0) = abs(SA_23(Idenom0) - SA_34(Idenom0)).^m + abs(SA_12(Idenom0) - SA_23(Idenom0)).^m ;
    end
    [Idenom0] = find(CTref_denom == 0);
    if ~isempty(Idenom0)
        CT_23(Idenom0) = CT_23(Idenom0) + 1e-6;
        CTref_num(Idenom0) = (abs(CT_23(Idenom0) - CT_34(Idenom0)).^m).*CT_12(Idenom0) + (abs(CT_12(Idenom0) - CT_23(Idenom0)).^m).*CT_34(Idenom0);
        CTref_denom(Idenom0) = abs(CT_23(Idenom0) - CT_34(Idenom0)).^m + abs(CT_12(Idenom0) - CT_23(Idenom0)).^m ;
    end
    
    SA_ref = 0.5*(SA_23 + (SAref_num./SAref_denom));
    CT_ref = 0.5*(CT_23 + (CTref_num./CTref_denom));
    
    %eqn (3c)
    gamma1_23 = ((p_i(Ip_icentral) - p(Ip_2)).*(p_i(Ip_icentral) - p(Ip_3)))./ ...
        ((p(Ip_1) - p(Ip_2)).*(p(Ip_1) - p(Ip_3)));
    gamma2_31 = ((p_i(Ip_icentral) - p(Ip_3)).*(p_i(Ip_icentral) - p(Ip_1)))./ ...
        ((p(Ip_2) - p(Ip_3)).*(p(Ip_2) - p(Ip_1)));
    gamma3_12 = ((p_i(Ip_icentral) - p(Ip_1)).*(p_i(Ip_icentral) - p(Ip_2)))./ ...
        ((p(Ip_3) - p(Ip_1)).*(p(Ip_3) - p(Ip_2)));
    
    gamma2_34 = ((p_i(Ip_icentral) - p(Ip_3)).*(p_i(Ip_icentral) - p(Ip_4)))./ ...
        ((p(Ip_2) - p(Ip_3)).*(p(Ip_2) - p(Ip_4)));
    gamma3_42 = ((p_i(Ip_icentral) - p(Ip_4)).*(p_i(Ip_icentral) - p(Ip_2)))./ ...
        ((p(Ip_3) - p(Ip_4)).*(p(Ip_3) - p(Ip_2)));
    gamma4_23 = ((p_i(Ip_icentral) - p(Ip_2)).*(p_i(Ip_icentral) - p(Ip_3)))./ ...
        ((p(Ip_4) - p(Ip_2)).*(p(Ip_4) - p(Ip_3)));
    
    %eqn (3b)
    SAP1 = gamma1_23.*SA(Ip_1) + gamma2_31.*SA(Ip_2) + gamma3_12.*SA(Ip_3);
    CTP1 = gamma1_23.*CT(Ip_1) + gamma2_31.*CT(Ip_2) + gamma3_12.*CT(Ip_3);
    
    SAP2 = gamma2_34.*SA(Ip_2) + gamma3_42.*SA(Ip_3) + gamma4_23.*SA(Ip_4);
    CTP2 = gamma2_34.*CT(Ip_2) + gamma3_42.*CT(Ip_3) + gamma4_23.*CT(Ip_4);
    
    %eqn (3)
    SA_ref_minus_SAP1 = SA_ref - SAP1;
    SA_ref_minus_SAP2 = SA_ref - SAP2;
    Inodiff = find(SA_ref_minus_SAP1 == 0 & SA_ref_minus_SAP2 == 0);
    if ~isempty(Inodiff)
        SA_ref(Inodiff) = SA_ref(Inodiff) + 1e-6;
        SA_ref_minus_SAP1(Inodiff) = SA_ref(Inodiff) - SAP1(Inodiff);
        SA_ref_minus_SAP2(Inodiff) = SA_ref(Inodiff) - SAP2(Inodiff);
    end
        
    CT_ref_minus_CTP1 = CT_ref - CTP1;
    CT_ref_minus_CTP2 = CT_ref - CTP2;
    Inodiff = find(CT_ref_minus_CTP1 == 0 & CT_ref_minus_CTP2 == 0);
    if ~isempty(Inodiff)
        CT_ref(Inodiff) = CT_ref(Inodiff) + 1e-6;
        CT_ref_minus_CTP1(Inodiff) = CT_ref(Inodiff) - CTP1(Inodiff);
        CT_ref_minus_CTP2(Inodiff) = CT_ref(Inodiff) - CTP2(Inodiff);
    end
    
    SA_i(Ip_icentral) = (abs(SA_ref_minus_SAP1).*SAP2 + abs(SA_ref_minus_SAP2).*SAP1) ./ ...
        (abs(SA_ref_minus_SAP1) + abs(SA_ref_minus_SAP2));
    CT_i(Ip_icentral) = (abs(CT_ref_minus_CTP1).*CTP2 + abs(CT_ref_minus_CTP2).*CTP1) ./ ...
        (abs(CT_ref_minus_CTP1) + abs(CT_ref_minus_CTP2));
    
    % in some cases the spacing is such that SA_ref_minus_SAP1 and
    % SA_ref_minus_SAP2 are both zero, in these cases we will assign the
    % reference cast bottles to the interpolated cast bottels.
    Inan= find(isnan(SA_i(Ip_icentral)));
    SA_i(Ip_icentral(Inan)) = SA_ref(Inan);
    Inan= find(isnan(CT_i(Ip_icentral)));
    CT_i(Ip_icentral(Inan)) = CT_ref(Inan);
end

if ~isempty(Ip_ishallow)
    %shallow
    Ip_shallow = interp1q(p,Ip,p_i(Ip_ishallow));
    
    Ip_1 = floor(Ip_shallow);
    Ip_2 = ceil(Ip_shallow);
    %eqn (3d)
    %SA_12 = SA(Ip_1) + ((SA(Ip_2) - SA(Ip_1)).*(p_i(Ip_ishallow) - p(Ip_1))./(p(Ip_2) - p(Ip_1))));
    %CT_12 = CT(Ip_1) + ((CT(Ip_2) - CT(Ip_1)).*(p_i(Ip_ishallow) - p(Ip_1))./(p(Ip_2) - p(Ip_1))));
    % The follwing line produces the same result as the two commented lines above.
    [SA_12, CT_12] = gsw_linear_interp_SA_CT(SA,CT,p,p_i(Ip_ishallow));
    SA_12 = SA_12(:);
    CT_12 = CT_12(:);
    
    Ip_3 = Ip_2 + 1;
    Ip_4 = Ip_3 + 1;
    SA_13 = SA(Ip_1) + ((SA(Ip_3) - SA(Ip_1)).*(p_i(Ip_ishallow) - p(Ip_1))./(p(Ip_3) - p(Ip_1)));
    CT_13 = CT(Ip_1) + ((CT(Ip_3) - CT(Ip_1)).*(p_i(Ip_ishallow) - p(Ip_1))./(p(Ip_3) - p(Ip_1)));
    
    SA_23 = SA(Ip_2) + ((SA(Ip_3) - SA(Ip_2)).*(p_i(Ip_ishallow) - p(Ip_2))./(p(Ip_3) - p(Ip_2)));
    CT_23 = CT(Ip_2) + ((CT(Ip_3) - CT(Ip_2)).*(p_i(Ip_ishallow) - p(Ip_2))./(p(Ip_3) - p(Ip_2)));
    
    SA_34 = SA(Ip_3) + ((SA(Ip_4) - SA(Ip_3)).*(p_i(Ip_ishallow) - p(Ip_3))./(p(Ip_4) - p(Ip_3)));
    CT_34 = CT(Ip_3) + ((CT(Ip_4) - CT(Ip_3)).*(p_i(Ip_ishallow) - p(Ip_3))./(p(Ip_4) - p(Ip_3)));
    
    %eqn (3a')
    SAref_num = (abs(SA_12 - SA_23).^m).*SA_34 + (abs(SA_12 - SA_13).^m).*SA_23 ;
    CTref_num = (abs(CT_12 - CT_23).^m).*CT_34 + (abs(CT_12 - CT_13).^m).*CT_23;
    
    SAref_denom = abs(SA_12 - SA_23).^m + abs(SA_12 - SA_13).^m ;
    CTref_denom = abs(CT_12 - CT_23).^m + abs(CT_12 - CT_13).^m ;
    
    [Idenom0] = find(SAref_denom == 0);
    if ~isempty(Idenom0)
        SA_23(Idenom0) = SA_23(Idenom0) + 1e-6;
        SAref_num(Idenom0) = (abs(SA_12(Idenom0) - SA_23(Idenom0)).^m).*SA_34(Idenom0) + (abs(SA_12(Idenom0) - SA_13(Idenom0)).^m).*SA_23(Idenom0);
        SAref_denom(Idenom0) = abs(SA_12(Idenom0) - SA_23(Idenom0)).^m + abs(SA_12(Idenom0) - SA_13(Idenom0)).^m ;
    end
    [Idenom0] = find(CTref_denom == 0);
    if ~isempty(Idenom0)
        CT_23(Idenom0) = CT_23(Idenom0) + 1e-6;
        CTref_num(Idenom0) = (abs(CT_12(Idenom0) - CT_23(Idenom0)).^m).*CT_34(Idenom0) + (abs(CT_12(Idenom0) - CT_13(Idenom0)).^m).*CT_23(Idenom0);
        CTref_denom(Idenom0) = abs(CT_12(Idenom0) - CT_23(Idenom0)).^m + abs(CT_12(Idenom0) - CT_13(Idenom0)).^m ;
    end
    
    SA_ref = 0.5*(SA_12 + (SAref_num./SAref_denom));
    CT_ref = 0.5*(CT_12 + (CTref_num./CTref_denom));
    
    %eqn (3c)
    gamma1_23 = ((p_i(Ip_ishallow) - p(Ip_2)).*(p_i(Ip_ishallow) - p(Ip_3)))./ ...
        ((p(Ip_1) - p(Ip_2)).*(p(Ip_1) - p(Ip_3)));
    gamma2_31 = ((p_i(Ip_ishallow) - p(Ip_3)).*(p_i(Ip_ishallow) - p(Ip_1)))./ ...
        ((p(Ip_2) - p(Ip_3)).*(p(Ip_2) - p(Ip_1)));
    gamma3_12 = ((p_i(Ip_ishallow) - p(Ip_1)).*(p_i(Ip_ishallow) - p(Ip_2)))./ ...
        ((p(Ip_3) - p(Ip_1)).*(p(Ip_3) - p(Ip_2)));
    
    gamma1_24 = ((p_i(Ip_ishallow) - p(Ip_2)).*(p_i(Ip_ishallow) - p(Ip_4)))./ ...
        ((p(Ip_1) - p(Ip_2)).*(p(Ip_1) - p(Ip_4)));
    gamma2_41 = ((p_i(Ip_ishallow) - p(Ip_4)).*(p_i(Ip_ishallow) - p(Ip_1)))./ ...
        ((p(Ip_2) - p(Ip_4)).*(p(Ip_2) - p(Ip_1)));
    gamma4_12 = ((p_i(Ip_ishallow) - p(Ip_1)).*(p_i(Ip_ishallow) - p(Ip_2)))./ ...
        ((p(Ip_4) - p(Ip_1)).*(p(Ip_4) - p(Ip_2)));
    
    %eqn (3b')
    SAP1 = gamma1_23.*SA(Ip_1) + gamma2_31.*SA(Ip_2) + gamma3_12.*SA(Ip_3);
    CTP1 = gamma1_23.*CT(Ip_1) + gamma2_31.*CT(Ip_2) + gamma3_12.*CT(Ip_3);
    
    SAP2 = gamma1_24.*SA(Ip_1) + gamma2_41.*SA(Ip_2) + gamma4_12.*SA(Ip_4);
    CTP2 = gamma1_24.*CT(Ip_1) + gamma2_41.*CT(Ip_2) + gamma4_12.*CT(Ip_4);
    
    %eqn (3)
    
    SA_ref_minus_SAP1 = SA_ref - SAP1;
    SA_ref_minus_SAP2 = SA_ref - SAP2;
    Inodiff = find(SA_ref_minus_SAP1 == 0 & SA_ref_minus_SAP2 == 0);
    if ~isempty(Inodiff)
        SA_ref(Inodiff) = SA_ref(Inodiff) + 1e-6;
        SA_ref_minus_SAP1(Inodiff) = SA_ref(Inodiff) - SAP1(Inodiff);
        SA_ref_minus_SAP2(Inodiff) = SA_ref(Inodiff) - SAP2(Inodiff);
    end
    
    CT_ref_minus_CTP1 = CT_ref - CTP1;
    CT_ref_minus_CTP2 = CT_ref - CTP2;
    Inodiff = find(CT_ref_minus_CTP1 == 0 & CT_ref_minus_CTP2 == 0);
    if ~isempty(Inodiff)
        CT_ref(Inodiff) = CT_ref(Inodiff) + 1e-6;
        CT_ref_minus_CTP1(Inodiff) = CT_ref(Inodiff) - CTP1(Inodiff);
        CT_ref_minus_CTP2(Inodiff) = CT_ref(Inodiff) - CTP2(Inodiff);
    end
    
    SA_i(Ip_ishallow) = (abs(SA_ref_minus_SAP1).*SAP2 + abs(SA_ref_minus_SAP2).*SAP1) ./ ...
        (abs(SA_ref_minus_SAP1) + abs(SA_ref_minus_SAP2));
    CT_i(Ip_ishallow) = (abs(CT_ref_minus_CTP1).*CTP2 + abs(CT_ref_minus_CTP2).*CTP1) ./ ...
        (abs(CT_ref_minus_CTP1) + abs(CT_ref_minus_CTP2));
    
end

if ~isempty(Ip_ideep)
    %deep
    Ip_deep = interp1q(p,Ip,p_i(Ip_ideep));
    
    Ip_2 = floor(Ip_deep);
    Ip_1 = ceil(Ip_deep);
    %eqn (3d)
    %SA_12 = SA(Ip_1) + ((SA(Ip_2) - SA(Ip_1)).*(p_i(Ip_ideep) - p(Ip_1))./(p(Ip_2) - p(Ip_1))));
    %CT_12 = CT(Ip_1) + ((CT(Ip_2) - CT(Ip_1)).*(p_i(Ip_ideep) - p(Ip_1))./(p(Ip_2) - p(Ip_1))));;
    % The follwing line produces the same result as the two commented lines above.
    [SA_12, CT_12] = gsw_linear_interp_SA_CT(SA,CT,p,p_i(Ip_ideep));
    SA_12 = SA_12(:);
    CT_12 = CT_12(:);
    
    Ip_3 = Ip_2 - 1;
    Ip_4 = Ip_3 - 1;
    SA_13 = SA(Ip_1) + ((SA(Ip_3) - SA(Ip_1)).*(p_i(Ip_ideep) - p(Ip_1))./(p(Ip_3) - p(Ip_1)));
    CT_13 = CT(Ip_1) + ((CT(Ip_3) - CT(Ip_1)).*(p_i(Ip_ideep) - p(Ip_1))./(p(Ip_3) - p(Ip_1)));
    
    SA_23 = SA(Ip_2) + ((SA(Ip_3) - SA(Ip_2)).*(p_i(Ip_ideep) - p(Ip_2))./(p(Ip_3) - p(Ip_2)));
    CT_23 = CT(Ip_2) + ((CT(Ip_3) - CT(Ip_2)).*(p_i(Ip_ideep) - p(Ip_2))./(p(Ip_3) - p(Ip_2)));
    
    SA_34 = SA(Ip_3) + ((SA(Ip_4) - SA(Ip_3)).*(p_i(Ip_ideep) - p(Ip_3))./(p(Ip_4) - p(Ip_3)));
    CT_34 = CT(Ip_3) + ((CT(Ip_4) - CT(Ip_3)).*(p_i(Ip_ideep) - p(Ip_3))./(p(Ip_4) - p(Ip_3)));
    
    %eqn (3a')
    SAref_num = (abs(SA_12 - SA_23).^m).*SA_34 + (abs(SA_12 - SA_13).^m).*SA_23 ;
    CTref_num = (abs(CT_12 - CT_23).^m).*CT_34 + (abs(CT_12 - CT_13).^m).*CT_23;
    
    SAref_denom = abs(SA_12 - SA_23).^m + abs(SA_12 - SA_13).^m ;
    CTref_denom = abs(CT_12 - CT_23).^m + abs(CT_12 - CT_13).^m ;
    
    [Idenom0] = find(SAref_denom == 0);
    if ~isempty(Idenom0)
        SA_23(Idenom0) = SA_23(Idenom0) + 1e-6;
        SAref_num(Idenom0) = (abs(SA_12(Idenom0) - SA_23(Idenom0)).^m).*SA_34(Idenom0) + (abs(SA_12(Idenom0) - SA_13(Idenom0)).^m).*SA_23(Idenom0);
        SAref_denom(Idenom0) = abs(SA_12(Idenom0) - SA_23(Idenom0)).^m + abs(SA_12(Idenom0) - SA_13(Idenom0)).^m ;
    end
    [Idenom0] = find(CTref_denom == 0);
    if ~isempty(Idenom0)
        CT_23(Idenom0) = CT_23(Idenom0) + 1e-6;
        CTref_num(Idenom0) = (abs(CT_12(Idenom0) - CT_23(Idenom0)).^m).*CT_34(Idenom0) + (abs(CT_12(Idenom0) - CT_13(Idenom0)).^m).*CT_23(Idenom0);
        CTref_denom(Idenom0) = abs(CT_12(Idenom0) - CT_23(Idenom0)).^m + abs(CT_12(Idenom0) - CT_13(Idenom0)).^m ;
    end
    
    SA_ref = 0.5*(SA_12 + (SAref_num./SAref_denom));
    CT_ref = 0.5*(CT_12 + (CTref_num./CTref_denom));
    
    %eqn (3c)
    gamma1_23 = ((p_i(Ip_ideep) - p(Ip_2)).*(p_i(Ip_ideep) - p(Ip_3)))./ ...
        ((p(Ip_1) - p(Ip_2)).*(p(Ip_1) - p(Ip_3)));
    gamma2_31 = ((p_i(Ip_ideep) - p(Ip_3)).*(p_i(Ip_ideep) - p(Ip_1)))./ ...
        ((p(Ip_2) - p(Ip_3)).*(p(Ip_2) - p(Ip_1)));
    gamma3_12 = ((p_i(Ip_ideep) - p(Ip_1)).*(p_i(Ip_ideep) - p(Ip_2)))./ ...
        ((p(Ip_3) - p(Ip_1)).*(p(Ip_3) - p(Ip_2)));
    
    gamma1_24 = ((p_i(Ip_ideep) - p(Ip_2)).*(p_i(Ip_ideep) - p(Ip_4)))./ ...
        ((p(Ip_1) - p(Ip_2)).*(p(Ip_1) - p(Ip_4)));
    gamma2_41 = ((p_i(Ip_ideep) - p(Ip_4)).*(p_i(Ip_ideep) - p(Ip_1)))./ ...
        ((p(Ip_2) - p(Ip_4)).*(p(Ip_2) - p(Ip_1)));
    gamma4_12 = ((p_i(Ip_ideep) - p(Ip_1)).*(p_i(Ip_ideep) - p(Ip_2)))./ ...
        ((p(Ip_4) - p(Ip_1)).*(p(Ip_4) - p(Ip_2)));
    
    %eqn (3b')
    SAP1 = gamma1_23.*SA(Ip_1) + gamma2_31.*SA(Ip_2) + gamma3_12.*SA(Ip_3);
    CTP1 = gamma1_23.*CT(Ip_1) + gamma2_31.*CT(Ip_2) + gamma3_12.*CT(Ip_3);
    
    SAP2 = gamma1_24.*SA(Ip_1) + gamma2_41.*SA(Ip_2) + gamma4_12.*SA(Ip_4);
    CTP2 = gamma1_24.*CT(Ip_1) + gamma2_41.*CT(Ip_2) + gamma4_12.*CT(Ip_4);
    
    %eqn (3)
    SA_ref_minus_SAP1 = SA_ref - SAP1;
    SA_ref_minus_SAP2 = SA_ref - SAP2;
    Inodiff = find(SA_ref_minus_SAP1 == 0 & SA_ref_minus_SAP2 == 0);
    if ~isempty(Inodiff)
        SA_ref(Inodiff) = SA_ref(Inodiff) + 1e-6;
        SA_ref_minus_SAP1(Inodiff) = SA_ref(Inodiff) - SAP1(Inodiff);
        SA_ref_minus_SAP2(Inodiff) = SA_ref(Inodiff) - SAP2(Inodiff);
    end
    
    CT_ref_minus_CTP1 = CT_ref - CTP1;
    CT_ref_minus_CTP2 = CT_ref - CTP2;
    Inodiff = find(CT_ref_minus_CTP1 == 0 & CT_ref_minus_CTP2 == 0);
    if ~isempty(Inodiff)
        CT_ref(Inodiff) = CT_ref(Inodiff) + 1e-6;
        CT_ref_minus_CTP1(Inodiff) = CT_ref(Inodiff) - CTP1(Inodiff);
        CT_ref_minus_CTP2(Inodiff) = CT_ref(Inodiff) - CTP2(Inodiff);
    end
    
    SA_i(Ip_ideep) = (abs(SA_ref_minus_SAP1).*SAP2 + abs(SA_ref_minus_SAP2).*SAP1) ./ ...
        (abs(SA_ref_minus_SAP1) + abs(SA_ref_minus_SAP2));
    CT_i(Ip_ideep) = (abs(CT_ref_minus_CTP1).*CTP2 + abs(CT_ref_minus_CTP2).*CTP1) ./ ...
        (abs(CT_ref_minus_CTP1) + abs(CT_ref_minus_CTP2));

end    

% Insert any observed bottles that are at the required interpolated
% pressures
SA_i(Ip_i_at_p) = SA(Ip_at_p_i);
CT_i(Ip_i_at_p) = CT(Ip_at_p_i);

% set interpolated bottles that are shallower than the shallowest observed bottle
% to the shallowest observed bottle.
if min(p) < 30
    [Itooshallow] = find(p_i < min(p));
    if ~isempty(Itooshallow)
        SA_i(Itooshallow) = SA(1);
        CT_i(Itooshallow) = CT(1);
    end
end

end

