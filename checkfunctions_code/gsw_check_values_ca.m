% gsw_data = 'gsw_data_v3_0.mat';
% gsw_data_file = which(gsw_data);
% load (gsw_data_file,'gsw_cv');

%load library\gsw_chck_vals_file.mat gsw_cv
load gsw_chck_vals_file_v3_06_12.mat gsw_cv

SA_chck_cast = gsw_cv.SA_chck_cast;
t_chck_cast = gsw_cv.t_chck_cast;
p_chck_cast = gsw_cv.p_chck_cast;
lat_chck_cast = gsw_cv.lat_chck_cast;
long_chck_cast = gsw_cv.long_chck_cast;
Rt_chck_cast = gsw_cv.Rt_chck_cast;

SA_Arctic = gsw_cv.SA_Arctic;
t_Arctic = gsw_cv.t_Arctic;
p_Arctic = gsw_cv.p_Arctic;
SA_seaice = gsw_cv.SA_seaice;
t_ice = gsw_cv.t_seaice;
% p_ice = gsw_cv.p_Arctic;
w_ice = gsw_cv.w_seaice;

pr = 0;
pr_05 = 500;

for Iit = 1:4
    if Iit == 2
        SA_chck_cast(:,4:6) = SA_chck_cast(:,1:3) + 1.3e-10;
        t_chck_cast(:,4:6) = t_chck_cast(:,1:3);
        p_chck_cast(:,4:6) = p_chck_cast(:,1:3);
        lat_chck_cast(4:6) = lat_chck_cast(1:3);
        long_chck_cast(4:6) = long_chck_cast(1:3);
        Rt_chck_cast(4:6) = Rt_chck_cast(1:3);  
        SA_Arctic(:,4:6) = SA_Arctic(:,1:3) + 1.3e-10;
        t_Arctic(:,4:6) = t_Arctic(:,1:3);
        p_Arctic(:,4:6) = p_Arctic(:,1:3);
        SA_seaice(:,4:6) = SA_seaice(:,1:3) + 1.3e-10;
        t_ice(:,4:6) = t_ice(:,1:3);
%         p_ice(:,4:6) = p_ice(:,1:3);
        w_ice(:,4:6) = w_ice(:,1:3);
    elseif Iit == 3
        SA_chck_cast(:,7:9) = gsw_cv.SA_chck_cast(:,1:3);
        t_chck_cast(:,7:9) = gsw_cv.t_chck_cast(:,1:3) + 6e-10;
        p_chck_cast(:,7:9) = p_chck_cast(:,1:3);
        lat_chck_cast(7:9) = lat_chck_cast(1:3);
        long_chck_cast(7:9) = long_chck_cast(1:3);
        Rt_chck_cast(7:9) = Rt_chck_cast(1:3);
        SA_Arctic(:,7:9) = SA_Arctic(:,1:3);
        t_Arctic(:,7:9) = t_Arctic(:,1:3) + 6e-10;
        p_Arctic(:,7:9) = p_Arctic(:,1:3);
        SA_seaice(:,7:9) = SA_seaice(:,1:3);
        t_ice(:,7:9) = t_ice(:,1:3) + 6e-10;
%         p_ice(:,7:9) = p_ice(:,1:3);
        w_ice(:,7:9) = w_ice(:,1:3);
    elseif Iit == 4
        SA_chck_cast(:,10:12) = SA_chck_cast(:,1:3);
        t_chck_cast(:,10:12) = t_chck_cast(:,1:3);
        p_chck_cast(:,10:12) = p_chck_cast(:,1:3) + 2.3e-8;
        lat_chck_cast(10:12) = lat_chck_cast(1:3);
        long_chck_cast(10:12) = long_chck_cast(1:3);
        Rt_chck_cast(10:12) = Rt_chck_cast(1:3);
        SA_Arctic(:,10:12) = SA_Arctic(:,1:3);
        t_Arctic(:,10:12) = t_Arctic(:,1:3);
        p_Arctic(:,10:12) = p_Arctic(:,1:3) + 2.3e-8;
        SA_seaice(:,10:12) = SA_seaice(:,1:3);
        t_ice(:,10:12) = t_ice(:,1:3);
%         p_ice(:,10:12) = p_ice(:,1:3) + 2.3e-8;
        w_ice(:,10:12) = w_ice(:,1:3);
    end
end
CT_Arctic = gsw_CT_from_t(SA_Arctic,t_Arctic,p_Arctic);
t_seaice = t_ice;
w_seaice = w_ice;
SA_bulk = (1 - w_ice).*SA_Arctic;
h_bulk = (1 - w_ice).*gsw_enthalpy_CT_exact(SA_Arctic,CT_Arctic,p_Arctic) + w_ice.*gsw_enthalpy_ice(t_ice,p_Arctic);
h_pot_bulk = (1 - w_ice)*gsw_cp0.*gsw_CT_from_t(SA_Arctic,CT_Arctic,p_Arctic) + w_ice.*gsw_enthalpy_ice(t_ice,p_Arctic);

clear gsw_cv.SP_chck_cast

SP_chck_cast = gsw_SP_from_SA(SA_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
p_chck_cast_shallow = p_chck_cast;
p_chck_cast_deep = p_chck_cast + 10;

delta_p_chck_cast = nan(size(SA_chck_cast));
delta_p_chck_cast(1,:) = p_chck_cast(1,:) - 0;
for I = 2:45
    delta_p_chck_cast(I,:) = p_chck_cast(I,:) - p_chck_cast(I-1,:);
end
delta_p_chck_cast(:,10:12) = delta_p_chck_cast(:,10:12) + 2.3e-8;
% delta_p_chck_cast = delta_p_chck_cast;

Neutral_Density = 26.8*ones(1,length(lat_chck_cast));
% Neutral_Density = Neutral_Density;

p_Neutral_Density = [200 550 91 200 550 91 200 550 91 200 550 91];
% gsw_cv.p_Neutral_Density = p_Neutral_Density;

SA_ref = 30;
CT_ref = 5;
p_i = ([0:100:5000]).';

%% Practical Salinity (SP):- PSS-78 

C = gsw_C_from_SP(SP_chck_cast,t_chck_cast,p_chck_cast);
C_from_SP_error = nan(45,9);
for I = 1:9
    C_from_SP_error(:,I) = abs(C(:,I) - C(:,I+3));
end
gsw_cv.C_from_SP_ca = nanmax(nanmax(C_from_SP_error));
  
SP_from_C = gsw_SP_from_C(C,t_chck_cast,p_chck_cast);
SP_from_C_error = nan(45,9);
for I = 1:9
    SP_from_C_error(:,I) = abs(SP_from_C(:,I) - SP_from_C(:,I+3));
end
gsw_cv.SP_from_C_ca = nanmax(nanmax(SP_from_C_error));

R = gsw_R_from_SP(SP_chck_cast,t_chck_cast,p_chck_cast);
R_from_SP_error = nan(45,9);
for I = 1:9
    R_from_SP_error(:,I) = abs(R(:,I) - R(:,I+3));
end
gsw_cv.R_from_SP_ca = nanmax(nanmax(R_from_SP_error));
  
SP_from_R = gsw_SP_from_C(R,t_chck_cast,p_chck_cast);
SP_from_R_error = nan(45,9);
for I = 1:9
    SP_from_R_error(:,I) = abs(SP_from_R(:,I) - SP_from_R(:,I+3));
end
gsw_cv.SP_from_R_ca = nanmax(nanmax(SP_from_R_error));

Rt = gsw_Rt_from_SP(SP_chck_cast,t_chck_cast,p_chck_cast);
SP_salinometer = gsw_SP_salinometer(Rt,t_chck_cast);
SP_salinometer_error = nan(45,9);
for I = 1:9
    SP_salinometer_error(:,I) = abs(SP_salinometer(:,I) - SP_salinometer(:,I+3));
end
gsw_cv.SP_salinometer_ca = nanmax(nanmax(SP_salinometer_error));

SP_from_SK = gsw_SP_from_SK(SP_chck_cast);
SP_from_SK_error = nan(45,9);
for I = 1:9
    SP_from_SK_error(:,I) = abs(SP_from_SK(:,I) - SP_from_SK(:,I+3));
end
gsw_cv.SP_from_SK_ca = nanmax(nanmax(SP_from_SK_error));

%% Absolute Salinity (SA) and Preformed Salinity (Sstar) 
 
SA_from_SP = gsw_SA_from_SP(SP_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
SA_from_SP_error = nan(45,9);
for I = 1:9
    SA_from_SP_error(:,I) = abs(SA_from_SP(:,I) - SA_from_SP(:,I+3));
end
gsw_cv.SA_from_SP_ca = nanmax(nanmax(SA_from_SP_error));
 
Sstar_from_SP = gsw_Sstar_from_SP(SP_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
Sstar_from_SP_error = nan(45,9);
for I = 1:9
    Sstar_from_SP_error(:,I) = abs(Sstar_from_SP(:,I) - Sstar_from_SP(:,I+3));
end
gsw_cv.Sstar_from_SP_ca = nanmax(nanmax(Sstar_from_SP_error));
 
%% Conservative Temperature (CT) 
  
CT_from_t = gsw_CT_from_t(SA_chck_cast,t_chck_cast,p_chck_cast);
CT_from_t_error = nan(45,9);
for I = 1:9
    CT_from_t_error(:,I) = abs(CT_from_t(:,I) - CT_from_t(:,I+3));
end
gsw_cv.CT_from_t_ca = nanmax(nanmax(CT_from_t_error));
CT_chck_cast = gsw_CT_from_t(SA_chck_cast,t_chck_cast,p_chck_cast);

%% other conversions between temperatures, salinities, pressure and height 

deltaSA_from_SP = gsw_deltaSA_from_SP(SP_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
deltaSA_from_SP_error = nan(45,9);
for I = 1:9
    deltaSA_from_SP_error(:,I) = abs(deltaSA_from_SP(:,I) - deltaSA_from_SP(:,I+3));
end
gsw_cv.deltaSA_from_SP_ca = nanmax(nanmax(deltaSA_from_SP_error));

[SA_SA_Sstar_from_SP, Sstar_SA_Sstar_from_SP] = gsw_SA_Sstar_from_SP(SP_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
SA_SA_Sstar_from_SP_error = nan(45,9);
Sstar_SA_Sstar_from_SP_error = nan(45,9);
for I = 1:9
    SA_SA_Sstar_from_SP_error(:,I) = abs(SA_SA_Sstar_from_SP(:,I) - SA_SA_Sstar_from_SP(:,I+3));
end
gsw_cv.SA_SA_Sstar_from_SP_ca = nanmax(nanmax(SA_SA_Sstar_from_SP_error));
 for I = 1:9
    Sstar_SA_Sstar_from_SP_error(:,I) = abs(Sstar_SA_Sstar_from_SP(:,I) - Sstar_SA_Sstar_from_SP(:,I+3));
end
gsw_cv.Sstar_SA_Sstar_from_SP_ca = nanmax(nanmax(Sstar_SA_Sstar_from_SP_error));

SR_from_SP = gsw_SR_from_SP(SP_chck_cast);
SR_from_SP_error = nan(45,9);
for I = 1:9
    SR_from_SP_error(:,I) = abs(SR_from_SP(:,I) - SR_from_SP(:,I+3));
end
gsw_cv.SR_from_SP_ca = nanmax(nanmax(SR_from_SP_error));

SP_from_SR = gsw_SP_from_SR(SR_from_SP);
SP_from_SR_error = nan(45,9);
for I = 1:9
    SP_from_SR_error(:,I) = abs(SP_from_SR(:,I) - SP_from_SR(:,I+3));
end
gsw_cv.SP_from_SR_ca = nanmax(nanmax(SP_from_SR_error));

SP_from_SA = gsw_SP_from_SA(SA_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
SP_from_SA_error = nan(45,9);
for I = 1:9
    SP_from_SA_error(:,I) = abs(SP_from_SA(:,I) - SP_from_SA(:,I+3));
end
gsw_cv.SP_from_SA_ca = nanmax(nanmax(SP_from_SA_error));

Sstar_from_SA = gsw_Sstar_from_SA(SA_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
Sstar_from_SA_error = nan(45,9);
for I = 1:9
    Sstar_from_SA_error(:,I) = abs(Sstar_from_SA(:,I) - Sstar_from_SA(:,I+3));
end
gsw_cv.Sstar_from_SA_ca = nanmax(nanmax(Sstar_from_SA_error));
 
SA_from_Sstar = gsw_SA_from_Sstar(Sstar_from_SP,p_chck_cast,long_chck_cast,lat_chck_cast);
SA_from_Sstar_error = nan(45,9);
for I = 1:9
    SA_from_Sstar_error(:,I) = abs(SA_from_Sstar(:,I) - SA_from_Sstar(:,I+3));
end
gsw_cv.SA_from_Sstar_ca = nanmax(nanmax(SA_from_Sstar_error));
 
SP_from_Sstar = gsw_SP_from_Sstar(Sstar_from_SP,p_chck_cast,long_chck_cast,lat_chck_cast);
SP_from_Sstar_error = nan(45,9);
for I = 1:9
    SP_from_Sstar_error(:,I) = abs(SP_from_Sstar(:,I) - SP_from_Sstar(:,I+3));
end
gsw_cv.SP_from_Sstar_ca = nanmax(nanmax(SP_from_Sstar_error));

t_from_CT = gsw_t_from_CT(SA_chck_cast,CT_chck_cast,p_chck_cast);
t_from_CT_error = nan(45,9);
for I = 1:9
    t_from_CT_error(:,I) = abs(t_from_CT(:,I) - t_from_CT(:,I+3));
end
gsw_cv.t_from_CT_ca = nanmax(nanmax(t_from_CT_error));

pt = gsw_pt_from_CT(SA_chck_cast,CT_chck_cast);
pt_error = nan(45,9);
for I = 1:9
    pt_error(:,I) = abs(pt(:,I) - pt(:,I+3));
end
gsw_cv.pt_from_CT_ca = nanmax(nanmax(pt_error));

CT_from_pt = gsw_CT_from_pt(SA_chck_cast,pt);
CT_from_pt_error = nan(45,9);
for I = 1:9
    CT_from_pt_error(:,I) = abs(CT_from_pt(:,I) - CT_from_pt(:,I+3));
end
gsw_cv.CT_from_pt_ca = nanmax(nanmax(CT_from_pt_error));

pot_enthalpy = gsw_pot_enthalpy_from_pt(SA_chck_cast,pt);
pot_enthalpy_error = nan(45,9);
for I = 1:9
    pot_enthalpy_error(:,I) = abs(pot_enthalpy(:,I) - pot_enthalpy(:,I+3));
end
gsw_cv.pot_enthalpy_from_pt_ca = nanmax(nanmax(pot_enthalpy_error));

pt0 = gsw_pt0_from_t(SA_chck_cast,t_chck_cast,p_chck_cast);
pt0_error = nan(45,9);
for I = 1:9
    pt0_error(:,I) = abs(pt0(:,I) - pt0(:,I+3));
end
gsw_cv.pt0_from_t_ca = nanmax(nanmax(pt0_error));

pt_from_t = gsw_pt_from_t(SA_chck_cast,t_chck_cast,p_chck_cast,pr);
pt_from_t_error = nan(45,9);
for I = 1:9
    pt_from_t_error(:,I) = abs(pt_from_t(:,I) - pt_from_t(:,I+3));
end
gsw_cv.pt_from_t_ca = nanmax(nanmax(pt_from_t_error));
   
t90_from_t68 = gsw_t90_from_t68(t_chck_cast);
t90_from_t68_error = nan(45,9);
for I = 1:9
    t90_from_t68_error(:,I) = abs(t90_from_t68(:,I) - t90_from_t68(:,I+3));
end
gsw_cv.t90_from_t68_ca = nanmax(nanmax(t90_from_t68_error));
  
t90_from_t48 = gsw_t90_from_t48(t_chck_cast);
t90_from_t48_error = nan(45,9);
for I = 1:9
    t90_from_t48_error(:,I) = abs(t90_from_t48(:,I) - t90_from_t48(:,I+3));
end
gsw_cv.t90_from_t48_ca = nanmax(nanmax(t90_from_t48_error));

z_from_p = gsw_z_from_p(p_chck_cast,lat_chck_cast);
z_from_p_error = nan(45,9);
for I = 1:9
    z_from_p_error(:,I) = abs(z_from_p(:,I) - z_from_p(:,I+3));
end
gsw_cv.z_from_p_ca = nanmax(nanmax(z_from_p_error));
  
p_from_z = gsw_p_from_z(z_from_p,lat_chck_cast);
p_from_z_error = nan(45,9);
for I = 1:9
    p_from_z_error(:,I) = abs(p_from_z(:,I) - p_from_z(:,I+3));
end
gsw_cv.p_from_z_ca = nanmax(nanmax(p_from_z_error));
    
depth_from_z = gsw_depth_from_z(z_from_p);
depth_from_z_error = nan(45,9);
for I = 1:9
    depth_from_z_error(:,I) = abs(depth_from_z(:,I) - depth_from_z(:,I+3));
end
gsw_cv.depth_from_z_ca = nanmax(nanmax(depth_from_z_error));

z_from_depth = gsw_z_from_depth(depth_from_z);
z_from_depth_error = nan(45,9);
for I = 1:9
    z_from_depth_error(:,I) = abs(z_from_depth(:,I) - z_from_depth(:,I+3));
end
gsw_cv.z_from_depth_ca = nanmax(nanmax(z_from_depth_error));

Abs_Pressure = gsw_Abs_Pressure_from_p(p_chck_cast);
Abs_Pressure_error = nan(45,8);
for I = 1:8
    Abs_Pressure_error(:,I) = abs(Abs_Pressure(:,I) - Abs_Pressure(:,I+3));
end
gsw_cv.Abs_Pressure_from_p_ca = nanmax(nanmax(Abs_Pressure_error));

p_from_Abs_Pressure = gsw_p_from_Abs_Pressure(Abs_Pressure);
p_from_Abs_Pressure_error = nan(45,8);
for I = 1:8
    p_from_Abs_Pressure_error(:,I) = abs(p_from_Abs_Pressure(:,I) - p_from_Abs_Pressure(:,I+3));
end
gsw_cv.p_from_Abs_Pressure_ca = nanmax(nanmax(p_from_Abs_Pressure_error));

entropy_from_CT = gsw_entropy_from_CT(SA_chck_cast,CT_chck_cast);
entropy_from_CT_error = nan(45,9);
for I = 1:9
    entropy_from_CT_error(:,I) = abs(entropy_from_CT(:,I) - entropy_from_CT(:,I+3));
end
gsw_cv.entropy_from_CT_ca = nanmax(nanmax(entropy_from_CT_error));
 
CT_from_entropy = gsw_CT_from_entropy(SA_chck_cast,entropy_from_CT);
CT_from_entropy_error = nan(45,9);
for I = 1:9
    CT_from_entropy_error(:,I) = abs(CT_from_entropy(:,I) - CT_from_entropy(:,I+3));
end
gsw_cv.CT_from_entropy_ca = nanmax(nanmax(CT_from_entropy_error));

entropy_from_pt = gsw_entropy_from_pt(SA_chck_cast,pt);
entropy_from_pt_error = nan(45,9);
for I = 1:9
    entropy_from_pt_error(:,I) = abs(entropy_from_pt(:,I) - entropy_from_pt(:,I+3));
end
gsw_cv.entropy_from_pt_ca = nanmax(nanmax(entropy_from_pt_error));

pt_from_entropy = gsw_pt_from_entropy(SA_chck_cast,entropy_from_pt);
pt_from_entropy_error = nan(45,9);
for I = 1:9
    pt_from_entropy_error(:,I) = abs(pt_from_entropy(:,I) - pt_from_entropy(:,I+3));
end
gsw_cv.pt_from_entropy_ca = nanmax(nanmax(pt_from_entropy_error));

entropy_from_t = gsw_entropy_from_t(SA_chck_cast,t_chck_cast,p_chck_cast);
entropy_from_t_error = nan(45,9);
for I = 1:9
    entropy_from_t_error(:,I) = abs(entropy_from_t(:,I) - entropy_from_t(:,I+3));
end
gsw_cv.entropy_from_t_ca = nanmax(nanmax(entropy_from_t_error));

t_from_entropy = gsw_t_from_entropy(SA_chck_cast,entropy_from_t,p_chck_cast);
t_from_entropy_error = nan(45,9);
for I = 1:9
    t_from_entropy_error(:,I) = abs(t_from_entropy(:,I) - t_from_entropy(:,I+3));
end
gsw_cv.t_from_entropy_ca = nanmax(nanmax(t_from_entropy_error));

adiabatic_lapse_rate_from_CT = gsw_adiabatic_lapse_rate_from_CT(SA_chck_cast,CT_chck_cast,p_chck_cast);
adiabatic_lapse_rate_from_CT_error = nan(45,9);
for I = 1:9
    adiabatic_lapse_rate_from_CT_error(:,I) = abs(adiabatic_lapse_rate_from_CT(:,I) - adiabatic_lapse_rate_from_CT(:,I+3));
end
gsw_cv.adiabatic_lapse_rate_from_CT_ca = nanmax(nanmax(adiabatic_lapse_rate_from_CT_error));

adiabatic_lapse_rate_from_t = gsw_adiabatic_lapse_rate_from_t(SA_chck_cast,t_chck_cast,p_chck_cast);
adiabatic_lapse_rate_from_t_error = nan(45,9);
for I = 1:9
    adiabatic_lapse_rate_from_t_error(:,I) = abs(adiabatic_lapse_rate_from_t(:,I) - adiabatic_lapse_rate_from_t(:,I+3));
end
gsw_cv.adiabatic_lapse_rate_from_t_ca = nanmax(nanmax(adiabatic_lapse_rate_from_t_error));

molality = gsw_molality_from_SA(SA_chck_cast);
molality_error = nan(45,9);
for I = 1:9
    molality_error(:,I) = abs(molality(:,I) - molality(:,I+3));
end
gsw_cv.molality_from_SA_ca = nanmax(nanmax(molality_error));
  
ionic_strength = gsw_ionic_strength_from_SA(SA_chck_cast);
ionic_strength_error = nan(45,9);
for I = 1:9
    ionic_strength_error(:,I) = abs(ionic_strength(:,I) - ionic_strength(:,I+3));
end
gsw_cv.ionic_strength_from_SA_ca = nanmax(nanmax(ionic_strength_error));
  
%% specific volume, density and enthalpy 

specvol = gsw_specvol(SA_chck_cast,CT_chck_cast,p_chck_cast);
specvol_error = nan(45,9);
for I = 1:9
    specvol_error(:,I) = abs(specvol(:,I) - specvol(:,I+3));
end
gsw_cv.specvol_ca = nanmax(nanmax(specvol_error));
 
alpha = gsw_alpha(SA_chck_cast,CT_chck_cast,p_chck_cast);
alpha_error = nan(45,9);
for I = 1:9
    alpha_error(:,I) = abs(alpha(:,I) - alpha(:,I+3));
end
gsw_cv.alpha_ca = nanmax(nanmax(alpha_error));

beta = gsw_beta(SA_chck_cast,CT_chck_cast,p_chck_cast);
beta_error = nan(45,9);
for I = 1:9
    beta_error(:,I) = abs(beta(:,I) - beta(:,I+3));
end
gsw_cv.beta_ca = nanmax(nanmax(beta_error));

[v_vab, alpha_vab, beta_vab] = gsw_specvol_alpha_beta(SA_chck_cast,CT_chck_cast,p_chck_cast);
v_vab_error = nan(45,9);
alpha_vab_error = nan(45,9);
beta_vab_error = nan(45,9);
for I = 1:9
    v_vab_error(:,I) = abs(v_vab(:,I) - v_vab(:,I+3));
end
gsw_cv.v_vab_ca = nanmax(nanmax(v_vab_error));
 for I = 1:9
    alpha_vab_error(:,I) = abs(alpha_vab(:,I) - alpha_vab(:,I+3));
end
gsw_cv.alpha_vab_ca = nanmax(nanmax(alpha_vab_error));
 for I = 1:9
    beta_vab_error(:,I) = abs(beta_vab(:,I) - beta_vab(:,I+3));
end
gsw_cv.beta_vab_ca = nanmax(nanmax(beta_vab_error));
  
alpha_on_beta = gsw_alpha_on_beta(SA_chck_cast,CT_chck_cast,p_chck_cast);
alpha_on_beta_error = nan(45,9);
for I = 1:9
    alpha_on_beta_error(:,I) = abs(alpha_on_beta(:,I) - alpha_on_beta(:,I+3));
end
gsw_cv.alpha_on_beta_ca = nanmax(nanmax(alpha_on_beta_error));

[v_SA, v_CT, v_P] = gsw_specvol_first_derivatives(SA_chck_cast,CT_chck_cast,p_chck_cast);
v_SA_error = nan(45,9);
v_CT_error = nan(45,9);
v_P_error = nan(45,9);
for I = 1:9
    v_SA_error(:,I) = abs(v_SA(:,I) - v_SA(:,I+3));
end
gsw_cv.v_SA_ca = nanmax(nanmax(v_SA_error));
 for I = 1:9
    v_CT_error(:,I) = abs(v_CT(:,I) - v_CT(:,I+3));
end
gsw_cv.v_CT_ca = nanmax(nanmax(v_CT_error));
 for I = 1:9
    v_P_error(:,I) = abs(v_P(:,I) - v_P(:,I+3));
end
gsw_cv.v_P_ca = nanmax(nanmax(v_P_error));

[v_SA_SA, v_SA_CT, v_CT_CT, v_SA_P, v_CT_P] = gsw_specvol_second_derivatives(SA_chck_cast,CT_chck_cast,p_chck_cast);
v_SA_SA_error = nan(45,9);
v_SA_CT_error = nan(45,9);
v_CT_CT_error = nan(45,9);
v_SA_P_error = nan(45,9);
v_CT_P_error = nan(45,9);
for I = 1:9
    v_SA_SA_error(:,I) = abs(v_SA_SA(:,I) - v_SA_SA(:,I+3));
end
gsw_cv.v_SA_SA_ca = nanmax(nanmax(v_SA_SA_error));
for I = 1:9
    v_SA_CT_error(:,I) = abs(v_SA_CT(:,I) - v_SA_CT(:,I+3));
end
gsw_cv.v_SA_CT_ca = nanmax(nanmax(v_SA_CT_error));
 for I = 1:9
    v_CT_CT_error(:,I) = abs(v_CT_CT(:,I) - v_CT_CT(:,I+3));
end
gsw_cv.v_CT_CT_ca = nanmax(nanmax(v_CT_CT_error));
 for I = 1:9
    v_SA_P_error(:,I) = abs(v_SA_P(:,I) - v_SA_P(:,I+3));
end
gsw_cv.v_SA_P_ca = nanmax(nanmax(v_SA_P_error));
for I = 1:9
    v_CT_P_error(:,I) = abs(v_CT_P(:,I) - v_CT_P(:,I+3));
end
gsw_cv.v_CT_P_ca = nanmax(nanmax(v_CT_P_error));

[v_SA_wrt_h, v_h] = gsw_specvol_first_derivatives_wrt_enthalpy(SA_chck_cast,CT_chck_cast,p_chck_cast);
v_SA_wrt_h_error = nan(45,9);
v_h_error = nan(45,9);
for I = 1:9
    v_SA_wrt_h_error(:,I) = abs(v_SA_wrt_h(:,I) - v_SA_wrt_h(:,I+3));
end
gsw_cv.v_SA_wrt_h_ca = nanmax(nanmax(v_SA_wrt_h_error));
 for I = 1:9
    v_h_error(:,I) = abs(v_h(:,I) - v_h(:,I+3));
end
gsw_cv.v_h_ca = nanmax(nanmax(v_h_error));

[v_SA_SA_wrt_h, v_SA_h, v_h_h] = gsw_specvol_second_derivatives_wrt_enthalpy(SA_chck_cast,CT_chck_cast,p_chck_cast);
v_SA_SA_wrt_h_error = nan(45,9);
v_SA_h_error = nan(45,9);
v_h_h_error = nan(45,9);
for I = 1:9
    v_SA_SA_wrt_h_error(:,I) = abs(v_SA_SA_wrt_h(:,I) - v_SA_SA_wrt_h(:,I+3));
end
gsw_cv.v_SA_SA_wrt_h_ca = nanmax(nanmax(v_SA_SA_wrt_h_error));
for I = 1:9
    v_SA_h_error(:,I) = abs(v_SA_h(:,I) - v_SA_h(:,I+3));
end
gsw_cv.v_SA_h_ca = nanmax(nanmax(v_SA_h_error));
 for I = 1:9
    v_h_h_error(:,I) = abs(v_h_h(:,I) - v_h_h(:,I+3));
end
gsw_cv.v_h_h_ca = nanmax(nanmax(v_h_h_error));
specvol_anom = gsw_specvol_anom(SA_chck_cast,CT_chck_cast,p_chck_cast,SA_ref,CT_ref);
specvol_anom_error = nan(45,9);
for I = 1:9
    specvol_anom_error(:,I) = abs(specvol_anom(:,I) - specvol_anom(:,I+3));
end
gsw_cv.specvol_anom_ca = nanmax(nanmax(specvol_anom_error));

specvol_anom_standard = gsw_specvol_anom_standard(SA_chck_cast,CT_chck_cast,p_chck_cast);
specvol_anom_standard_error = nan(45,9);
for I = 1:9
    specvol_anom_standard_error(:,I) = abs(specvol_anom_standard(:,I) - specvol_anom_standard(:,I+3));
end
gsw_cv.specvol_anom_standard_ca = nanmax(nanmax(specvol_anom_standard_error));
  
rho = gsw_rho(SA_chck_cast,CT_chck_cast,p_chck_cast);
rho_error = nan(45,9);
for I = 1:9
    rho_error(:,I) = abs(rho(:,I) - rho(:,I+3));
end
gsw_cv.rho_ca = nanmax(nanmax(rho_error));

[rho_rab, alpha_rab, beta_rab] = gsw_rho_alpha_beta(SA_chck_cast,CT_chck_cast,p_chck_cast);
rho_rab_error = nan(45,9);
alpha_rab_error = nan(45,9);
beta_rab_error = nan(45,9);
for I = 1:9
    rho_rab_error(:,I) = abs(rho_rab(:,I) - rho_rab(:,I+3));
end
gsw_cv.rho_rab_ca = nanmax(nanmax(rho_rab_error));
 for I = 1:9
    alpha_rab_error(:,I) = abs(alpha_rab(:,I) - alpha_rab(:,I+3));
end
gsw_cv.alpha_rab_ca = nanmax(nanmax(alpha_rab_error));
 for I = 1:9
    beta_rab_error(:,I) = abs(beta_rab(:,I) - beta_rab(:,I+3));
end
gsw_cv.beta_rab_ca = nanmax(nanmax(beta_rab_error));

[rho_SA, rho_CT, rho_P] = gsw_rho_first_derivatives(SA_chck_cast,CT_chck_cast,p_chck_cast);
rho_SA_error = nan(45,9);
rho_CT_error = nan(45,9);
rho_P_error = nan(45,9);
for I = 1:9
    rho_SA_error(:,I) = abs(rho_SA(:,I) - rho_SA(:,I+3));
end
gsw_cv.rho_SA_ca = nanmax(nanmax(rho_SA_error));
 for I = 1:9
    rho_CT_error(:,I) = abs(rho_CT(:,I) - rho_CT(:,I+3));
end
gsw_cv.rho_CT_ca = nanmax(nanmax(rho_CT_error));
 for I = 1:9
    rho_P_error(:,I) = abs(rho_P(:,I) - rho_P(:,I+3));
end
gsw_cv.rho_P_ca = nanmax(nanmax(rho_P_error));

[rho_SA_SA, rho_SA_CT, rho_CT_CT, rho_SA_P, rho_CT_P] = gsw_rho_second_derivatives(SA_chck_cast,CT_chck_cast,p_chck_cast);
rho_SA_SA_error = nan(45,9);
rho_SA_CT_error = nan(45,9);
rho_CT_CT_error = nan(45,9);
rho_SA_P_error = nan(45,9);
rho_CT_P_error = nan(45,9);
for I = 1:9
    rho_SA_SA_error(:,I) = abs(rho_SA_SA(:,I) - rho_SA_SA(:,I+3));
end
gsw_cv.rho_SA_SA_ca = nanmax(nanmax(rho_SA_SA_error));
for I = 1:9
    rho_SA_CT_error(:,I) = abs(rho_SA_CT(:,I) - rho_SA_CT(:,I+3));
end
gsw_cv.rho_SA_CT_ca = nanmax(nanmax(rho_SA_CT_error));
 for I = 1:9
    rho_CT_CT_error(:,I) = abs(rho_CT_CT(:,I) - rho_CT_CT(:,I+3));
end
gsw_cv.rho_CT_CT_ca = nanmax(nanmax(rho_CT_CT_error));
 for I = 1:9
    rho_SA_P_error(:,I) = abs(rho_SA_P(:,I) - rho_SA_P(:,I+3));
end
gsw_cv.rho_SA_P_ca = nanmax(nanmax(rho_SA_P_error));
for I = 1:9
    rho_CT_P_error(:,I) = abs(rho_CT_P(:,I) - rho_CT_P(:,I+3));
end
gsw_cv.rho_CT_P_ca = nanmax(nanmax(rho_CT_P_error));

[rho_SA_wrt_h, rho_h] = gsw_rho_first_derivatives_wrt_enthalpy(SA_chck_cast,CT_chck_cast,p_chck_cast);
rho_SA_wrt_h_error = nan(45,9);
rho_h_error = nan(45,9);
for I = 1:9
    rho_SA_wrt_h_error(:,I) = abs(rho_SA_wrt_h(:,I) - rho_SA_wrt_h(:,I+3));
end
gsw_cv.rho_SA_wrt_h_ca = nanmax(nanmax(rho_SA_wrt_h_error));
 for I = 1:9
    rho_h_error(:,I) = abs(rho_h(:,I) - rho_h(:,I+3));
end
gsw_cv.rho_h_ca = nanmax(nanmax(rho_h_error));

[rho_SA_SA_wrt_h, rho_SA_h, rho_h_h] = gsw_rho_second_derivatives_wrt_enthalpy(SA_chck_cast,CT_chck_cast,p_chck_cast);
rho_SA_SA_wrt_h_error = nan(45,9);
rho_SA_h_error = nan(45,9);
rho_h_h_error = nan(45,9);
for I = 1:9
    rho_SA_SA_wrt_h_error(:,I) = abs(rho_SA_SA_wrt_h(:,I) - rho_SA_SA_wrt_h(:,I+3));
end
gsw_cv.rho_SA_SA_wrt_h_ca = nanmax(nanmax(rho_SA_SA_wrt_h_error));
for I = 1:9
    rho_SA_h_error(:,I) = abs(rho_SA_h(:,I) - rho_SA_h(:,I+3));
end
gsw_cv.rho_SA_h_ca = nanmax(nanmax(rho_SA_h_error));
 for I = 1:9
    rho_h_h_error(:,I) = abs(rho_h_h(:,I) - rho_h_h(:,I+3));
end
gsw_cv.rho_h_h_ca = nanmax(nanmax(rho_h_h_error));

sigma0 = gsw_sigma0(SA_chck_cast,CT_chck_cast);
sigma0_error = nan(45,9);
for I = 1:9
    sigma0_error(:,I) = abs(sigma0(:,I) - sigma0(:,I+3));
end
gsw_cv.sigma0_ca = nanmax(nanmax(sigma0_error));
 
sigma1 = gsw_sigma1(SA_chck_cast,CT_chck_cast);
sigma1_error = nan(45,9);
for I = 1:9
    sigma1_error(:,I) = abs(sigma1(:,I) - sigma1(:,I+3));
end
gsw_cv.sigma1_ca = nanmax(nanmax(sigma1_error));
 
sigma2 = gsw_sigma2(SA_chck_cast,CT_chck_cast);
sigma2_error = nan(45,9);
for I = 1:9
    sigma2_error(:,I) = abs(sigma2(:,I) - sigma2(:,I+3));
end
gsw_cv.sigma2_ca = nanmax(nanmax(sigma2_error));
 
sigma3 = gsw_sigma3(SA_chck_cast,CT_chck_cast);
sigma3_error = nan(45,9);
for I = 1:9
    sigma3_error(:,I) = abs(sigma3(:,I) - sigma3(:,I+3));
end
gsw_cv.sigma3_ca = nanmax(nanmax(sigma3_error));
 
sigma4 = gsw_sigma4(SA_chck_cast,CT_chck_cast);
sigma4_error = nan(45,9);
for I = 1:9
    sigma4_error(:,I) = abs(sigma4(:,I) - sigma4(:,I+3));
end
gsw_cv.sigma4_ca = nanmax(nanmax(sigma4_error));

cabbeling = gsw_cabbeling(SA_chck_cast,CT_chck_cast,p_chck_cast);
cabbeling_error = nan(45,9);
for I = 1:9
    cabbeling_error(:,I) = abs(cabbeling(:,I) - cabbeling(:,I+3));
end
gsw_cv.cabbeling_ca = nanmax(nanmax(cabbeling_error));
  
thermobaric = gsw_thermobaric(SA_chck_cast,CT_chck_cast,p_chck_cast);
thermobaric_error = nan(45,9);
for I = 1:9
    thermobaric_error(:,I) = abs(thermobaric(:,I) - thermobaric(:,I+3));
end
gsw_cv.thermobaric_ca = nanmax(nanmax(thermobaric_error));

enthalpy = gsw_enthalpy(SA_chck_cast,CT_chck_cast,p_chck_cast);
enthalpy_error = nan(45,9);
for I = 1:9
    enthalpy_error(:,I) = abs(enthalpy(:,I) - enthalpy(:,I+3));
end
gsw_cv.enthalpy_ca = nanmax(nanmax(enthalpy_error));
  
enthalpy_diff =  gsw_enthalpy_diff(SA_chck_cast,CT_chck_cast,p_chck_cast_shallow,p_chck_cast_deep);
enthalpy_diff_error = nan(45,9);
for I = 1:9
    enthalpy_diff_error(:,I) = abs(enthalpy_diff(:,I) - enthalpy_diff(:,I+3));
end
gsw_cv.enthalpy_diff_ca = nanmax(nanmax(enthalpy_diff_error));
  
dynamic_enthalpy =  gsw_dynamic_enthalpy(SA_chck_cast,CT_chck_cast,p_chck_cast);
dynamic_enthalpy_error = nan(45,9);
for I = 1:9
    dynamic_enthalpy_error(:,I) = abs(dynamic_enthalpy(:,I) - dynamic_enthalpy(:,I+3));
end
gsw_cv.dynamic_enthalpy_ca = nanmax(nanmax(dynamic_enthalpy_error));

[h_SA, h_CT] = gsw_enthalpy_first_derivatives(SA_chck_cast,CT_chck_cast,p_chck_cast);
h_SA_error = nan(45,9);
h_CT_error = nan(45,9);
for I = 1:9
    h_SA_error(:,I) = abs(h_SA(:,I) - h_SA(:,I+3));
end
gsw_cv.h_SA_ca = nanmax(nanmax(h_SA_error));
 for I = 1:9
    h_CT_error(:,I) = abs(h_CT(:,I) - h_CT(:,I+3));
end
gsw_cv.h_CT_ca = nanmax(nanmax(h_CT_error));

[h_SA_SA, h_SA_CT, h_CT_CT] = gsw_enthalpy_second_derivatives(SA_chck_cast,CT_chck_cast,p_chck_cast);
h_SA_SA_error = nan(45,9);
h_SA_CT_error = nan(45,9);
h_CT_CT_error = nan(45,9);
for I = 1:9
    h_SA_SA_error(:,I) = abs(h_SA_SA(:,I) - h_SA_SA(:,I+3));
end
gsw_cv.h_SA_SA_ca = nanmax(nanmax(h_SA_SA_error));
 for I = 1:9
    h_SA_CT_error(:,I) = abs(h_SA_CT(:,I) - h_SA_CT(:,I+3));
end
gsw_cv.h_SA_CT_ca = nanmax(nanmax(h_SA_CT_error));
 for I = 1:9
    h_CT_CT_error(:,I) = abs(h_CT_CT(:,I) - h_CT_CT(:,I+3));
end
gsw_cv.h_CT_CT_ca = nanmax(nanmax(h_CT_CT_error));

sound_speed = gsw_sound_speed(SA_chck_cast,CT_chck_cast,p_chck_cast);
sound_speed_error = nan(45,9);
for I = 1:9
    sound_speed_error(:,I) = abs(sound_speed(:,I) - sound_speed(:,I+3));
end
gsw_cv.sound_speed_ca = nanmax(nanmax(sound_speed_error));

kappa = gsw_kappa(SA_chck_cast,CT_chck_cast,p_chck_cast);
kappa_error = nan(45,9);
for I = 1:9
    kappa_error(:,I) = abs(kappa(:,I) - kappa(:,I+3));
end
gsw_cv.kappa_ca = nanmax(nanmax(kappa_error));

internal_energy = gsw_internal_energy(SA_chck_cast,CT_chck_cast,p_chck_cast);
internal_energy_error = nan(45,9);
for I = 1:9
    internal_energy_error(:,I) = abs(internal_energy(:,I) - internal_energy(:,I+3));
end
gsw_cv.internal_energy_ca = nanmax(nanmax(internal_energy_error));

[u_SA, u_CT, u_P] = gsw_internal_energy_first_derivatives(SA_chck_cast,CT_chck_cast,p_chck_cast);
u_SA_error = nan(45,9);
u_CT_error = nan(45,9);
u_P_error = nan(45,9);
for I = 1:9
    u_SA_error(:,I) = abs(u_SA(:,I) - u_SA(:,I+3));
end
gsw_cv.u_SA_ca = nanmax(nanmax(u_SA_error));
 for I = 1:9
    u_CT_error(:,I) = abs(u_CT(:,I) - u_CT(:,I+3));
end
gsw_cv.u_CT_ca = nanmax(nanmax(u_CT_error));
 for I = 1:9
    u_P_error(:,I) = abs(u_P(:,I) - u_P(:,I+3));
end
gsw_cv.u_P_ca = nanmax(nanmax(u_P_error));

[u_SA_SA, u_SA_CT, u_CT_CT, u_SA_P, u_CT_P] = gsw_internal_energy_second_derivatives(SA_chck_cast,CT_chck_cast,p_chck_cast);
u_SA_SA_error = nan(45,9);
u_SA_CT_error = nan(45,9);
u_CT_CT_error = nan(45,9);
u_SA_P_error = nan(45,9);
u_CT_P_error = nan(45,9);
for I = 1:9
    u_SA_SA_error(:,I) = abs(u_SA_SA(:,I) - u_SA_SA(:,I+3));
end
gsw_cv.u_SA_SA_ca = nanmax(nanmax(u_SA_SA_error));
for I = 1:9
    u_SA_CT_error(:,I) = abs(u_SA_CT(:,I) - u_SA_CT(:,I+3));
end
gsw_cv.u_SA_CT_ca = nanmax(nanmax(u_SA_CT_error));
 for I = 1:9
    u_CT_CT_error(:,I) = abs(u_CT_CT(:,I) - u_CT_CT(:,I+3));
end
gsw_cv.u_CT_CT_ca = nanmax(nanmax(u_CT_CT_error));
 for I = 1:9
    u_SA_P_error(:,I) = abs(u_SA_P(:,I) - u_SA_P(:,I+3));
end
gsw_cv.u_SA_P_ca = nanmax(nanmax(u_SA_P_error));
for I = 1:9
    u_CT_P_error(:,I) = abs(u_CT_P(:,I) - u_CT_P(:,I+3));
end
gsw_cv.u_CT_P_ca = nanmax(nanmax(u_CT_P_error));

CT_from_enthalpy = gsw_CT_from_enthalpy(SA_chck_cast,enthalpy,p_chck_cast);
CT_from_enthalpy_error = nan(45,9);
for I = 1:9
    CT_from_enthalpy_error(:,I) = abs(CT_from_enthalpy(:,I) - CT_from_enthalpy(:,I+3));
end
gsw_cv.CT_from_enthalpy_ca = nanmax(nanmax(CT_from_enthalpy_error));

SA_from_rho_CT = gsw_SA_from_rho(rho,CT_chck_cast,p_chck_cast);
SA_from_rho_CT_error = nan(45,9);
for I = 1:9
    SA_from_rho_CT_error(:,I) = abs(SA_from_rho_CT(:,I) - SA_from_rho_CT(:,I+3));
end
gsw_cv.SA_from_rho_ca = nanmax(nanmax(SA_from_rho_CT_error));

CT_from_rho = gsw_CT_from_rho(rho,SA_chck_cast,p_chck_cast);
CT_from_rho_error = nan(45,9);
for I = 1:9
    CT_from_rho_error(:,I) = abs(CT_from_rho(:,I) - CT_from_rho(:,I+3));
end
gsw_cv.CT_from_rho_ca = nanmax(nanmax(CT_from_rho_error));
CT_maxdensity = gsw_CT_maxdensity(SA_chck_cast,p_chck_cast);
CT_maxdensity_error = nan(45,9);
 for I = 1:9
    CT_maxdensity_error(:,I) = abs(CT_maxdensity(:,I) - CT_maxdensity(:,I+3));
end
gsw_cv.CT_maxdensity_ca = nanmax(nanmax(CT_maxdensity_error));

%% vertical stability and interpolation

[Tu, Rsubrho, p_mid_TuRsr] = gsw_Turner_Rsubrho(SA_chck_cast,CT_chck_cast,p_chck_cast);
Tu_error = nan(44,9);
Rsubrho_error = nan(44,9);
p_mid_TuRsr_error= nan(44,9);
for I = 1:9
    Tu_error(:,I) = abs(Tu(:,I) - Tu(:,I+3));
end
gsw_cv.Tu_ca = nanmax(nanmax(Tu_error));
 for I = 1:9
    Rsubrho_error(:,I) = abs(Rsubrho(:,I) - Rsubrho(:,I+3));
end
gsw_cv.Rsubrho_ca = nanmax(nanmax(Rsubrho_error));
 for I = 1:9
    p_mid_TuRsr_error(:,I) = abs(p_mid_TuRsr(:,I) - p_mid_TuRsr(:,I+3));
end
gsw_cv.p_mid_TuRsr_ca = nanmax(nanmax(p_mid_TuRsr_error));

[n2, p_mid_n2] = gsw_Nsquared(SA_chck_cast,CT_chck_cast,p_chck_cast,lat_chck_cast);
n2_error = nan(44,9);
p_mid_n2_error = nan(44,9);
for I = 1:9
    n2_error(:,I) = abs(n2(:,I) - n2(:,I+3));
end
gsw_cv.n2_ca = nanmax(nanmax(n2_error));
 for I = 1:9
    p_mid_n2_error(:,I) = abs(p_mid_n2(:,I) - p_mid_n2(:,I+3));
end
gsw_cv.p_mid_n2_ca = nanmax(nanmax(p_mid_n2_error));
  
[n2min, n2min_pmid, n2min_specvol, n2min_alpha, n2min_beta, n2min_dsa, n2min_dct, n2min_dp] = gsw_Nsquared_min(SA_chck_cast,CT_chck_cast,p_chck_cast,lat_chck_cast);
n2min_error = nan(44,9);
n2min_pmid_error = nan(44,9);
n2min_specvol_error = nan(44,9);
n2min_alpha_error = nan(44,9);
n2min_beta_error = nan(44,9);
n2min_dsa_error = nan(44,9);
n2min_dct_error = nan(44,9);
n2min_dp_error = nan(44,9);
for I = 1:9
    n2min_error(:,I) = abs(n2min(:,I) - n2min(:,I+3));
end
gsw_cv.n2min_ca = nanmax(nanmax(n2min_error));
 for I = 1:9
    n2min_pmid_error(:,I) = abs(n2min_pmid(:,I) - n2min_pmid(:,I+3));
end
gsw_cv.n2min_pmid_ca = nanmax(nanmax(n2min_pmid_error));
for I = 1:9
    n2min_specvol_error(:,I) = abs(n2min_specvol(:,I) - n2min_specvol(:,I+3));
end
gsw_cv.n2min_specvol_ca = nanmax(nanmax(n2min_specvol_error));
 for I = 1:9
    n2min_alpha_error(:,I) = abs(n2min_alpha(:,I) - n2min_alpha(:,I+3));
end
gsw_cv.n2min_alpha_ca = nanmax(nanmax(n2min_alpha_error));
 for I = 1:9
    n2min_beta_error(:,I) = abs(n2min_beta(:,I) - n2min_beta(:,I+3));
end
gsw_cv.n2min_beta_ca = nanmax(nanmax(n2min_beta_error));
 for I = 1:9
    n2min_dsa_error(:,I) = abs(n2min_dsa(:,I) - n2min_dsa(:,I+3));
end
gsw_cv.n2min_dsa_ca = nanmax(nanmax(n2min_dsa_error));
 for I = 1:9
    n2min_dct_error(:,I) = abs(n2min_dct(:,I) - n2min_dct(:,I+3));
end
gsw_cv.n2min_dct_ca = nanmax(nanmax(n2min_dct_error));
 for I = 1:9
    n2min_dp_error(:,I) = abs(n2min_dp(:,I) - n2min_dp(:,I+3));
end
gsw_cv.n2min_dp_ca = nanmax(nanmax(n2min_dp_error));

mlp = gsw_mlp(SA_chck_cast,CT_chck_cast,p_chck_cast);
mlp_error = nan(9,1);
for I = 1:9
    mlp_error(I) = abs(mlp(I) - mlp(I+3));
end
gsw_cv.mlp_ca = nanmax(nanmax(mlp_error));

n2_lowerlimit = gsw_Nsquared_lowerlimit(p_chck_cast,long_chck_cast,lat_chck_cast);
n2_lowerlimit_error = nan(45,9);
for I = 1:9
    n2_lowerlimit_error(:,I) = abs(n2_lowerlimit(:,I) - n2_lowerlimit(:,I+3));
end
gsw_cv.n2_lowerlimit_ca = nanmax(nanmax(n2_lowerlimit_error));
  

[SAi_SACTinterp, CTi_SACTinterp] = gsw_SA_CT_interp(SA_chck_cast,CT_chck_cast,p_chck_cast,p_i);
SAi_SACTinterp_error = nan(51,9);
CTi_SACTinterp_error = nan(51,9);
for I = 1:9
    SAi_SACTinterp_error(:,I) = abs(SAi_SACTinterp(:,I) - SAi_SACTinterp(:,I+3));
end
gsw_cv.SAi_SACTinterp_ca = nanmax(nanmax(SAi_SACTinterp_error));
 for I = 1:9
    CTi_SACTinterp_error(:,I) = abs(CTi_SACTinterp(:,I) - CTi_SACTinterp(:,I+3));
end
gsw_cv.CTi_SACTinterp_ca = nanmax(nanmax(CTi_SACTinterp_error));

ti_tinterp = gsw_t_interp(CT_chck_cast,p_chck_cast,p_i);
ti_tinterp_error = nan(51,9);
for I = 1:9
    ti_tinterp_error(:,I) = abs(ti_tinterp(:,I) - ti_tinterp(:,I+3));
end
gsw_cv.ti_tinterp_ca = nanmax(nanmax(ti_tinterp_error));

[traceri_tracerCTinterp, CTi_tracerCTinterp] = gsw_tracer_CT_interp(SA_chck_cast,CT_chck_cast,p_chck_cast,p_i,9);
traceri_tracerCTinterp_error = nan(51,9);
CTi_tracerCTinterp_error = nan(51,9);
for I = 1:9
    traceri_tracerCTinterp_error(:,I) = abs(traceri_tracerCTinterp(:,I) - traceri_tracerCTinterp(:,I+3));
end
gsw_cv.traceri_tracerCTinterp_ca = nanmax(nanmax(traceri_tracerCTinterp_error));
 for I = 1:9
    CTi_tracerCTinterp_error(:,I) = abs(CTi_tracerCTinterp(:,I) - CTi_tracerCTinterp(:,I+3));
end
gsw_cv.CTi_tracerCTinterp_ca = nanmax(nanmax(CTi_tracerCTinterp_error));

traceri_tracerinterp = gsw_tracer_interp(CT_chck_cast,p_chck_cast,p_i);
traceri_tracerinterp_error = nan(51,9);
for I = 1:9
    traceri_tracerinterp_error(:,I) = abs(traceri_tracerinterp(:,I) - traceri_tracerinterp(:,I+3));
end
gsw_cv.traceri_tracerinterp_ca = nanmax(nanmax(traceri_tracerinterp_error));


[IPVfN2, p_mid_IPVfN2] = gsw_IPV_vs_fNsquared_ratio(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
IPVfN2_error = nan(44,9);
p_mid_IPVfN2_error = nan(44,9);
for I = 1:9
    IPVfN2_error(:,I) = abs(IPVfN2(:,I) - IPVfN2(:,I+3));
end
gsw_cv.IPVfN2_ca = nanmax(nanmax(IPVfN2_error));
 for I = 1:9
    p_mid_IPVfN2_error(:,I) = abs(p_mid_IPVfN2(:,I) - p_mid_IPVfN2(:,I+3));
end
gsw_cv.p_mid_IPVfN2_ca = nanmax(nanmax(p_mid_IPVfN2_error));
  
%% geostrophic streamfunctions, travel time and Geostrophic velocity

geo_strf_dyn_height = gsw_geo_strf_dyn_height(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
geo_strf_dyn_height_error = nan(45,9);
for I = 1:9
    geo_strf_dyn_height_error(:,I) = abs(geo_strf_dyn_height(:,I) - geo_strf_dyn_height(:,I+3));
end
gsw_cv.geo_strf_dyn_height_ca = nanmax(nanmax(geo_strf_dyn_height_error));
 
[geo_strf_dyn_height_pc, geo_strf_dyn_height_pc_p_mid] = gsw_geo_strf_dyn_height_pc(SA_chck_cast,CT_chck_cast,delta_p_chck_cast);
geo_strf_dyn_height_pc_error = nan(45,9);
geo_strf_dyn_height_pc_p_mid_error = nan(45,9);
for I = 1:9
    geo_strf_dyn_height_pc_error(:,I) = abs(geo_strf_dyn_height_pc(:,I) - geo_strf_dyn_height_pc(:,I+3));
end
gsw_cv.geo_strf_dyn_height_pc_ca = nanmax(nanmax(geo_strf_dyn_height_pc_error));
 for I = 1:9
    geo_strf_dyn_height_pc_p_mid_error(:,I) = abs(geo_strf_dyn_height_pc_p_mid(:,I) - geo_strf_dyn_height_pc_p_mid(:,I+3));
end
gsw_cv.geo_strf_dyn_height_pc_p_mid_ca = 1e-15;
  
geo_strf_isopycnal = gsw_geo_strf_isopycnal(SA_chck_cast,CT_chck_cast,p_chck_cast,pr,Neutral_Density,p_Neutral_Density);
geo_strf_isopycnal_error = nan(45,9);
for I = 1:9
    geo_strf_isopycnal_error(:,I) = abs(geo_strf_isopycnal(:,I) - geo_strf_isopycnal(:,I+3));
end
gsw_cv.geo_strf_isopycnal_ca = nanmax(nanmax(geo_strf_isopycnal_error));
  
[geo_strf_isopycnal_pc, geo_strf_isopycnal_pc_p_mid] = gsw_geo_strf_isopycnal_pc(SA_chck_cast,CT_chck_cast,delta_p_chck_cast,26.8,3);
geo_strf_isopycnal_pc_error = nan(45,9);
mk_p_mid_error = nan(45,9);
for I = 1:9
    geo_strf_isopycnal_pc_error(:,I) = abs(geo_strf_isopycnal_pc(:,I) - geo_strf_isopycnal_pc(:,I+3));
end
gsw_cv.geo_strf_isopycnal_pc_ca = nanmax(nanmax(geo_strf_isopycnal_pc_error));
 for I = 1:9
    mk_p_mid_error(:,I) = abs(geo_strf_isopycnal_pc_p_mid(:,I) - geo_strf_isopycnal_pc_p_mid(:,I+3));
end
gsw_cv.geo_strf_isopycnal_pc_p_mid_ca = 1e-15;
  
geo_strf_Montgomery = gsw_geo_strf_Montgomery(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
geo_strf_Montgomery_error = nan(45,9);
for I = 1:9
    geo_strf_Montgomery_error(:,I) = abs(geo_strf_Montgomery(:,I) - geo_strf_Montgomery(:,I+3));
end
gsw_cv.geo_strf_Montgomery_ca = nanmax(nanmax(geo_strf_Montgomery_error));
  
geo_strf_Cunningham = gsw_geo_strf_Cunningham(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
geo_strf_Cunningham_error = nan(45,9);
for I = 1:9
    geo_strf_Cunningham_error(:,I) = abs(geo_strf_Cunningham(:,I) - geo_strf_Cunningham(:,I+3));
end
gsw_cv.geo_strf_Cunningham_ca = nanmax(nanmax(geo_strf_Cunningham_error));
  
steric_height = gsw_geo_strf_steric_height(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
steric_height_error = nan(45,8);
for I = 1:8
    steric_height_error(:,I) = abs(steric_height(:,I) - steric_height(:,I+3));
end
gsw_cv.geo_strf_steric_height_ca = nanmax(nanmax(steric_height_error));

PISH = gsw_geo_strf_PISH(SA_chck_cast,CT_chck_cast,p_chck_cast,pr_05);
PISH_error = nan(45,8);
for I = 1:8
    PISH_error(:,I) = abs(PISH(:,I) - PISH(:,I+3));
end
gsw_cv.geo_strf_PISH_ca = nanmax(nanmax(PISH_error));

travel_time = gsw_travel_time(SA_chck_cast,CT_chck_cast,p_chck_cast,lat_chck_cast(1));
travel_time_error = nan(45,8);
for I = 1:8
    travel_time_error(:,I) = abs(travel_time(:,I) - travel_time(:,I+3));
end
gsw_cv.travel_time_ca = nanmax(nanmax(travel_time_error));

[geo_strf_velocity, geo_strf_velocity_mid_lat, geo_strf_velocity_mid_long] = gsw_geostrophic_velocity(geo_strf_dyn_height,long_chck_cast,lat_chck_cast,p_chck_cast);
geo_strf_velocity_error = nan(45,8);
geo_strf_velocity_mid_lat_error = nan(45,8);
geo_strf_velocity_mid_long_error = nan(45,8);
for I = 1:8
    geo_strf_velocity_error(:,I) = abs(geo_strf_velocity(:,I) - geo_strf_velocity(:,I+3));
end
gsw_cv.geo_strf_velocity_ca = nanmax(nanmax(geo_strf_velocity_error));
 for I = 1:8
    geo_strf_velocity_mid_lat_error(:,I) = abs(geo_strf_velocity_mid_lat(:,I) - geo_strf_velocity_mid_lat(:,I+3));
end
gsw_cv.geo_strf_velocity_mid_lat_ca = 1e-15;
 for I = 1:8
    geo_strf_velocity_mid_long_error(:,I) = abs(geo_strf_velocity_mid_long(:,I) - geo_strf_velocity_mid_long(:,I+3));
end
gsw_cv.geo_strf_velocity_mid_long_ca = 1e-15;
        
%% neutral versus isopycnal slopes and ratios 

isopycnal_slope_ratio = gsw_isopycnal_slope_ratio(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
isopycnal_slope_ratio_error = nan(45,9);
for I = 1:9
    isopycnal_slope_ratio_error(:,I) = abs(isopycnal_slope_ratio(:,I) - isopycnal_slope_ratio(:,I+3));
end
gsw_cv.isopycnal_slope_ratio_ca = nanmax(nanmax(isopycnal_slope_ratio_error));
  
[G_CT, p_mid_G_CT] = gsw_isopycnal_vs_ntp_CT_ratio(SA_chck_cast,CT_chck_cast,p_chck_cast,pr);
G_CT_error = nan(44,9);
p_mid_G_CT_error = nan(44,9);
for I = 1:9
    G_CT_error(:,I) = abs(G_CT(:,I) - G_CT(:,I+3));
end
gsw_cv.G_CT_ca = nanmax(nanmax(G_CT_error));
 for I = 1:9
    p_mid_G_CT_error(:,I) = abs(p_mid_G_CT(:,I) - p_mid_G_CT(:,I+3));
end
gsw_cv.p_mid_G_CT_ca = nanmax(nanmax(p_mid_G_CT_error));
  
ntpptCT = gsw_ntp_pt_vs_CT_ratio(SA_chck_cast,CT_chck_cast,p_chck_cast);
ntpptCT_error = nan(45,9);
for I = 1:9
    ntpptCT_error(:,I) = abs(ntpptCT(:,I) - ntpptCT(:,I+3));
end
gsw_cv.ntpptCT_ca = nanmax(nanmax(ntpptCT_error));
  
%% derivatives of entropy, CT and pt

[CT_SA, CT_pt] = gsw_CT_first_derivatives(SA_chck_cast,pt);
CT_SA_error = nan(45,9);
CT_pt_error = nan(45,9);
for I = 1:9
    CT_SA_error(:,I) = abs(CT_SA(:,I) - CT_SA(:,I+3));
end
gsw_cv.CT_SA_ca = nanmax(nanmax(CT_SA_error));
 for I = 1:9
    CT_pt_error(:,I) = abs(CT_pt(:,I) - CT_pt(:,I+3));
end
gsw_cv.CT_pt_ca = nanmax(nanmax(CT_pt_error));
  
[CT_SA_SA, CT_SA_pt, CT_pt_pt] = gsw_CT_second_derivatives(SA_chck_cast,pt);
CT_SA_SA_error = nan(45,9);
CT_SA_pt_error = nan(45,9);
CT_pt_pt_error = nan(45,9);
for I = 1:9
    CT_SA_SA_error(:,I) = abs(CT_SA_SA(:,I) - CT_SA_SA(:,I+3));
end
gsw_cv.CT_SA_SA_ca = nanmax(nanmax(CT_SA_SA_error));
 for I = 1:9
    CT_SA_pt_error(:,I) = abs(CT_SA_pt(:,I) - CT_SA_pt(:,I+3));
end
gsw_cv.CT_SA_pt_ca = nanmax(nanmax(CT_SA_pt_error));
 for I = 1:9
    CT_pt_pt_error(:,I) = abs(CT_pt_pt(:,I) - CT_pt_pt(:,I+3));
end
gsw_cv.CT_pt_pt_ca = nanmax(nanmax(CT_pt_pt_error));
  
[eta_SA, eta_CT] = gsw_entropy_first_derivatives(SA_chck_cast,CT_chck_cast);
eta_SA_error = nan(45,9);
eta_CT_error = nan(45,9);
for I = 1:9
    eta_SA_error(:,I) = abs(eta_SA(:,I) - eta_SA(:,I+3));
end
gsw_cv.eta_SA_ca = nanmax(nanmax(eta_SA_error));
 for I = 1:9
    eta_CT_error(:,I) = abs(eta_CT(:,I) - eta_CT(:,I+3));
end
gsw_cv.eta_CT_ca = nanmax(nanmax(eta_CT_error));
  
[eta_SA_SA, eta_SA_CT, eta_CT_CT] = gsw_entropy_second_derivatives(SA_chck_cast,CT_chck_cast);
eta_SA_SA_error = nan(45,9);
eta_SA_CT_error = nan(45,9);
eta_CT_CT_error = nan(45,9);
for I = 1:9
    eta_SA_SA_error(:,I) = abs(eta_SA_SA(:,I) - eta_SA_SA(:,I+3));
end
gsw_cv.eta_SA_SA_ca = nanmax(nanmax(eta_SA_SA_error));
 for I = 1:9
    eta_SA_CT_error(:,I) = abs(eta_SA_CT(:,I) - eta_SA_CT(:,I+3));
end
gsw_cv.eta_SA_CT_ca = nanmax(nanmax(eta_SA_CT_error));
 for I = 1:9
    eta_CT_CT_error(:,I) = abs(eta_CT_CT(:,I) - eta_CT_CT(:,I+3));
end
gsw_cv.eta_CT_CT_ca = nanmax(nanmax(eta_CT_CT_error));
 
[pt_SA, pt_CT] = gsw_pt_first_derivatives(SA_chck_cast,CT_chck_cast);
pt_SA_error = nan(45,9);
pt_CT_error = nan(45,9);
for I = 1:9
    pt_SA_error(:,I) = abs(pt_SA(:,I) - pt_SA(:,I+3));
end
gsw_cv.pt_SA_ca = nanmax(nanmax(pt_SA_error));
 for I = 1:9
    pt_CT_error(:,I) = abs(pt_CT(:,I) - pt_CT(:,I+3));
end
gsw_cv.pt_CT_ca = nanmax(nanmax(pt_CT_error));
  
[pt_SA_SA, pt_SA_CT, pt_CT_CT] = gsw_pt_second_derivatives(SA_chck_cast,CT_chck_cast);
pt_SA_SA_error = nan(45,9);
pt_SA_CT_error = nan(45,9);
pt_CT_CT_error = nan(45,9);
for I = 1:9
    pt_SA_SA_error(:,I) = abs(pt_SA_SA(:,I) - pt_SA_SA(:,I+3));
end
gsw_cv.pt_SA_SA_ca = nanmax(nanmax(pt_SA_SA_error));
 for I = 1:9
    pt_SA_CT_error(:,I) = abs(pt_SA_CT(:,I) - pt_SA_CT(:,I+3));
end
gsw_cv.pt_SA_CT_ca = nanmax(nanmax(pt_SA_CT_error));
 for I = 1:9
    pt_CT_CT_error(:,I) = abs(pt_CT_CT(:,I) - pt_CT_CT(:,I+3));
end
gsw_cv.pt_CT_CT_ca = nanmax(nanmax(pt_CT_CT_error));
  
%% seawater properties at freezing temperatures

CT_freezing = gsw_CT_freezing(SA_chck_cast,p_chck_cast,0.5);
CT_freezing_error = nan(45,9);
for I = 1:9
    CT_freezing_error(:,I) = abs(CT_freezing(:,I) - CT_freezing(:,I+3));
end
gsw_cv.CT_freezing_ca =  nanmax(nanmax(CT_freezing_error));

CT_freezing_poly = gsw_CT_freezing_poly(SA_chck_cast,p_chck_cast,0.5);
CT_freezing_poly_error = nan(45,9);
for I = 1:9
    CT_freezing_poly_error(:,I) = abs(CT_freezing_poly(:,I) - CT_freezing_poly(:,I+3));
end
gsw_cv.CT_freezing_poly_ca =  nanmax(nanmax(CT_freezing_poly_error));

t_freezing = gsw_t_freezing(SA_chck_cast,p_chck_cast,0.5);
t_freezing_error = nan(45,9);
for I = 1:9
    t_freezing_error(:,I) = abs(t_freezing(:,I) - t_freezing(:,I+3));
end
gsw_cv.t_freezing_ca =  nanmax(nanmax(t_freezing_error));

t_freezing_poly = gsw_t_freezing_poly(SA_chck_cast,p_chck_cast,0.5);
t_freezing_poly_error = nan(45,9);
for I = 1:9
    t_freezing_poly_error(:,I) = abs(t_freezing_poly(:,I) - t_freezing_poly(:,I+3));
end
gsw_cv.t_freezing_poly_ca =  nanmax(nanmax(t_freezing_poly_error));

pot_enthalpy_ice_freezing = gsw_pot_enthalpy_ice_freezing(SA_chck_cast,p_chck_cast);
pot_enthalpy_ice_freezing_error = nan(45,9);
for I = 1:9
    pot_enthalpy_ice_freezing_error(:,I) = abs(pot_enthalpy_ice_freezing(:,I) - pot_enthalpy_ice_freezing(:,I+3));
end
gsw_cv.pot_enthalpy_ice_freezing_ca =  nanmax(nanmax(pot_enthalpy_ice_freezing_error));

pot_enthalpy_ice_freezing_poly = gsw_pot_enthalpy_ice_freezing_poly(SA_chck_cast,p_chck_cast);
pot_enthalpy_ice_freezing_poly_error = nan(45,9);
for I = 1:9
    pot_enthalpy_ice_freezing_poly_error(:,I) = abs(pot_enthalpy_ice_freezing_poly(:,I) - pot_enthalpy_ice_freezing_poly(:,I+3));
end
gsw_cv.pot_enthalpy_ice_freezing_poly_ca =  nanmax(nanmax(pot_enthalpy_ice_freezing_poly_error));

SA_freezing_from_CT = gsw_SA_freezing_from_CT(CT_freezing,p_chck_cast,0.5);
SA_freezing_from_CT_error = nan(45,9);
for I = 1:9
    SA_freezing_from_CT_error(:,I) = abs(SA_freezing_from_CT(:,I) - SA_freezing_from_CT(:,I+3));
end
gsw_cv.SA_freezing_from_CT_ca =  nanmax(nanmax(SA_freezing_from_CT_error));

SA_freezing_from_CT_poly = gsw_SA_freezing_from_CT_poly(CT_freezing_poly,p_chck_cast,0.5);
SA_freezing_from_CT_poly_error = nan(45,9);
for I = 1:9
    SA_freezing_from_CT_poly_error(:,I) = abs(SA_freezing_from_CT_poly(:,I) - SA_freezing_from_CT_poly(:,I+3));
end
gsw_cv.SA_freezing_from_CT_poly_ca =  nanmax(nanmax(SA_freezing_from_CT_poly_error));

SA_freezing_from_t = gsw_SA_freezing_from_t(t_freezing,p_chck_cast,0.5);
SA_freezing_from_t_error = nan(45,9);
for I = 1:9
    SA_freezing_from_t_error(:,I) = abs(SA_freezing_from_t(:,I) - SA_freezing_from_t(:,I+3));
end
gsw_cv.SA_freezing_from_t_ca = nanmax(nanmax(SA_freezing_from_t_error));

SA_freezing_from_t_poly = gsw_SA_freezing_from_t_poly(t_freezing_poly,p_chck_cast,0.5);
SA_freezing_from_t_poly_error = nan(45,9);
for I = 1:9
    SA_freezing_from_t_poly_error(:,I) = abs(SA_freezing_from_t_poly(:,I) - SA_freezing_from_t_poly(:,I+3));
end
gsw_cv.SA_freezing_from_t_poly_ca = nanmax(nanmax(SA_freezing_from_t_poly_error));

pressure_freezing_CT = gsw_pressure_freezing_CT(SA_Arctic,(CT_Arctic - 1),0.5);
pressure_freezing_CT_error = nan(36,9);
for I = 1:9
    pressure_freezing_CT_error(:,I) = abs(pressure_freezing_CT(:,I) - pressure_freezing_CT(:,I+3));
end
gsw_cv.pressure_freezing_CT_ca = nanmax(nanmax(pressure_freezing_CT_error));

[CTfreezing_SA, CTfreezing_P] = gsw_CT_freezing_first_derivatives(SA_chck_cast,p_chck_cast,0.5);
CTfreezing_SA_error = nan(45,9);
for I = 1:9
    CTfreezing_SA_error(:,I) = abs(CTfreezing_SA(:,I) - CTfreezing_SA(:,I+3));
end
gsw_cv.CTfreezing_SA_ca = nanmax(nanmax(CTfreezing_SA_error));
CTfreezing_P_error = nan(45,9);
for I = 1:9
    CTfreezing_P_error(:,I) = abs(CTfreezing_P(:,I) - CTfreezing_P(:,I+3));
end
gsw_cv.CTfreezing_P_ca = nanmax(nanmax(CTfreezing_P_error));

[CTfreezing_SA_poly, CTfreezing_P_poly] = gsw_CT_freezing_first_derivatives_poly(SA_chck_cast,p_chck_cast,0.5);
CTfreezing_SA_poly_error = nan(45,9);
for I = 1:9
    CTfreezing_SA_poly_error(:,I) = abs(CTfreezing_SA_poly(:,I) - CTfreezing_SA_poly(:,I+3));
end
gsw_cv.CTfreezing_SA_poly_ca = nanmax(nanmax(CTfreezing_SA_poly_error));
CTfreezing_P_poly_error = nan(45,9);
for I = 1:9
    CTfreezing_P_poly_error(:,I) = abs(CTfreezing_P_poly(:,I) - CTfreezing_P_poly(:,I+3));
end
gsw_cv.CTfreezing_P_poly_ca = nanmax(nanmax(CTfreezing_P_poly_error));

[tfreezing_SA, tfreezing_P] = gsw_t_freezing_first_derivatives(SA_chck_cast,p_chck_cast,0.5);
tfreezing_SA_error = nan(45,9);
for I = 1:9
    tfreezing_SA_error(:,I) = abs(tfreezing_SA(:,I) - tfreezing_SA(:,I+3));
end
gsw_cv.tfreezing_SA_ca = nanmax(nanmax(tfreezing_SA_error));
tfreezing_P_error = nan(45,9);
for I = 1:9
    tfreezing_P_error(:,I) = abs(tfreezing_P(:,I) - tfreezing_P(:,I+3));
end
gsw_cv.tfreezing_P_ca = nanmax(nanmax(tfreezing_P_error));

[tfreezing_SA_poly, tfreezing_P_poly] = gsw_t_freezing_first_derivatives_poly(SA_chck_cast,p_chck_cast,0.5);
tfreezing_SA_poly_error = nan(45,9);
for I = 1:9
    tfreezing_SA_poly_error(:,I) = abs(tfreezing_SA_poly(:,I) - tfreezing_SA_poly(:,I+3));
end
gsw_cv.tfreezing_SA_poly_ca = nanmax(nanmax(tfreezing_SA_poly_error));
tfreezing_P_poly_error = nan(45,9);
for I = 1:9
    tfreezing_P_poly_error(:,I) = abs(tfreezing_P_poly(:,I) - tfreezing_P_poly(:,I+3));
end
gsw_cv.tfreezing_P_poly_ca = nanmax(nanmax(tfreezing_P_poly_error));

[pot_enthalpy_ice_freezing_SA, pot_enthalpy_ice_freezing_P] = gsw_pot_enthalpy_ice_freezing_first_derivatives(SA_chck_cast,p_chck_cast);
pot_enthalpy_ice_freezing_SA_error = nan(45,9);
for I = 1:9
    pot_enthalpy_ice_freezing_SA_error(:,I) = abs(pot_enthalpy_ice_freezing_SA(:,I) - pot_enthalpy_ice_freezing_SA(:,I+3));
end
gsw_cv.pot_enthalpy_ice_freezing_SA_ca = nanmax(nanmax(pot_enthalpy_ice_freezing_SA_error));
pot_enthalpy_ice_freezing_P_error = nan(45,9);
for I = 1:9
    pot_enthalpy_ice_freezing_P_error(:,I) = abs(pot_enthalpy_ice_freezing_P(:,I) - pot_enthalpy_ice_freezing_P(:,I+3));
end
gsw_cv.pot_enthalpy_ice_freezing_P_ca = nanmax(nanmax(pot_enthalpy_ice_freezing_P_error));

[pot_enthalpy_ice_freezing_SA_poly, pot_enthalpy_ice_freezing_P_poly] = gsw_pot_enthalpy_ice_freezing_first_derivatives_poly(SA_chck_cast,p_chck_cast);
pot_enthalpy_ice_freezing_SA_poly_error = nan(45,9);
for I = 1:9
    pot_enthalpy_ice_freezing_SA_poly_error(:,I) = abs(pot_enthalpy_ice_freezing_SA_poly(:,I) - pot_enthalpy_ice_freezing_SA_poly(:,I+3));
end
gsw_cv.pot_enthalpy_ice_freezing_SA_poly_ca = nanmax(nanmax(pot_enthalpy_ice_freezing_SA_poly_error));
pot_enthalpy_ice_freezing_P_poly_error = nan(45,9);
for I = 1:9
    pot_enthalpy_ice_freezing_P_poly_error(:,I) = abs(pot_enthalpy_ice_freezing_P_poly(:,I) - pot_enthalpy_ice_freezing_P_poly(:,I+3));
end
gsw_cv.pot_enthalpy_ice_freezing_P_poly_ca = nanmax(nanmax(pot_enthalpy_ice_freezing_P_poly_error));

latentheat_melting = gsw_latentheat_melting(SA_chck_cast,p_chck_cast);
latentheat_melting_error = nan(45,9);
for I = 1:9
    latentheat_melting_error(:,I) = abs(latentheat_melting(:,I) - latentheat_melting(:,I+3));
end
gsw_cv.latentheat_melting_ca =  nanmax(nanmax(latentheat_melting_error));


%%  thermodynamic interaction between ice and seawater

CT_Arctic = gsw_CT_from_t(SA_Arctic,t_Arctic,p_Arctic);

melting_ice_SA_CT_ratio = gsw_melting_ice_SA_CT_ratio(SA_Arctic,CT_Arctic,p_Arctic,t_seaice);
melting_ice_SA_CT_ratio_error = nan(36,9);
for I = 1:9
    melting_ice_SA_CT_ratio_error(:,I) = abs(melting_ice_SA_CT_ratio(:,I) - melting_ice_SA_CT_ratio(:,I+3));
end
gsw_cv.melting_ice_SA_CT_ratio_ca = nanmax(nanmax(melting_ice_SA_CT_ratio_error));

melting_ice_SA_CT_ratio_poly = gsw_melting_ice_SA_CT_ratio_poly(SA_Arctic,CT_Arctic,p_Arctic,t_seaice);
melting_ice_SA_CT_ratio_poly_error = nan(36,9);
for I = 1:9
    melting_ice_SA_CT_ratio_poly_error(:,I) = abs(melting_ice_SA_CT_ratio_poly(:,I) - melting_ice_SA_CT_ratio_poly(:,I+3));
end
gsw_cv.melting_ice_SA_CT_ratio_poly_ca = nanmax(nanmax(melting_ice_SA_CT_ratio_poly_error));

melting_ice_equilibrium_SA_CT_ratio = gsw_melting_ice_equilibrium_SA_CT_ratio(SA_Arctic,p_Arctic);
melting_ice_equilibrium_SA_CT_ratio_error = nan(36,9);
for I = 1:9
    melting_ice_equilibrium_SA_CT_ratio_error(:,I) = abs(melting_ice_equilibrium_SA_CT_ratio(:,I) - melting_ice_equilibrium_SA_CT_ratio(:,I+3));
end
gsw_cv.melting_ice_equilibrium_SA_CT_ratio_ca = nanmax(nanmax(melting_ice_equilibrium_SA_CT_ratio_error));

melting_ice_equilibrium_SA_CT_ratio_poly = gsw_melting_ice_equilibrium_SA_CT_ratio_poly(SA_Arctic,p_Arctic);
melting_ice_equilibrium_SA_CT_ratio_poly_error = nan(36,9);
for I = 1:9
    melting_ice_equilibrium_SA_CT_ratio_poly_error(:,I) = abs(melting_ice_equilibrium_SA_CT_ratio_poly(:,I) - melting_ice_equilibrium_SA_CT_ratio_poly(:,I+3));
end
gsw_cv.melting_ice_equilibrium_SA_CT_ratio_poly_ca = nanmax(nanmax(melting_ice_equilibrium_SA_CT_ratio_poly_error));

[melting_ice_into_seawater_SA_final, melting_ice_into_seawater_CT_final] = gsw_melting_ice_into_seawater(SA_Arctic,CT_Arctic+0.1,p_Arctic,w_seaice,t_seaice);
melting_ice_into_seawater_SA_final_error = nan(36,9);
for I = 1:9
    melting_ice_into_seawater_SA_final_error(:,I) = abs(melting_ice_into_seawater_SA_final(:,I) - melting_ice_into_seawater_SA_final(:,I+3));
end
gsw_cv.melting_ice_into_seawater_SA_final_ca = nanmax(nanmax(melting_ice_into_seawater_SA_final_error));
melting_ice_into_seawater_CT_final_error = nan(36,9);
for I = 1:9
    melting_ice_into_seawater_CT_final_error(:,I) = abs(melting_ice_into_seawater_CT_final(:,I) - melting_ice_into_seawater_CT_final(:,I+3));
end
gsw_cv.melting_ice_into_seawater_CT_final_ca = nanmax(nanmax(melting_ice_into_seawater_CT_final_error));

[ice_fraction_to_freeze_seawater_SA_freeze, ice_fraction_to_freeze_seawater_CT_freeze, ice_fraction_to_freeze_seawater_w_Ih] = gsw_ice_fraction_to_freeze_seawater(SA_Arctic,CT_Arctic,p_Arctic,t_seaice);
ice_fraction_to_freeze_seawater_SA_freeze_error = nan(36,9);
for I = 1:9
    ice_fraction_to_freeze_seawater_SA_freeze_error(:,I) = abs(ice_fraction_to_freeze_seawater_SA_freeze(:,I) - ice_fraction_to_freeze_seawater_SA_freeze(:,I+3));
end
gsw_cv.ice_fraction_to_freeze_seawater_SA_freeze_ca = nanmax(nanmax(ice_fraction_to_freeze_seawater_SA_freeze_error));
ice_fraction_to_freeze_seawater_CT_freeze_error = nan(36,9);
for I = 1:9
    ice_fraction_to_freeze_seawater_CT_freeze_error(:,I) = abs(ice_fraction_to_freeze_seawater_CT_freeze(:,I) - ice_fraction_to_freeze_seawater_CT_freeze(:,I+3));
end
gsw_cv.ice_fraction_to_freeze_seawater_CT_freeze_ca = nanmax(nanmax(ice_fraction_to_freeze_seawater_CT_freeze_error));
ice_fraction_to_freeze_seawater_w_Ih_error = nan(36,9);
for I = 1:9
    ice_fraction_to_freeze_seawater_w_Ih_error(:,I) = abs(ice_fraction_to_freeze_seawater_w_Ih(:,I) - ice_fraction_to_freeze_seawater_w_Ih(:,I+3));
end
gsw_cv.ice_fraction_to_freeze_seawater_w_Ih_ca = nanmax(nanmax(ice_fraction_to_freeze_seawater_w_Ih_error));
 
[dSA_dCT_frazil, dSA_dP_frazil, dCT_dP_frazil] = gsw_frazil_ratios_adiabatic(SA_Arctic,p_Arctic,w_seaice);
dSA_dCT_frazil_error = nan(36,9);
for I = 1:9
    dSA_dCT_frazil_error(:,I) = abs(dSA_dCT_frazil(:,I) - dSA_dCT_frazil(:,I+3));
end
gsw_cv.dSA_dCT_frazil_ca = nanmax(nanmax(dSA_dCT_frazil_error));
dSA_dP_frazil_error = nan(36,9);
for I = 1:9
    dSA_dP_frazil_error(:,I) = abs(dSA_dP_frazil(:,I) - dSA_dP_frazil(:,I+3));
end
gsw_cv.dSA_dP_frazil_ca = nanmax(nanmax(dSA_dP_frazil_error));
dCT_dP_frazil_error = nan(36,9);
for I = 1:9
    dCT_dP_frazil_error(:,I) = abs(dCT_dP_frazil(:,I) - dCT_dP_frazil(:,I+3));
end
gsw_cv.dCT_dP_frazil_ca = nanmax(nanmax(dCT_dP_frazil_error));

[dSA_dCT_frazil_poly, dSA_dP_frazil_poly, dCT_dP_frazil_poly] = gsw_frazil_ratios_adiabatic_poly(SA_Arctic,p_Arctic,w_seaice);
dSA_dCT_frazil_poly_error = nan(36,9);
for I = 1:9
    dSA_dCT_frazil_poly_error(:,I) = abs(dSA_dCT_frazil_poly(:,I) - dSA_dCT_frazil_poly(:,I+3));
end
gsw_cv.dSA_dCT_frazil_poly_ca = nanmax(nanmax(dSA_dCT_frazil_poly_error));
dSA_dP_frazil_poly_error = nan(36,9);
for I = 1:9
    dSA_dP_frazil_poly_error(:,I) = abs(dSA_dP_frazil_poly(:,I) - dSA_dP_frazil_poly(:,I+3));
end
gsw_cv.dSA_dP_frazil_poly_ca = nanmax(nanmax(dSA_dP_frazil_poly_error));
dCT_dP_frazil_poly_error = nan(36,9);
for I = 1:9
    dCT_dP_frazil_poly_error(:,I) = abs(dCT_dP_frazil_poly(:,I) - dCT_dP_frazil_poly(:,I+3));
end
gsw_cv.dCT_dP_frazil_poly_ca = nanmax(nanmax(dCT_dP_frazil_poly_error));

[frazil_properties_potential_SA_final, frazil_properties_potential_CT_final, frazil_properties_potential_w_Ih_final] = gsw_frazil_properties_potential(SA_bulk,h_pot_bulk,p_Arctic);
frazil_properties_potential_SA_final_error = nan(36,9);
for I = 1:9
    frazil_properties_potential_SA_final_error(:,I) = abs(frazil_properties_potential_SA_final(:,I) - frazil_properties_potential_SA_final(:,I+3));
end
gsw_cv.frazil_properties_potential_SA_final_ca = nanmax(nanmax(frazil_properties_potential_SA_final_error));
frazil_properties_potential_CT_final_error = nan(36,9);
for I = 1:9
    frazil_properties_potential_CT_final_error(:,I) = abs(frazil_properties_potential_CT_final(:,I) - frazil_properties_potential_CT_final(:,I+3));
end
gsw_cv.frazil_properties_potential_CT_final_ca = nanmax(nanmax(frazil_properties_potential_CT_final_error));
frazil_properties_potential_w_Ih_final_error = nan(36,9);
for I = 1:9
    frazil_properties_potential_w_Ih_final_error(:,I) = abs(frazil_properties_potential_w_Ih_final(:,I) - frazil_properties_potential_w_Ih_final(:,I+3));
end
gsw_cv.frazil_properties_potential_w_Ih_final_ca = nanmax(nanmax(frazil_properties_potential_w_Ih_final_error));

[frazil_properties_potential_poly_SA_final, frazil_properties_potential_poly_CT_final, frazil_properties_potential_poly_w_Ih_final] = gsw_frazil_properties_potential_poly(SA_bulk,h_pot_bulk,p_Arctic);
frazil_properties_potential_poly_SA_final_error = nan(36,9);
for I = 1:9
    frazil_properties_potential_poly_SA_final_error(:,I) = abs(frazil_properties_potential_poly_SA_final(:,I) - frazil_properties_potential_poly_SA_final(:,I+3));
end
gsw_cv.frazil_properties_potential_poly_SA_final_ca = nanmax(nanmax(frazil_properties_potential_poly_SA_final_error));
frazil_properties_potential_poly_CT_final_error = nan(36,9);
for I = 1:9
    frazil_properties_potential_poly_CT_final_error(:,I) = abs(frazil_properties_potential_poly_CT_final(:,I) - frazil_properties_potential_poly_CT_final(:,I+3));
end
gsw_cv.frazil_properties_potential_poly_CT_final_ca = nanmax(nanmax(frazil_properties_potential_poly_CT_final_error));
frazil_properties_potential_poly_w_Ih_final_error = nan(36,9);
for I = 1:9
    frazil_properties_potential_poly_w_Ih_final_error(:,I) = abs(frazil_properties_potential_poly_w_Ih_final(:,I) - frazil_properties_potential_poly_w_Ih_final(:,I+3));
end
gsw_cv.frazil_properties_potential_poly_w_Ih_final_ca = nanmax(nanmax(frazil_properties_potential_poly_w_Ih_final_error));

[frazil_properties_SA_final, frazil_properties_CT_final, frazil_properties_w_Ih_final] = gsw_frazil_properties(SA_bulk,h_bulk,p_Arctic);
frazil_properties_SA_final_error = nan(36,9);
for I = 1:9
    frazil_properties_SA_final_error(:,I) = abs(frazil_properties_SA_final(:,I) - frazil_properties_SA_final(:,I+3));
end
gsw_cv.frazil_properties_SA_final_ca = nanmax(nanmax(frazil_properties_SA_final_error));
frazil_properties_CT_final_error = nan(36,9);
for I = 1:9
    frazil_properties_CT_final_error(:,I) = abs(frazil_properties_CT_final(:,I) - frazil_properties_CT_final(:,I+3));
end
gsw_cv.frazil_properties_CT_final_ca = nanmax(nanmax(frazil_properties_CT_final_error));
frazil_properties_w_Ih_final_error = nan(36,9);
for I = 1:9
    frazil_properties_w_Ih_final_error(:,I) = abs(frazil_properties_w_Ih_final(:,I) - frazil_properties_w_Ih_final(:,I+3));
end
gsw_cv.frazil_properties_w_Ih_final_ca = nanmax(nanmax(frazil_properties_w_Ih_final_error));

%%  thermodynamic interaction between seaice and seawater

melting_seaice_SA_CT_ratio = gsw_melting_seaice_SA_CT_ratio(SA_Arctic,CT_Arctic,p_Arctic,SA_seaice,t_seaice);
melting_seaice_SA_CT_ratio_error = nan(36,9);
for I = 1:9
    melting_seaice_SA_CT_ratio_error(:,I) = abs(melting_seaice_SA_CT_ratio(:,I) - melting_seaice_SA_CT_ratio(:,I+3));
end
gsw_cv.melting_seaice_SA_CT_ratio_ca = nanmax(nanmax(melting_seaice_SA_CT_ratio_error));

melting_seaice_SA_CT_ratio_poly = gsw_melting_seaice_SA_CT_ratio_poly(SA_Arctic,CT_Arctic,p_Arctic,SA_seaice,t_seaice);
melting_seaice_SA_CT_ratio_poly_error = nan(36,9);
for I = 1:9
    melting_seaice_SA_CT_ratio_poly_error(:,I) = abs(melting_seaice_SA_CT_ratio_poly(:,I) - melting_seaice_SA_CT_ratio_poly(:,I+3));
end
gsw_cv.melting_seaice_SA_CT_ratio_poly_ca = nanmax(nanmax(melting_seaice_SA_CT_ratio_poly_error));

melting_seaice_equilibrium_SA_CT_ratio = gsw_melting_seaice_equilibrium_SA_CT_ratio(SA_Arctic,p_Arctic);
melting_seaice_equilibrium_SA_CT_ratio_error = nan(36,9);
for I = 1:9
    melting_seaice_equilibrium_SA_CT_ratio_error(:,I) = abs(melting_seaice_equilibrium_SA_CT_ratio(:,I) - melting_seaice_equilibrium_SA_CT_ratio(:,I+3));
end
gsw_cv.melting_seaice_equilibrium_SA_CT_ratio_ca = nanmax(nanmax(melting_seaice_equilibrium_SA_CT_ratio_error));

melting_seaice_equilibrium_SA_CT_ratio_poly = gsw_melting_seaice_equilibrium_SA_CT_ratio_poly(SA_Arctic,p_Arctic);
melting_seaice_equilibrium_SA_CT_ratio_poly_error = nan(36,9);
for I = 1:9
    melting_seaice_equilibrium_SA_CT_ratio_poly_error(:,I) = abs(melting_seaice_equilibrium_SA_CT_ratio_poly(:,I) - melting_seaice_equilibrium_SA_CT_ratio_poly(:,I+3));
end
gsw_cv.melting_seaice_equilibrium_SA_CT_ratio_poly_ca = nanmax(nanmax(melting_seaice_equilibrium_SA_CT_ratio_poly_error));

[melting_seaice_into_seawater_SA_final, melting_seaice_into_seawater_CT_final] = gsw_melting_seaice_into_seawater(SA_Arctic,CT_Arctic,p_Arctic,w_seaice,SA_seaice,t_seaice);
melting_seaice_into_seawater_SA_final_error = nan(36,9);
for I = 1:9
    melting_seaice_into_seawater_SA_final_error(:,I) = abs(melting_seaice_into_seawater_SA_final(:,I) - melting_seaice_into_seawater_SA_final(:,I+3));
end
gsw_cv.melting_seaice_into_seawater_SA_final_ca = nanmax(nanmax(melting_seaice_into_seawater_SA_final_error));
melting_seaice_into_seawater_CT_final_error = nan(36,9);
for I = 1:9
    melting_seaice_into_seawater_CT_final_error(:,I) = abs(melting_seaice_into_seawater_CT_final(:,I) - melting_seaice_into_seawater_CT_final(:,I+3));
end
gsw_cv.melting_seaice_into_seawater_CT_final_ca = nanmax(nanmax(melting_seaice_into_seawater_CT_final_error));

[seaice_fraction_to_freeze_seawater_SA_freeze, seaice_fraction_to_freeze_seawater_CT_freeze, seaice_fraction_to_freeze_seawater_w_Ih] = gsw_seaice_fraction_to_freeze_seawater(SA_Arctic,CT_Arctic,p_Arctic,SA_seaice,t_seaice);
seaice_fraction_to_freeze_seawater_SA_freeze_error = nan(36,9);
for I = 1:9
    seaice_fraction_to_freeze_seawater_SA_freeze_error(:,I) = abs(seaice_fraction_to_freeze_seawater_SA_freeze(:,I) - seaice_fraction_to_freeze_seawater_SA_freeze(:,I+3));
end
gsw_cv.seaice_fraction_to_freeze_seawater_SA_freeze_ca = nanmax(nanmax(seaice_fraction_to_freeze_seawater_SA_freeze_error));
seaice_fraction_to_freeze_seawater_CT_freeze_error = nan(36,9);
for I = 1:9
    seaice_fraction_to_freeze_seawater_CT_freeze_error(:,I) = abs(seaice_fraction_to_freeze_seawater_CT_freeze(:,I) - seaice_fraction_to_freeze_seawater_CT_freeze(:,I+3));
end
gsw_cv.seaice_fraction_to_freeze_seawater_CT_freeze_ca = nanmax(nanmax(seaice_fraction_to_freeze_seawater_CT_freeze_error));
seaice_fraction_to_freeze_seawater_w_Ih_error = nan(36,9);
for I = 1:9
    seaice_fraction_to_freeze_seawater_w_Ih_error(:,I) = abs(seaice_fraction_to_freeze_seawater_w_Ih(:,I) - seaice_fraction_to_freeze_seawater_w_Ih(:,I+3));
end
gsw_cv.seaice_fraction_to_freeze_seawater_w_Ih_ca = nanmax(nanmax(seaice_fraction_to_freeze_seawater_w_Ih_error));

%% themodynamic properties of ice Ih

rho_ice = gsw_rho_ice(t_ice,p_Arctic);
rho_ice_error = nan(36,9);
for I = 1:9
    rho_ice_error(:,I) = abs(rho_ice(:,I) - rho_ice(:,I+3));
end
gsw_cv.rho_ice_ca = nanmax(nanmax(rho_ice_error));

alpha_ice = gsw_alpha_wrt_t_ice(t_ice,p_Arctic);
alpha_ice_error = nan(36,9);
for I = 1:9
    alpha_ice_error(:,I) = abs(alpha_ice(:,I) - alpha_ice(:,I+3));
end
gsw_cv.alpha_wrt_t_ice_ca = nanmax(nanmax(alpha_ice_error));

specvol_ice = gsw_specvol_ice(t_ice,p_Arctic);
specvol_ice_error = nan(36,9);
for I = 1:9
    specvol_ice_error(:,I) = abs(specvol_ice(:,I) - specvol_ice(:,I+3));
end
gsw_cv.specvol_ice_ca = nanmax(nanmax(specvol_ice_error));

pressure_coefficient_ice  = gsw_pressure_coefficient_ice (t_ice,p_Arctic);
pressure_coefficient_ice_error = nan(36,9);
for I = 1:9
    pressure_coefficient_ice_error(:,I) = abs(pressure_coefficient_ice(:,I) - pressure_coefficient_ice(:,I+3));
end
gsw_cv.pressure_coefficient_ice_ca = nanmax(nanmax(pressure_coefficient_ice_error));

sound_speed_ice = gsw_sound_speed_ice(t_ice,p_Arctic);
sound_speed_ice_error = nan(36,9);
for I = 1:9
    sound_speed_ice_error(:,I) = abs(sound_speed_ice(:,I) - sound_speed_ice(:,I+3));
end
gsw_cv.sound_speed_ice_ca = nanmax(nanmax(sound_speed_ice_error));

kappa_ice = gsw_kappa_ice(t_ice,p_Arctic);
kappa_ice_error = nan(36,9);
for I = 1:9
    kappa_ice_error(:,I) = abs(kappa_ice(:,I) - kappa_ice(:,I+3));
end
gsw_cv.kappa_ice_ca = nanmax(nanmax(kappa_ice_error));

kappa_const_t_ice = gsw_kappa_const_t_ice(t_ice,p_Arctic);
kappa_const_t_ice_error = nan(36,9);
for I = 1:9
    kappa_const_t_ice_error(:,I) = abs(kappa_const_t_ice(:,I) - kappa_const_t_ice(:,I+3));
end
gsw_cv.kappa_const_t_ice_ca = nanmax(nanmax(kappa_const_t_ice_error));

internal_energy_ice = gsw_internal_energy_ice(t_ice,p_Arctic);
internal_energy_ice_error = nan(36,9);
for I = 1:9
    internal_energy_ice_error(:,I) = abs(internal_energy_ice(:,I) - internal_energy_ice(:,I+3));
end
gsw_cv.internal_energy_ice_ca = nanmax(nanmax(internal_energy_ice_error));

enthalpy_ice = gsw_enthalpy_ice(t_ice,p_Arctic);
enthalpy_ice_error = nan(36,9);
for I = 1:9
    enthalpy_ice_error(:,I) = abs(enthalpy_ice(:,I) - enthalpy_ice(:,I+3));
end
gsw_cv.enthalpy_ice_ca = nanmax(nanmax(enthalpy_ice_error));

entropy_ice = gsw_entropy_ice(t_ice,p_Arctic);
entropy_ice_error = nan(36,9);
for I = 1:9
    entropy_ice_error(:,I) = abs(entropy_ice(:,I) - entropy_ice(:,I+3));
end
gsw_cv.entropy_ice_ca = nanmax(nanmax(entropy_ice_error));

cp_ice = gsw_cp_ice(t_ice,p_Arctic);
cp_ice_error = nan(36,9);
for I = 1:9
    cp_ice_error(:,I) = abs(cp_ice(:,I) - cp_ice(:,I+3));
end
gsw_cv.cp_ice_ca = nanmax(nanmax(cp_ice_error));

chem_potential_water_ice =  gsw_chem_potential_water_ice(t_ice,p_Arctic);
chem_potential_water_ice_error = nan(36,9);
for I = 1:9
    chem_potential_water_ice_error(:,I) = abs(chem_potential_water_ice(:,I) - chem_potential_water_ice(:,I+3));
end
gsw_cv.chem_potential_water_ice_ca = nanmax(nanmax(chem_potential_water_ice_error));

Helmholtz_energy_ice = gsw_Helmholtz_energy_ice(t_ice,p_Arctic);
Helmholtz_energy_ice_error = nan(36,9);
for I = 1:9
    Helmholtz_energy_ice_error(:,I) = abs(Helmholtz_energy_ice(:,I) - Helmholtz_energy_ice(:,I+3));
end
gsw_cv.Helmholtz_energy_ice_ca = nanmax(nanmax(Helmholtz_energy_ice_error));

adiabatic_lapse_rate_ice = gsw_adiabatic_lapse_rate_ice(t_ice,p_Arctic);
adiabatic_lapse_rate_ice_error = nan(36,9);
for I = 1:9
    adiabatic_lapse_rate_ice_error(:,I) = abs(adiabatic_lapse_rate_ice(:,I) - adiabatic_lapse_rate_ice(:,I+3));
end
gsw_cv.adiabatic_lapse_rate_ice_ca = nanmax(nanmax(adiabatic_lapse_rate_ice_error));

pt0_from_t_ice = gsw_pt0_from_t_ice(t_ice,p_Arctic);
pt0_from_t_ice_error = nan(36,9);
for I = 1:9
    pt0_from_t_ice_error(:,I) = abs(pt0_from_t_ice(:,I) - pt0_from_t_ice(:,I+3));
end
gsw_cv.pt0_from_t_ice_ca = nanmax(nanmax(pt0_from_t_ice_error));

pt_from_t_ice = gsw_pt_from_t_ice(t_ice,p_Arctic,pr);
pt_from_t_ice_error = nan(36,9);
for I = 1:9
    pt_from_t_ice_error(:,I) = abs(pt_from_t_ice(:,I) - pt_from_t_ice(:,I+3));
end
gsw_cv.pt_from_t_ice_ca = nanmax(nanmax(pt_from_t_ice_error));

t_from_pt0_ice = gsw_t_from_pt0_ice(pt0_from_t_ice,p_Arctic);
t_from_pt0_ice_error = nan(36,9);
for I = 1:9
    t_from_pt0_ice_error(:,I) = abs(t_from_pt0_ice(:,I) - t_from_pt0_ice(:,I+3));
end
gsw_cv.t_from_pt0_ice_ca = nanmax(nanmax(t_from_pt0_ice_error));

t_from_rho_ice = gsw_t_from_rho_ice(rho_ice,p_Arctic);
t_from_rho_ice_error = nan(36,9);
for I = 1:9
    t_from_rho_ice_error(:,I) = abs(t_from_rho_ice(:,I) - t_from_rho_ice(:,I+3));
end
gsw_cv.t_from_rho_ice_ca = nanmax(nanmax(t_from_rho_ice_error));

pot_enthalpy_from_pt_ice = gsw_pot_enthalpy_from_pt_ice(pt_from_t_ice);
pot_enthalpy_from_pt_ice_error = nan(36,9);
for I = 1:9
    pot_enthalpy_from_pt_ice_error(:,I) = abs(pot_enthalpy_from_pt_ice(:,I) - pot_enthalpy_from_pt_ice(:,I+3));
end
gsw_cv.pot_enthalpy_from_pt_ice_ca = nanmax(nanmax(pot_enthalpy_from_pt_ice_error));

pt_from_pot_enthalpy_ice = gsw_pt_from_pot_enthalpy_ice(pot_enthalpy_from_pt_ice);
pt_from_pot_enthalpy_ice_error = nan(36,9);
for I = 1:9
    pt_from_pot_enthalpy_ice_error(:,I) = abs(pt_from_pot_enthalpy_ice(:,I) - pt_from_pot_enthalpy_ice(:,I+3));
end
gsw_cv.pt_from_pot_enthalpy_ice_ca = nanmax(nanmax(pt_from_pot_enthalpy_ice_error));

pot_enthalpy_from_pt_ice_poly = gsw_pot_enthalpy_from_pt_ice_poly(pt_from_t_ice);
pot_enthalpy_from_pt_ice_poly_error = nan(36,9);
for I = 1:9
    pot_enthalpy_from_pt_ice_poly_error(:,I) = abs(pot_enthalpy_from_pt_ice_poly(:,I) - pot_enthalpy_from_pt_ice_poly(:,I+3));
end
gsw_cv.pot_enthalpy_from_pt_ice_poly_ca = nanmax(nanmax(pot_enthalpy_from_pt_ice_poly_error));

pt_from_pot_enthalpy_ice_poly = gsw_pt_from_pot_enthalpy_ice_poly(pot_enthalpy_from_pt_ice_poly);
pt_from_pot_enthalpy_ice_poly_error = nan(36,9);
for I = 1:9
    pt_from_pot_enthalpy_ice_poly_error(:,I) = abs(pt_from_pot_enthalpy_ice_poly(:,I) - pt_from_pot_enthalpy_ice_poly(:,I+3));
end
gsw_cv.pt_from_pot_enthalpy_ice_poly_ca = nanmax(nanmax(pt_from_pot_enthalpy_ice_poly_error));

pot_enthalpy_from_specvol_ice = gsw_pot_enthalpy_from_specvol_ice(specvol_ice,p_Arctic);
pot_enthalpy_from_specvol_ice_error = nan(36,9);
for I = 1:9
    pot_enthalpy_from_specvol_ice_error(:,I) = abs(pot_enthalpy_from_specvol_ice(:,I) - pot_enthalpy_from_specvol_ice(:,I+3));
end
gsw_cv.pot_enthalpy_from_specvol_ice_ca = nanmax(nanmax(pot_enthalpy_from_specvol_ice_error));

specvol_from_pot_enthalpy_ice = gsw_specvol_from_pot_enthalpy_ice(pot_enthalpy_from_specvol_ice,p_Arctic);
specvol_from_pot_enthalpy_ice_error = nan(36,9);
for I = 1:9
    specvol_from_pot_enthalpy_ice_error(:,I) = abs(specvol_from_pot_enthalpy_ice(:,I) - specvol_from_pot_enthalpy_ice(:,I+3));
end
gsw_cv.specvol_from_pot_enthalpy_ice_ca = nanmax(nanmax(specvol_from_pot_enthalpy_ice_error));

pot_enthalpy_from_specvol_ice_poly = gsw_pot_enthalpy_from_specvol_ice_poly(specvol_ice,p_Arctic);
pot_enthalpy_from_specvol_ice_poly_error = nan(36,9);
for I = 1:9
    pot_enthalpy_from_specvol_ice_poly_error(:,I) = abs(pot_enthalpy_from_specvol_ice_poly(:,I) - pot_enthalpy_from_specvol_ice_poly(:,I+3));
end
gsw_cv.pot_enthalpy_from_specvol_ice_poly_ca = nanmax(nanmax(pot_enthalpy_from_specvol_ice_poly_error));

specvol_from_pot_enthalpy_ice_poly = gsw_specvol_from_pot_enthalpy_ice_poly(pot_enthalpy_from_specvol_ice_poly,p_Arctic);
specvol_from_pot_enthalpy_ice_poly_error = nan(36,9);
for I = 1:9
    specvol_from_pot_enthalpy_ice_poly_error(:,I) = abs(specvol_from_pot_enthalpy_ice_poly(:,I) - specvol_from_pot_enthalpy_ice_poly(:,I+3));
end
gsw_cv.specvol_from_pot_enthalpy_ice_poly_ca = nanmax(nanmax(specvol_from_pot_enthalpy_ice_poly_error));

%% isobaric evaporation enthalpy 

latentheat_evap_CT = gsw_latentheat_evap_CT(SA_chck_cast,CT_chck_cast);
latentheat_evap_CT_error = nan(45,9);
for I = 1:9
    latentheat_evap_CT_error(:,I) = abs(latentheat_evap_CT(:,I) - latentheat_evap_CT(:,I+3));
end
gsw_cv.latentheat_evap_CT_ca =  nanmax(nanmax(latentheat_evap_CT_error));

latentheat_evap_t = gsw_latentheat_evap_t(SA_chck_cast,t_chck_cast);
latentheat_evap_t_error = nan(45,9);
for I = 1:9
    latentheat_evap_t_error(:,I) = abs(latentheat_evap_t(:,I) - latentheat_evap_t(:,I+3));
end
gsw_cv.latentheat_evap_t_ca =  nanmax(nanmax(latentheat_evap_t_error));

%% spiciness

spiciness0 = gsw_spiciness0(SA_chck_cast,CT_chck_cast);
spiciness0_error = nan(45,9);
for I = 1:9
    spiciness0_error(:,I) = abs(spiciness0(:,I) - spiciness0(:,I+3));
end
gsw_cv.spiciness0_ca =  nanmax(nanmax(spiciness0_error));

spiciness1 = gsw_spiciness1(SA_chck_cast,CT_chck_cast);
spiciness1_error = nan(45,9);
for I = 1:9
    spiciness1_error(:,I) = abs(spiciness1(:,I) - spiciness1(:,I+3));
end
gsw_cv.spiciness1_ca =  nanmax(nanmax(spiciness1_error));

spiciness2 = gsw_spiciness2(SA_chck_cast,CT_chck_cast);
spiciness2_error = nan(45,9);
for I = 1:9
    spiciness2_error(:,I) = abs(spiciness2(:,I) - spiciness2(:,I+3));
end
gsw_cv.spiciness2_ca =  nanmax(nanmax(spiciness2_error));

[SA_spicsig0, CT_spicsig0] = gsw_SA_CT_from_sigma0_spiciness0(sigma0,spiciness0);
SA_spicsig0_error = nan(45,9);
for I = 1:9
    SA_spicsig0_error(:,I) = abs(SA_spicsig0(:,I) - SA_spicsig0(:,I+3));
end
gsw_cv.SA_spicsig0_ca =  nanmax(nanmax(SA_spicsig0_error));
CT_spicsig0_error = nan(45,9);
for I = 1:9
    CT_spicsig0_error(:,I) = abs(CT_spicsig0(:,I) - CT_spicsig0(:,I+3));
end
gsw_cv.CT_spicsig0_ca =  nanmax(nanmax(CT_spicsig0_error));

[SA_spicsig1, CT_spicsig1] = gsw_SA_CT_from_sigma1_spiciness1(sigma1,spiciness1);
SA_spicsig1_error = nan(45,9);
for I = 1:9
    SA_spicsig1_error(:,I) = abs(SA_spicsig1(:,I) - SA_spicsig1(:,I+3));
end
gsw_cv.SA_spicsig1_ca =  nanmax(nanmax(SA_spicsig1_error));
CT_spicsig1_error = nan(45,9);
for I = 1:9
    CT_spicsig1_error(:,I) = abs(CT_spicsig1(:,I) - CT_spicsig1(:,I+3));
end
gsw_cv.CT_spicsig1_ca =  nanmax(nanmax(CT_spicsig1_error));

[SA_spicsig2, CT_spicsig2] = gsw_SA_CT_from_sigma2_spiciness2(sigma2,spiciness2);
SA_spicsig2_error = nan(45,9);
for I = 1:9
    SA_spicsig2_error(:,I) = abs(SA_spicsig2(:,I) - SA_spicsig2(:,I+3));
end
gsw_cv.SA_spicsig2_ca =  nanmax(nanmax(SA_spicsig2_error));
CT_spicsig2_error = nan(45,9);
for I = 1:9
    CT_spicsig2_error(:,I) = abs(CT_spicsig2(:,I) - CT_spicsig2(:,I+3));
end
gsw_cv.CT_spicsig2_ca =  nanmax(nanmax(CT_spicsig2_error));

%% planet earth properties 

f = gsw_f(lat_chck_cast);
f_error = nan(45,9);
for I = 1:9
    f_error(:,I) = abs(f(:,I) - f(:,I+3));
end
gsw_cv.f_ca = 1e-15;
  
grav = gsw_grav(lat_chck_cast,p_chck_cast);
grav_error = nan(45,9);
for I = 1:9
    grav_error(:,I) = abs(grav(:,I) - grav(:,I+3));
end
gsw_cv.grav_ca = nanmax(nanmax(grav_error));
  
distance = gsw_distance(long_chck_cast,lat_chck_cast,p_chck_cast);
distance_error = nan(45,8);
for I = 1:8
    distance_error(:,I) = abs(distance(:,I) - distance(:,I+3));
end
gsw_cv.distance_ca = nanmax(nanmax(distance_error));


%% TEOS-10 constants

% gsw_cv.t0_ca = 0;
% 
% gsw_cv.p0_ca = 0;
% 
% gsw_cv.SSO_ca = 0;
% 
% gsw_cv.uPS_ca = 0;
% 
% gsw_cv.cp0_ca = 0;
% 
% gsw_cv.C3515_ca = 0;
% 
% gsw_cv.SonCl_ca = 0;
% 
% gsw_cv.valence_factor_ca = 0;
% 
% gsw_cv.atomic_weight_ca = 0;

%% dissolved gasses

Arsol = gsw_Arsol(SA_chck_cast,CT_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
Arsol_error = nan(45,9);
for I = 1:9
    Arsol_error(:,I) = abs(Arsol(:,I) - Arsol(:,I+3));
end
gsw_cv.Arsol_ca = nanmax(nanmax(Arsol_error));

Arsol_SP_pt = gsw_Arsol_SP_pt(SP_chck_cast,pt0);
Arsol_SP_pt_error = nan(45,9);
for I = 1:9
    Arsol_SP_pt_error(:,I) = abs(Arsol_SP_pt(:,I) - Arsol_SP_pt(:,I+3));
end
gsw_cv.Arsol_SP_pt_ca = nanmax(nanmax(Arsol_SP_pt_error));

Hesol = gsw_Hesol(SA_chck_cast,CT_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
Hesol_error = nan(45,9);
for I = 1:9
    Hesol_error(:,I) = abs(Hesol(:,I) - Hesol(:,I+3));
end
gsw_cv.Hesol_ca = nanmax(nanmax(Hesol_error));

Hesol_SP_pt = gsw_Hesol_SP_pt(SP_chck_cast,pt0);
Hesol_SP_pt_error = nan(45,9);
for I = 1:9
    Hesol_SP_pt_error(:,I) = abs(Hesol_SP_pt(:,I) - Hesol_SP_pt(:,I+3));
end
gsw_cv.Hesol_SP_pt_ca = nanmax(nanmax(Hesol_SP_pt_error));

Krsol = gsw_Krsol(SA_chck_cast,CT_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
Krsol_error = nan(45,9);
for I = 1:9
    Krsol_error(:,I) = abs(Krsol(:,I) - Krsol(:,I+3));
end
gsw_cv.Krsol_ca = nanmax(nanmax(Krsol_error));

Krsol_SP_pt = gsw_Krsol_SP_pt(SP_chck_cast,pt0);
Krsol_SP_pt_error = nan(45,9);
for I = 1:9
    Krsol_SP_pt_error(:,I) = abs(Krsol_SP_pt(:,I) - Arsol_SP_pt(:,I+3));
end
gsw_cv.Krsol_SP_pt_ca = nanmax(nanmax(Krsol_SP_pt_error));

N2Osol = gsw_N2Osol(SA_chck_cast,CT_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
N2Osol_error = nan(45,9);
for I = 1:9
    N2Osol_error(:,I) = abs(N2Osol(:,I) - N2Osol(:,I+3));
end
gsw_cv.N2Osol_ca = nanmax(nanmax(N2Osol_error));

N2Osol_SP_pt = gsw_N2Osol_SP_pt(SP_chck_cast,pt0);
N2Osol_SP_pt_error = nan(45,9);
for I = 1:9
    N2Osol_SP_pt_error(:,I) = abs(N2Osol_SP_pt(:,I) - N2Osol_SP_pt(:,I+3));
end
gsw_cv.N2Osol_SP_pt_ca = nanmax(nanmax(N2Osol_SP_pt_error));

N2sol = gsw_N2sol(SA_chck_cast,CT_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
N2sol_error = nan(45,9);
for I = 1:9
    N2sol_error(:,I) = abs(N2sol(:,I) - N2sol(:,I+3));
end
gsw_cv.N2sol_ca = nanmax(nanmax(N2sol_error));

N2sol_SP_pt = gsw_N2sol_SP_pt(SP_chck_cast,pt0);
N2sol_SP_pt_error = nan(45,9);
for I = 1:9
    N2sol_SP_pt_error(:,I) = abs(N2sol_SP_pt(:,I) - N2sol_SP_pt(:,I+3));
end
gsw_cv.N2sol_SP_pt_ca = nanmax(nanmax(N2sol_SP_pt_error));

Nesol = gsw_Nesol(SA_chck_cast,CT_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
Nesol_error = nan(45,9);
for I = 1:9
    Nesol_error(:,I) = abs(Nesol(:,I) - Nesol(:,I+3));
end
gsw_cv.Nesol_ca = nanmax(nanmax(Nesol_error));

Nesol_SP_pt = gsw_Nesol_SP_pt(SP_chck_cast,pt0);
Nesol_SP_pt_error = nan(45,9);
for I = 1:9
    Nesol_SP_pt_error(:,I) = abs(Nesol_SP_pt(:,I) - Nesol_SP_pt(:,I+3));
end
gsw_cv.Nesol_SP_pt_ca = nanmax(nanmax(Nesol_SP_pt_error));

O2sol = gsw_O2sol(SA_chck_cast,CT_chck_cast,p_chck_cast,long_chck_cast,lat_chck_cast);
O2sol_error = nan(45,9);
for I = 1:9
    O2sol_error(:,I) = abs(O2sol(:,I) - O2sol(:,I+3));
end
gsw_cv.O2sol_ca = nanmax(nanmax(O2sol_error));

O2sol_SP_pt = gsw_O2sol_SP_pt(SP_chck_cast,pt0);
O2sol_SP_pt_error = nan(45,9);
for I = 1:9
    O2sol_SP_pt_error(:,I) = abs(O2sol_SP_pt(:,I) - O2sol_SP_pt(:,I+3));
end
gsw_cv.O2sol_SP_pt_ca = nanmax(nanmax(O2sol_SP_pt_error));


%% density and enthalpy in terms of CT, derived from the exact Gibbs function

specvol_CT_exact = gsw_specvol_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
specvol_CT_exact_error = nan(45,9);
for I = 1:9
    specvol_CT_exact_error(:,I) = abs(specvol_CT_exact(:,I) - specvol_CT_exact(:,I+3));
end
gsw_cv.specvol_CT_exact_ca = nanmax(nanmax(specvol_CT_exact_error));
 
alpha_CT_exact = gsw_alpha_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
alpha_CT_exact_error = nan(45,9);
for I = 1:9
    alpha_CT_exact_error(:,I) = abs(alpha_CT_exact(:,I) - alpha_CT_exact(:,I+3));
end
gsw_cv.alpha_CT_exact_ca = nanmax(nanmax(alpha_CT_exact_error));

beta_CT_exact = gsw_beta_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
beta_CT_exact_error = nan(45,9);
for I = 1:9
    beta_CT_exact_error(:,I) = abs(beta_CT_exact(:,I) - beta_CT_exact(:,I+3));
end
gsw_cv.beta_CT_exact_ca = nanmax(nanmax(beta_CT_exact_error));

[v_vab_CT_exact, alpha_vab_CT_exact, beta_vab_CT_exact] = gsw_specvol_alpha_beta_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
v_vab_CT_exact_error = nan(45,9);
alpha_vab_CT_exact_error = nan(45,9);
beta_vab_CT_exact_error = nan(45,9);
for I = 1:9
    v_vab_CT_exact_error(:,I) = abs(v_vab_CT_exact(:,I) - v_vab_CT_exact(:,I+3));
end
gsw_cv.v_vab_CT_exact_ca = nanmax(nanmax(v_vab_CT_exact_error));
 for I = 1:9
    alpha_vab_CT_exact_error(:,I) = abs(alpha_vab_CT_exact(:,I) - alpha_vab_CT_exact(:,I+3));
end
gsw_cv.alpha_vab_CT_exact_ca = nanmax(nanmax(alpha_vab_CT_exact_error));
 for I = 1:9
    beta_vab_CT_exact_error(:,I) = abs(beta_vab_CT_exact(:,I) - beta_vab_CT_exact(:,I+3));
end
gsw_cv.beta_vab_CT_exact_ca = nanmax(nanmax(beta_vab_CT_exact_error));
  
alpha_on_beta_CT_exact = gsw_alpha_on_beta_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
alpha_on_beta_CT_exact_error = nan(45,9);
for I = 1:9
    alpha_on_beta_CT_exact_error(:,I) = abs(alpha_on_beta_CT_exact(:,I) - alpha_on_beta_CT_exact(:,I+3));
end
gsw_cv.alpha_on_beta_CT_exact_ca = nanmax(nanmax(alpha_on_beta_CT_exact_error));

[v_SA_CT_exact, v_CT_CT_exact, v_P_CT_exact] = gsw_specvol_first_derivatives_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
v_SA_CT_exact_error = nan(45,9);
v_CT_CT_exact_error = nan(45,9);
v_P_CT_exact_error = nan(45,9);
for I = 1:9
    v_SA_CT_exact_error(:,I) = abs(v_SA_CT_exact(:,I) - v_SA_CT_exact(:,I+3));
end
gsw_cv.v_SA_CT_exact_ca = nanmax(nanmax(v_SA_CT_exact_error));
 for I = 1:9
    v_CT_CT_exact_error(:,I) = abs(v_CT_CT_exact(:,I) - v_CT_CT_exact(:,I+3));
end
gsw_cv.v_CT_CT_exact_ca = nanmax(nanmax(v_CT_CT_exact_error));
 for I = 1:9
    v_P_CT_exact_error(:,I) = abs(v_P_CT_exact(:,I) - v_P_CT_exact(:,I+3));
end
gsw_cv.v_P_CT_exact_ca = nanmax(nanmax(v_P_CT_exact_error));

[v_SA_SA_CT_exact, v_SA_CT_CT_exact, v_CT_CT_CT_exact, v_SA_P_CT_exact, v_CT_P_CT_exact] = gsw_specvol_second_derivatives_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
v_SA_SA_CT_exact_error = nan(45,9);
v_SA_CT_CT_exact_error = nan(45,9);
v_CT_CT_CT_exact_error = nan(45,9);
v_SA_P_CT_exact_error = nan(45,9);
v_CT_P_CT_exact_error = nan(45,9);
for I = 1:9
    v_SA_SA_CT_exact_error(:,I) = abs(v_SA_SA_CT_exact(:,I) - v_SA_SA_CT_exact(:,I+3));
end
gsw_cv.v_SA_SA_CT_exact_ca = nanmax(nanmax(v_SA_SA_CT_exact_error));
for I = 1:9
    v_SA_CT_CT_exact_error(:,I) = abs(v_SA_CT_CT_exact(:,I) - v_SA_CT_CT_exact(:,I+3));
end
gsw_cv.v_SA_CT_CT_exact_ca = nanmax(nanmax(v_SA_CT_CT_exact_error));
 for I = 1:9
    v_CT_CT_CT_exact_error(:,I) = abs(v_CT_CT_CT_exact(:,I) - v_CT_CT_CT_exact(:,I+3));
end
gsw_cv.v_CT_CT_CT_exact_ca = nanmax(nanmax(v_CT_CT_error));
 for I = 1:9
    v_SA_P_CT_exact_error(:,I) = abs(v_SA_P_CT_exact(:,I) - v_SA_P_CT_exact(:,I+3));
end
gsw_cv.v_SA_P_CT_exact_ca = nanmax(nanmax(v_SA_P_CT_exact_error));
for I = 1:9
    v_CT_P_CT_exact_error(:,I) = abs(v_CT_P_CT_exact(:,I) - v_CT_P_CT_exact(:,I+3));
end
gsw_cv.v_CT_P_CT_exact_ca = nanmax(nanmax(v_CT_P_CT_exact_error));

[v_SA_wrt_h_CT_exact, v_h_CT_exact] = gsw_specvol_first_derivatives_wrt_enthalpy_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
v_SA_wrt_h_CT_exact_error = nan(45,9);
v_h_CT_exact_error = nan(45,9);
for I = 1:9
    v_SA_wrt_h_CT_exact_error(:,I) = abs(v_SA_wrt_h_CT_exact(:,I) - v_SA_wrt_h_CT_exact(:,I+3));
end
gsw_cv.v_SA_wrt_h_CT_exact_ca = nanmax(nanmax(v_SA_wrt_h_CT_exact_error));
 for I = 1:9
    v_h_CT_exact_error(:,I) = abs(v_h_CT_exact(:,I) - v_h_CT_exact(:,I+3));
end
gsw_cv.v_h_CT_exact_ca = nanmax(nanmax(v_h_CT_exact_error));

[v_SA_SA_wrt_h_CT_exact, v_SA_h_CT_exact, v_h_h_CT_exact] = gsw_specvol_second_derivatives_wrt_enthalpy_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
v_SA_SA_wrt_h_CT_exact_error = nan(45,9);
v_SA_h_CT_exact_error = nan(45,9);
v_h_h_CT_exact_error = nan(45,9);
for I = 1:9
    v_SA_SA_wrt_h_CT_exact_error(:,I) = abs(v_SA_SA_wrt_h_CT_exact(:,I) - v_SA_SA_wrt_h_CT_exact(:,I+3));
end
gsw_cv.v_SA_SA_wrt_h_CT_exact_ca = nanmax(nanmax(v_SA_SA_wrt_h_CT_exact_error));
for I = 1:9
    v_SA_h_CT_exact_error(:,I) = abs(v_SA_h_CT_exact(:,I) - v_SA_h_CT_exact(:,I+3));
end
gsw_cv.v_SA_h_CT_exact_ca = nanmax(nanmax(v_SA_h_CT_exact_error));
 for I = 1:9
    v_h_h_CT_exact_error(:,I) = abs(v_h_h_CT_exact(:,I) - v_h_h_CT_exact(:,I+3));
end
gsw_cv.v_h_h_CT_exact_ca = nanmax(nanmax(v_h_h_CT_exact_error));

specvol_anom_CT_exact = gsw_specvol_anom_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast,SA_ref,CT_ref);
specvol_anom_CT_exact_error = nan(45,9);
for I = 1:9
    specvol_anom_CT_exact_error(:,I) = abs(specvol_anom_CT_exact(:,I) - specvol_anom_CT_exact(:,I+3));
end
gsw_cv.specvol_anom_CT_exact_ca = nanmax(nanmax(specvol_anom_CT_exact_error));

specvol_anom_standard_CT_exact = gsw_specvol_anom_standard_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
specvol_anom_standard_CT_exact_error = nan(45,9);
for I = 1:9
    specvol_anom_standard_CT_exact_error(:,I) = abs(specvol_anom_CT_exact(:,I) - specvol_anom_CT_exact(:,I+3));
end
gsw_cv.specvol_anom_standard_CT_exact_ca = nanmax(nanmax(specvol_anom_standard_CT_exact_error));
  
rho_CT_exact = gsw_rho_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
rho_CT_exact_error = nan(45,9);
for I = 1:9
    rho_CT_exact_error(:,I) = abs(rho_CT_exact(:,I) - rho_CT_exact(:,I+3));
end
gsw_cv.rho_CT_exact_ca = nanmax(nanmax(rho_CT_exact_error));

[rho_rab_CT_exact, alpha_rab_CT_exact, beta_rab_CT_exact] = gsw_rho_alpha_beta_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
rho_rab_CT_exact_error = nan(45,9);
alpha_rab_CT_exact_error = nan(45,9);
beta_rab_CT_exact_error = nan(45,9);
for I = 1:9
    rho_rab_CT_exact_error(:,I) = abs(rho_rab_CT_exact(:,I) - rho_rab_CT_exact(:,I+3));
end
gsw_cv.rho_rab_CT_exact_ca = nanmax(nanmax(rho_rab_CT_exact_error));
 for I = 1:9
    alpha_rab_CT_exact_error(:,I) = abs(alpha_rab_CT_exact(:,I) - alpha_rab_CT_exact(:,I+3));
end
gsw_cv.alpha_rab_CT_exact_ca = nanmax(nanmax(alpha_rab_CT_exact_error));
 for I = 1:9
    beta_rab_CT_exact_error(:,I) = abs(beta_rab_CT_exact(:,I) - beta_rab_CT_exact(:,I+3));
end
gsw_cv.beta_rab_CT_exact_ca = nanmax(nanmax(beta_rab_CT_exact_error));

[rho_SA_CT_exact, rho_CT_CT_exact, rho_P_CT_exact] = gsw_rho_first_derivatives_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
rho_SA_CT_exact_error = nan(45,9);
rho_CT_CT_exact_error = nan(45,9);
rho_P_CT_exact_error = nan(45,9);
for I = 1:9
    rho_SA_CT_exact_error(:,I) = abs(rho_SA_CT_exact(:,I) - rho_SA_CT_exact(:,I+3));
end
gsw_cv.rho_SA_CT_exact_ca = nanmax(nanmax(rho_SA_CT_exact_error));
 for I = 1:9
    rho_CT_CT_exact_error(:,I) = abs(rho_CT_CT_exact(:,I) - rho_CT_CT_exact(:,I+3));
end
gsw_cv.rho_CT_CT_exact_ca = nanmax(nanmax(rho_CT_CT_exact_error));
 for I = 1:9
    rho_P_CT_exact_error(:,I) = abs(rho_P_CT_exact(:,I) - rho_P_CT_exact(:,I+3));
end
gsw_cv.rho_P_CT_exact_ca = nanmax(nanmax(rho_P_CT_exact_error));

[rho_SA_SA_CT_exact, rho_SA_CT_CT_exact, rho_CT_CT_CT_exact, rho_SA_P_CT_exact, rho_CT_P_CT_exact] = gsw_rho_second_derivatives_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
rho_SA_SA_CT_exact_error = nan(45,9);
rho_SA_CT_CT_exact_error = nan(45,9);
rho_CT_CT_CT_exact_error = nan(45,9);
rho_SA_P_CT_exact_error = nan(45,9);
rho_CT_P_CT_exact_error = nan(45,9);
for I = 1:9
    rho_SA_SA_CT_exact_error(:,I) = abs(rho_SA_SA_CT_exact(:,I) - rho_SA_SA_CT_exact(:,I+3));
end
gsw_cv.rho_SA_SA_CT_exact_ca = nanmax(nanmax(rho_SA_SA_CT_exact_error));
for I = 1:9
    rho_SA_CT_CT_exact_error(:,I) = abs(rho_SA_CT_CT_exact(:,I) - rho_SA_CT_CT_exact(:,I+3));
end
gsw_cv.rho_SA_CT_CT_exact_ca = nanmax(nanmax(rho_SA_CT_CT_exact_error));
 for I = 1:9
    rho_CT_CT_CT_exact_error(:,I) = abs(rho_CT_CT_CT_exact(:,I) - rho_CT_CT_CT_exact(:,I+3));
end
gsw_cv.rho_CT_CT_CT_exact_ca = nanmax(nanmax(rho_CT_CT_CT_exact_error));
 for I = 1:9
    rho_SA_P_CT_exact_error(:,I) = abs(rho_SA_P_CT_exact(:,I) - rho_SA_P_CT_exact(:,I+3));
end
gsw_cv.rho_SA_P_CT_exact_ca = nanmax(nanmax(rho_SA_P_CT_exact_error));
for I = 1:9
    rho_CT_P_CT_exact_error(:,I) = abs(rho_CT_P_CT_exact(:,I) - rho_CT_P_CT_exact(:,I+3));
end
gsw_cv.rho_CT_P_CT_exact_ca = nanmax(nanmax(rho_CT_P_CT_exact_error));

[rho_SA_wrt_h_CT_exact, rho_h_CT_exact] = gsw_rho_first_derivatives_wrt_enthalpy_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
rho_SA_wrt_h_CT_exact_error = nan(45,9);
rho_h_CT_exact_error = nan(45,9);
for I = 1:9
    rho_SA_wrt_h_CT_exact_error(:,I) = abs(rho_SA_wrt_h_CT_exact(:,I) - rho_SA_wrt_h_CT_exact(:,I+3));
end
gsw_cv.rho_SA_wrt_h_CT_exact_ca = nanmax(nanmax(rho_SA_wrt_h_CT_exact_error));
 for I = 1:9
    rho_h_CT_exact_error(:,I) = abs(rho_h_CT_exact(:,I) - rho_h_CT_exact(:,I+3));
end
gsw_cv.rho_h_CT_exact_ca = nanmax(nanmax(rho_h_CT_exact_error));

[rho_SA_SA_wrt_h_CT_exact, rho_SA_h_CT_exact, rho_h_h_CT_exact] = gsw_rho_second_derivatives_wrt_enthalpy_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
rho_SA_SA_wrt_h_CT_exact_error = nan(45,9);
rho_SA_h_CT_exact_error = nan(45,9);
rho_h_h_CT_exact_error = nan(45,9);
for I = 1:9
    rho_SA_SA_wrt_h_CT_exact_error(:,I) = abs(rho_SA_SA_wrt_h_CT_exact(:,I) - rho_SA_SA_wrt_h_CT_exact(:,I+3));
end
gsw_cv.rho_SA_SA_wrt_h_CT_exact_ca = nanmax(nanmax(rho_SA_SA_wrt_h_CT_exact_error));
for I = 1:9
    rho_SA_h_CT_exact_error(:,I) = abs(rho_SA_h_CT_exact(:,I) - rho_SA_h_CT_exact(:,I+3));
end
gsw_cv.rho_SA_h_CT_exact_ca = nanmax(nanmax(rho_SA_h_CT_exact_error));
 for I = 1:9
    rho_h_h_CT_exact_error(:,I) = abs(rho_h_h_CT_exact(:,I) - rho_h_h_CT_exact(:,I+3));
end
gsw_cv.rho_h_h_CT_exact_ca = nanmax(nanmax(rho_h_h_CT_exact_error));

sigma0_CT_exact = gsw_sigma0_CT_exact(SA_chck_cast,CT_chck_cast);
sigma0_CT_exact_error = nan(45,9);
for I = 1:9
    sigma0_CT_exact_error(:,I) = abs(sigma0_CT_exact(:,I) - sigma0_CT_exact(:,I+3));
end
gsw_cv.sigma0_CT_exact_ca = nanmax(nanmax(sigma0_CT_exact_error));
  
sigma1_CT_exact = gsw_sigma1_CT_exact(SA_chck_cast,CT_chck_cast);
sigma1_CT_exact_error = nan(45,9);
for I = 1:9
    sigma1_CT_exact_error(:,I) = abs(sigma1_CT_exact(:,I) - sigma1_CT_exact(:,I+3));
end
gsw_cv.sigma1_CT_exact_ca = nanmax(nanmax(sigma1_CT_exact_error));
 
sigma2_CT_exact = gsw_sigma2_CT_exact(SA_chck_cast,CT_chck_cast);
sigma2_CT_exact_error = nan(45,9);
for I = 1:9
    sigma2_CT_exact_error(:,I) = abs(sigma2_CT_exact(:,I) - sigma2_CT_exact(:,I+3));
end
gsw_cv.sigma2_CT_exact_ca = nanmax(nanmax(sigma2_CT_exact_error));
 
sigma3_CT_exact = gsw_sigma3_CT_exact(SA_chck_cast,CT_chck_cast);
sigma3_CT_exact_error = nan(45,9);
for I = 1:9
    sigma3_CT_exact_error(:,I) = abs(sigma3_CT_exact(:,I) - sigma3_CT_exact(:,I+3));
end
gsw_cv.sigma3_CT_exact_ca = nanmax(nanmax(sigma3_CT_exact_error));
 
sigma4_CT_exact = gsw_sigma4_CT_exact(SA_chck_cast,CT_chck_cast);
sigma4_CT_exact_error = nan(45,9);
for I = 1:9
    sigma4_CT_exact_error(:,I) = abs(sigma4_CT_exact(:,I) - sigma4_CT_exact(:,I+3));
end
gsw_cv.sigma4_CT_exact_ca = nanmax(nanmax(sigma4_CT_exact_error));

cabbeling_CT_exact = gsw_cabbeling_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
cabbeling_CT_exact_error = nan(45,9);
for I = 1:9
    cabbeling_CT_exact_error(:,I) = abs(cabbeling_CT_exact(:,I) - cabbeling_CT_exact(:,I+3));
end
gsw_cv.cabbeling_CT_exact_ca = nanmax(nanmax(cabbeling_CT_exact_error));
  
thermobaric_CT_exact = gsw_thermobaric_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
thermobaric_CT_exact_error = nan(45,9);
for I = 1:9
    thermobaric_CT_exact_error(:,I) = abs(thermobaric_CT_exact(:,I) - thermobaric_CT_exact(:,I+3));
end
gsw_cv.thermobaric_CT_exact_ca = nanmax(nanmax(thermobaric_CT_exact_error));

enthalpy_CT_exact =  gsw_enthalpy_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
enthalpy_CT_exact_error = nan(45,9);
for I = 1:9
    enthalpy_CT_exact_error(:,I) = abs(enthalpy_CT_exact(:,I) - enthalpy_CT_exact(:,I+3));
end
gsw_cv.enthalpy_CT_exact_ca = nanmax(nanmax(enthalpy_CT_exact_error));

enthalpy_diff_CT_exact =  gsw_enthalpy_diff_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast_shallow,p_chck_cast_deep);
enthalpy_diff_CT_exact_error = nan(45,9);
for I = 1:9
    enthalpy_diff_CT_exact_error(:,I) = abs(enthalpy_diff_CT_exact(:,I) - enthalpy_diff_CT_exact(:,I+3));
end
gsw_cv.enthalpy_diff_CT_exact_ca = nanmax(nanmax(enthalpy_diff_CT_exact_error));
  
dynamic_enthalpy_CT_exact =  gsw_dynamic_enthalpy_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
dynamic_enthalpy_CT_exact_error = nan(45,9);
for I = 1:9
    dynamic_enthalpy_CT_exact_error(:,I) = abs(dynamic_enthalpy_CT_exact(:,I) - dynamic_enthalpy_CT_exact(:,I+3));
end
gsw_cv.dynamic_enthalpy_CT_exact_ca = nanmax(nanmax(dynamic_enthalpy_CT_exact_error));

[h_SA_CT_exact, h_CT_CT_exact] = gsw_enthalpy_first_derivatives_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
h_SA_CT_exact_error = nan(45,9);
h_CT_CT_exact_error = nan(45,9);
for I = 1:9
    h_SA_CT_exact_error(:,I) = abs(h_SA_CT_exact(:,I) - h_SA_CT_exact(:,I+3));
end
gsw_cv.h_SA_CT_exact_ca = nanmax(nanmax(h_SA_CT_exact_error));
 for I = 1:9
    h_CT_CT_exact_error(:,I) = abs(h_CT_CT_exact(:,I) - h_CT_CT_exact(:,I+3));
end
gsw_cv.h_CT_CT_exact_ca = nanmax(nanmax(h_CT_CT_exact_error));
  
[h_SA_SA_CT_exact, h_SA_CT_CT_exact, h_CT_CT_CT_exact] = gsw_enthalpy_second_derivatives_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
h_SA_SA_CT_exact_error = nan(45,9);
h_SA_CT_CT_exact_error = nan(45,9);
h_CT_CT_CT_exact_error = nan(45,9);
for I = 1:9
    h_SA_SA_CT_exact_error(:,I) = abs(h_SA_SA_CT_exact(:,I) - h_SA_SA_CT_exact(:,I+3));
end
gsw_cv.h_SA_SA_CT_exact_ca = nanmax(nanmax(h_SA_SA_CT_exact_error));
 for I = 1:9
    h_SA_CT_CT_exact_error(:,I) = abs(h_SA_CT_CT_exact(:,I) - h_SA_CT_CT_exact(:,I+3));
end
gsw_cv.h_SA_CT_CT_exact_ca = nanmax(nanmax(h_SA_CT_CT_exact_error));
 for I = 1:9
    h_CT_CT_CT_exact_error(:,I) = abs(h_CT_CT_CT_exact(:,I) - h_CT_CT_CT_exact(:,I+3));
end
gsw_cv.h_CT_CT_CT_exact_ca = nanmax(nanmax(h_CT_CT_CT_exact_error));

sound_speed_CT_exact = gsw_sound_speed_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
sound_speed_CT_exact_error = nan(45,9);
for I = 1:9
    sound_speed_CT_exact_error(:,I) = abs(sound_speed_CT_exact(:,I) - sound_speed_CT_exact(:,I+3));
end
gsw_cv.sound_speed_CT_exact_ca = nanmax(nanmax(sound_speed_CT_exact_error));

kappa_CT_exact = gsw_kappa_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
kappa_CT_exact_error = nan(45,9);
for I = 1:9
    kappa_CT_exact_error(:,I) = abs(kappa_CT_exact(:,I) - kappa_CT_exact(:,I+3));
end
gsw_cv.kappa_CT_exact_ca = nanmax(nanmax(kappa_CT_exact_error));

internal_energy_CT_exact = gsw_internal_energy_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
internal_energy_CT_exact_error = nan(45,9);
for I = 1:9
    internal_energy_CT_exact_error(:,I) = abs(internal_energy_CT_exact(:,I) - internal_energy_CT_exact(:,I+3));
end
gsw_cv.internal_energy_CT_exact_ca = nanmax(nanmax(internal_energy_CT_exact_error));

[u_SA_CT_exact, u_CT_CT_exact, u_P_CT_exact] = gsw_internal_energy_first_derivatives_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
u_SA_CT_exact_error = nan(45,9);
u_CT_CT_exact_error = nan(45,9);
u_P_CT_exact_error = nan(45,9);
for I = 1:9
    u_SA_CT_exact_error(:,I) = abs(u_SA_CT_exact(:,I) - u_SA_CT_exact(:,I+3));
end
gsw_cv.u_SA_CT_exact_ca = nanmax(nanmax(u_SA_CT_exact_error));
 for I = 1:9
    u_CT_CT_exact_error(:,I) = abs(u_CT_CT_exact(:,I) - u_CT_CT_exact(:,I+3));
end
gsw_cv.u_CT_CT_exact_ca = nanmax(nanmax(u_CT_CT_exact_error));
 for I = 1:9
    u_P_CT_exact_error(:,I) = abs(u_P_CT_exact(:,I) - u_P_CT_exact(:,I+3));
end
gsw_cv.u_P_CT_exact_ca = nanmax(nanmax(u_P_CT_exact_error));

[u_SA_SA_CT_exact, u_SA_CT_CT_exact, u_CT_CT_CT_exact, u_SA_P_CT_exact, u_CT_P_CT_exact] = gsw_internal_energy_second_derivatives_CT_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
u_SA_SA_CT_exact_error = nan(45,9);
u_SA_CT_CT_exact_error = nan(45,9);
u_CT_CT_CT_exact_error = nan(45,9);
u_SA_P_CT_exact_error = nan(45,9);
u_CT_P_CT_exact_error = nan(45,9);
for I = 1:9
    u_SA_SA_CT_exact_error(:,I) = abs(u_SA_SA_CT_exact(:,I) - u_SA_SA_CT_exact(:,I+3));
end
gsw_cv.u_SA_SA_CT_exact_ca = nanmax(nanmax(u_SA_SA_CT_exact_error));
for I = 1:9
    u_SA_CT_CT_exact_error(:,I) = abs(u_SA_CT_CT_exact(:,I) - u_SA_CT_CT_exact(:,I+3));
end
gsw_cv.u_SA_CT_CT_exact_ca = nanmax(nanmax(u_SA_CT_CT_exact_error));
 for I = 1:9
    u_CT_CT_CT_exact_error(:,I) = abs(u_CT_CT_CT_exact(:,I) - u_CT_CT_CT_exact(:,I+3));
end
gsw_cv.u_CT_CT_CT_exact_ca = nanmax(nanmax(u_CT_CT_CT_exact_error));
 for I = 1:9
    u_SA_P_CT_exact_error(:,I) = abs(u_SA_P_CT_exact(:,I) - u_SA_P_CT_exact(:,I+3));
end
gsw_cv.u_SA_P_CT_exact_ca = nanmax(nanmax(u_SA_P_CT_exact_error));
for I = 1:9
    u_CT_P_CT_exact_error(:,I) = abs(u_CT_P_CT_exact(:,I) - u_CT_P_CT_exact(:,I+3));
end
gsw_cv.u_CT_P_CT_exact_ca = nanmax(nanmax(u_CT_P_CT_exact_error));

CT_from_enthalpy_exact = gsw_CT_from_enthalpy_exact(SA_chck_cast,enthalpy_CT_exact,p_chck_cast);
CT_from_enthalpy_exact_error = nan(45,9);
for I = 1:9
    CT_from_enthalpy_exact_error(:,I) = abs(CT_from_enthalpy_exact(:,I) - CT_from_enthalpy_exact(:,I+3));
end
gsw_cv.CT_from_enthalpy_exact_ca = nanmax(nanmax(CT_from_enthalpy_exact_error));

SA_from_rho_CT_exact = gsw_SA_from_rho_CT_exact(rho_CT_exact,CT_chck_cast,p_chck_cast);
SA_from_rho_CT_exact_error = nan(45,9);
for I = 1:9
    SA_from_rho_CT_exact_error(:,I) = abs(SA_from_rho_CT_exact(:,I) - SA_from_rho_CT_exact(:,I+3));
end
gsw_cv.SA_from_rho_CT_exact_ca = nanmax(nanmax(SA_from_rho_CT_exact_error));

CT_from_rho_exact = gsw_CT_from_rho_exact(rho_CT_exact,SA_chck_cast,p_chck_cast);
CT_from_rho_exact_error = nan(45,9);
for I = 1:9
    CT_from_rho_exact_error(:,I) = abs(CT_from_rho_exact(:,I) - CT_from_rho_exact(:,I+3));
end
gsw_cv.CT_from_rho_exact_ca = nanmax(nanmax(CT_from_rho_exact_error));

CT_maxdensity_exact = gsw_CT_maxdensity_exact(SA_chck_cast,p_chck_cast);
CT_maxdensity_exact_error = nan(45,9);
 for I = 1:9
    CT_maxdensity_exact_error(:,I) = abs(CT_maxdensity_exact(:,I) - CT_maxdensity_exact(:,I+3));
end
gsw_cv.CT_maxdensity_exact_ca = nanmax(nanmax(CT_maxdensity_exact_error));

%% Labrortory functions

rho_t_exact = gsw_rho_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
rho_t_exact_error = nan(45,9);
for I = 1:9
    rho_t_exact_error(:,I) = abs(rho_t_exact(:,I) - rho_t_exact(:,I+3));
end
gsw_cv.rho_t_exact_ca = nanmax(nanmax(rho_t_exact_error));

SA_from_rho_t_exact = gsw_SA_from_rho_t_exact(rho_t_exact,t_chck_cast,p_chck_cast);
SA_from_rho_t_exact_error = nan(45,9);
for I = 1:9
    SA_from_rho_t_exact_error(:,I) = abs(SA_from_rho_t_exact(:,I) - SA_from_rho_t_exact(:,I+3));
end
gsw_cv.SA_from_rho_t_exact_ca = nanmax(nanmax(SA_from_rho_t_exact_error));

deltaSA_from_rho_t_exact = gsw_deltaSA_from_rho_t_exact(rho_t_exact,SP_chck_cast,t_chck_cast,p_chck_cast);
deltaSA_from_rho_t_exact_error = nan(45,9);
for I = 1:9
    deltaSA_from_rho_t_exact_error(:,I) = abs(deltaSA_from_rho_t_exact(:,I) - deltaSA_from_rho_t_exact(:,I+3));
end
gsw_cv.deltaSA_from_rho_t_exact_ca = nanmax(nanmax(deltaSA_from_rho_t_exact_error));

%% basic thermodynamic properties interms of in-situ t, derived from the exact Gibbs function

% rho_t_exact = gsw_rho_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
% rho_t_exact_error = nan(45,9);
% for I = 1:9
%     rho_t_exact_error(:,I) = abs(rho_t_exact(:,I) - rho_t_exact(:,I+3));
% end
% gsw_cv.rho_t_exact_ca = nanmax(nanmax(rho_t_exact_error));

pot_rho_t_exact = gsw_pot_rho_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast,pr);
pot_rho_t_exact_error = nan(45,9);
for I = 1:9
    pot_rho_t_exact_error(:,I) = abs(pot_rho_t_exact(:,I) - pot_rho_t_exact(:,I+3));
end
gsw_cv.pot_rho_t_exact_ca = nanmax(nanmax(pot_rho_t_exact_error));
  
sigma0_pt0_exact = gsw_sigma0_pt0_exact(SA_chck_cast,pt0);
sigma0_pt0_exact_error = nan(45,9);
for I = 1:9
    sigma0_pt0_exact_error(:,I) = abs(sigma0_pt0_exact(:,I) - sigma0_pt0_exact(:,I+3));
end
gsw_cv.sigma0_pt0_exact_ca = nanmax(nanmax(sigma0_pt0_exact_error));

alpha_wrt_CT_t_exact = gsw_alpha_wrt_CT_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
alpha_wrt_CT_t_exact_error = nan(45,9);
for I = 1:9
    alpha_wrt_CT_t_exact_error(:,I) = abs(alpha_wrt_CT_t_exact(:,I) - alpha_wrt_CT_t_exact(:,I+3));
end
gsw_cv.alpha_wrt_CT_t_exact_ca = nanmax(nanmax(alpha_wrt_CT_t_exact_error));
  
alpha_wrt_pt_t_exact = gsw_alpha_wrt_pt_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
alpha_wrt_pt_t_exact_error = nan(45,9);
for I = 1:9
    alpha_wrt_pt_t_exact_error(:,I) = abs(alpha_wrt_pt_t_exact(:,I) - alpha_wrt_pt_t_exact(:,I+3));
end
gsw_cv.alpha_wrt_pt_t_exact_ca = nanmax(nanmax(alpha_wrt_pt_t_exact_error));
  
alpha_wrt_t_exact = gsw_alpha_wrt_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
alpha_wrt_t_exact_error = nan(45,9);
for I = 1:9
    alpha_wrt_t_exact_error(:,I) = abs(alpha_wrt_t_exact(:,I) - alpha_wrt_t_exact(:,I+3));
end
gsw_cv.alpha_wrt_t_exact_ca = nanmax(nanmax(alpha_wrt_t_exact_error));
  
beta_const_CT_t_exact = gsw_beta_const_CT_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
beta_const_CT_t_exact_error = nan(45,9);
for I = 1:9
    beta_const_CT_t_exact_error(:,I) = abs(beta_const_CT_t_exact(:,I) - beta_const_CT_t_exact(:,I+3));
end
gsw_cv.beta_const_CT_t_exact_ca = nanmax(nanmax(beta_const_CT_t_exact_error));
  
beta_const_pt_t_exact = gsw_beta_const_pt_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
beta_const_pt_t_exact_error = nan(45,9);
for I = 1:9
    beta_const_pt_t_exact_error(:,I) = abs(beta_const_pt_t_exact(:,I) - beta_const_pt_t_exact(:,I+3));
end
gsw_cv.beta_const_pt_t_exact_ca = nanmax(nanmax(beta_const_pt_t_exact_error));
  
beta_const_t_exact = gsw_beta_const_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
beta_const_t_exact_error = nan(45,9);
for I = 1:9
    beta_const_t_exact_error(:,I) = abs(beta_const_t_exact(:,I) - beta_const_t_exact(:,I+3));
end
gsw_cv.beta_const_t_exact_ca = nanmax(nanmax(beta_const_t_exact_error));
  
specvol_t_exact = gsw_specvol_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast); 
specvol_t_exact_error = nan(45,9);
for I = 1:9
    specvol_t_exact_error(:,I) = abs(specvol_t_exact(:,I) - specvol_t_exact(:,I+3));
end
gsw_cv.specvol_t_exact_ca = nanmax(nanmax(specvol_t_exact_error));
 
specvol_anom_standard_t_exact = gsw_specvol_anom_standard_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
specvol_anom_standard_t_exact_error = nan(45,9);
for I = 1:9
    specvol_anom_standard_t_exact_error(:,I) = abs(specvol_anom_standard_t_exact(:,I) - specvol_anom_standard_t_exact(:,I+3));
end
gsw_cv.specvol_anom_standard_t_exact_ca = nanmax(nanmax(specvol_anom_standard_t_exact_error));

sound_speed_t_exact = gsw_sound_speed_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
sound_speed_t_exact_error = nan(45,9);
for I = 1:9
    sound_speed_t_exact_error(:,I) = abs(sound_speed_t_exact(:,I) - sound_speed_t_exact(:,I+3));
end
gsw_cv.sound_speed_t_exact_ca = nanmax(nanmax(sound_speed_t_exact_error));
  
kappa_t_exact = gsw_kappa_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
kappa_t_exact_error = nan(45,9);
for I = 1:9
    kappa_t_exact_error(:,I) = abs(kappa_t_exact(:,I) - kappa_t_exact(:,I+3));
end
gsw_cv.kappa_t_exact_ca = nanmax(nanmax(kappa_t_exact_error));
 
kappa_const_t_exact = gsw_kappa_const_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
kappa_const_t_exact_error = nan(45,9);
for I = 1:9
    kappa_const_t_exact_error(:,I) = abs(kappa_const_t_exact(:,I) - kappa_const_t_exact(:,I+3));
end
gsw_cv.kappa_const_t_exact_ca = nanmax(nanmax(kappa_const_t_exact_error));

% SA_from_rho_t_exact = gsw_SA_from_rho_t_exact(rho_t_exact,t_chck_cast,p_chck_cast);
% SA_from_rho_t_exact_error = nan(45,9);
% for I = 1:9
%     SA_from_rho_t_exact_error(:,I) = abs(SA_from_rho_t_exact(:,I) - SA_from_rho_t_exact(:,I+3));
% end
% gsw_cv.SA_from_rho_t_exact_ca = nanmax(nanmax(SA_from_rho_t_exact_error));
  
t_from_rho_exact = gsw_t_from_rho_exact(rho_t_exact,SA_chck_cast,p_chck_cast);
t_from_rho_exact_error = nan(45,9);
for I = 1:9
    t_from_rho_exact_error(:,I) = abs(t_from_rho_exact(:,I) - t_from_rho_exact(:,I+3));
end
gsw_cv.t_from_rho_exact_ca = nanmax(nanmax(t_from_rho_exact_error));

t_maxdensity_exact = gsw_t_maxdensity_exact(SA_chck_cast,p_chck_cast);
t_maxdensity_exact_error = nan(45,9);
for I = 1:9
    t_maxdensity_exact_error(:,I) = abs(t_maxdensity_exact(:,I) - t_maxdensity_exact(:,I+3));
end
gsw_cv.t_maxdensity_exact_ca = nanmax(nanmax(t_maxdensity_exact_error));

internal_energy_t_exact = gsw_internal_energy_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
internal_energy_t_exact_error = nan(45,9);
for I = 1:9
    internal_energy_t_exact_error(:,I) = abs(internal_energy_t_exact(:,I) - internal_energy_t_exact(:,I+3));
end
gsw_cv.internal_energy_t_exact_ca = nanmax(nanmax(internal_energy_t_exact_error));

enthalpy_t_exact = gsw_enthalpy_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
enthalpy_t_exact_error = nan(45,9);
for I = 1:9
    enthalpy_t_exact_error(:,I) = abs(enthalpy_t_exact(:,I) - enthalpy_t_exact(:,I+3));
end
gsw_cv.enthalpy_t_exact_ca = nanmax(nanmax(enthalpy_t_exact_error));

dynamic_enthalpy_t_exact = gsw_dynamic_enthalpy_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
dynamic_enthalpy_t_exact_error = nan(45,9);
for I = 1:9
    dynamic_enthalpy_t_exact_error(:,I) = abs(dynamic_enthalpy_t_exact(:,I) - dynamic_enthalpy_t_exact(:,I+3));
end
gsw_cv.dynamic_enthalpy_t_exact_ca = nanmax(nanmax(dynamic_enthalpy_t_exact_error));

[CT_SA_wrt_t, CT_T_wrt_t, CT_P_wrt_t] = gsw_CT_first_derivatives_wrt_t_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
CT_SA_wrt_t_error = nan(45,9);
CT_T_wrt_t_error = nan(45,9);
CT_P_wrt_t_error = nan(45,9);
for I = 1:9
    CT_SA_wrt_t_error(:,I) = abs(CT_SA_wrt_t(:,I) - CT_SA_wrt_t(:,I+3));
end
gsw_cv.CT_SA_wrt_t_ca = nanmax(nanmax(CT_SA_wrt_t_error));
 for I = 1:9
    CT_T_wrt_t_error(:,I) = abs(CT_T_wrt_t(:,I) - CT_T_wrt_t(:,I+3));
end
gsw_cv.CT_T_wrt_t_ca = nanmax(nanmax(CT_T_wrt_t_error));
 for I = 1:9
    CT_P_wrt_t_error(:,I) = abs(CT_P_wrt_t(:,I) - CT_P_wrt_t(:,I+3));
end
gsw_cv.CT_P_wrt_t_ca = nanmax(nanmax(CT_P_wrt_t_error));

[h_SA_wrt_t, h_T_wrt_t, h_P_wrt_t] = gsw_enthalpy_first_derivatives_wrt_t_exact(SA_chck_cast,CT_chck_cast,p_chck_cast);
h_SA_wrt_t_error = nan(45,9);
h_T_wrt_t_error = nan(45,9);
h_P_wrt_t_error = nan(45,9);
for I = 1:9
    h_SA_wrt_t_error(:,I) = abs(h_SA_wrt_t(:,I) - h_SA_wrt_t(:,I+3));
end
gsw_cv.h_SA_wrt_t_ca = nanmax(nanmax(h_SA_wrt_t_error));
 for I = 1:9
    h_T_wrt_t_error(:,I) = abs(h_T_wrt_t(:,I) - h_T_wrt_t(:,I+3));
end
gsw_cv.h_T_wrt_t_ca = nanmax(nanmax(h_T_wrt_t_error));
 for I = 1:9
    h_P_wrt_t_error(:,I) = abs(h_P_wrt_t(:,I) - h_P_wrt_t(:,I+3));
end
gsw_cv.h_P_wrt_t_ca = nanmax(nanmax(h_P_wrt_t_error));

cp_t_exact = gsw_cp_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
cp_t_exact_error = nan(45,9);
for I = 1:9
    cp_t_exact_error(:,I) = abs(cp_t_exact(:,I) - cp_t_exact(:,I+3));
end
gsw_cv.cp_t_exact_ca = nanmax(nanmax(cp_t_exact_error));

isochoric_heat_cap_t_exact = gsw_isochoric_heat_cap_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
isochoric_heat_cap_t_exact_error = nan(45,9);
for I = 1:9
    isochoric_heat_cap_t_exact_error(:,I) = abs(isochoric_heat_cap_t_exact(:,I) - isochoric_heat_cap_t_exact(:,I+3));
end
gsw_cv.isochoric_heat_cap_t_exact_ca = nanmax(nanmax(isochoric_heat_cap_t_exact_error));

chem_potential_relative_t_exact =  gsw_chem_potential_relative_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
chem_potential_relative_t_exact_error = nan(45,9);
for I = 1:9
    chem_potential_relative_t_exact_error(:,I) = abs(chem_potential_relative_t_exact(:,I) - chem_potential_relative_t_exact(:,I+3));
end
gsw_cv.chem_potential_relative_t_exact_ca = nanmax(nanmax(chem_potential_relative_t_exact_error));
  
chem_potential_water_t_exact =  gsw_chem_potential_water_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
chem_potential_water_t_exact_error = nan(45,9);
for I = 1:9
    chem_potential_water_t_exact_error(:,I) = abs(chem_potential_water_t_exact(:,I) - chem_potential_water_t_exact(:,I+3));
end
gsw_cv.chem_potential_water_t_exact_ca = nanmax(nanmax(chem_potential_water_t_exact_error));
  
chem_potential_salt_t_exact =  gsw_chem_potential_salt_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
chem_potential_salt_t_exact_error = nan(45,9);
for I = 1:9
    chem_potential_salt_t_exact_error(:,I) = abs(chem_potential_salt_t_exact(:,I) - chem_potential_salt_t_exact(:,I+3));
end
gsw_cv.chem_potential_salt_t_exact_ca = nanmax(nanmax(chem_potential_salt_t_exact_error));
  
t_deriv_chem_potential_water_t_exact  =  gsw_t_deriv_chem_potential_water_t_exact (SA_chck_cast,t_chck_cast,p_chck_cast);
t_deriv_chem_potential_water_t_exact_error = nan(45,9);
for I = 1:9
    t_deriv_chem_potential_water_t_exact_error(:,I) = abs(t_deriv_chem_potential_water_t_exact(:,I) - t_deriv_chem_potential_water_t_exact(:,I+3));
end
gsw_cv.t_deriv_chem_potential_water_t_exact_ca = nanmax(nanmax(t_deriv_chem_potential_water_t_exact_error));

dilution_coefficient_t_exact = gsw_dilution_coefficient_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
dilution_coefficient_t_exact_error = nan(45,9);
for I = 1:9
    dilution_coefficient_t_exact_error(:,I) = abs(dilution_coefficient_t_exact(:,I) - dilution_coefficient_t_exact(:,I+3));
end
gsw_cv.dilution_coefficient_t_exact_ca = nanmax(nanmax(dilution_coefficient_t_exact_error));

Gibbs_energy_t_exact = gsw_Gibbs_energy_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
Gibbs_energy_t_exact_error = nan(45,9);
for I = 1:9
    Gibbs_energy_t_exact_error(:,I) = abs(Gibbs_energy_t_exact(:,I) - Gibbs_energy_t_exact(:,I+3));
end
gsw_cv.Gibbs_energy_t_exact_ca = nanmax(nanmax(Gibbs_energy_t_exact_error));

Helmholtz_energy_t_exact = gsw_Helmholtz_energy_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
Helmholtz_energy_t_exact_error = nan(45,9);
for I = 1:9
    Helmholtz_energy_t_exact_error(:,I) = abs(Helmholtz_energy_t_exact(:,I) - Helmholtz_energy_t_exact(:,I+3));
end
gsw_cv.Helmholtz_energy_t_exact_ca = nanmax(nanmax(Helmholtz_energy_t_exact_error));
  
osmotic_coefficient_t_exact = gsw_osmotic_coefficient_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
osmotic_coefficient_t_exact_error = nan(45,9);
for I = 1:9
    osmotic_coefficient_t_exact_error(:,I) = abs(osmotic_coefficient_t_exact(:,I) - osmotic_coefficient_t_exact(:,I+3));
end
gsw_cv.osmotic_coefficient_t_exact_ca = nanmax(nanmax(osmotic_coefficient_t_exact_error));
  
osmotic_pressure_t_exact = gsw_osmotic_pressure_t_exact(SA_chck_cast,t_chck_cast,p_chck_cast);
osmotic_pressure_t_exact_error = nan(45,9);
for I = 1:9
    osmotic_pressure_t_exact_error(:,I) = abs(osmotic_pressure_t_exact(:,I) - osmotic_pressure_t_exact(:,I+3));
end
gsw_cv.osmotic_pressure_t_exact_ca = nanmax(nanmax(osmotic_pressure_t_exact_error));
 


%% Library

Fdelta = gsw_Fdelta(p_chck_cast,long_chck_cast,lat_chck_cast);
Fdelta_error = nan(45,9);
for I = 1:9
    Fdelta_error(:,I) = abs(Fdelta(:,I) - Fdelta(:,I+3));
end
gsw_cv.Fdelta_ca = nanmax(nanmax(Fdelta_error));

lat_chck_cast_temp = nan(45,length(long_chck_cast(1,:)));
long_chck_cast_temp = lat_chck_cast_temp;

for I = 1:45
    long_chck_cast_temp(I,:) = long_chck_cast(1,:);
    lat_chck_cast_temp(I,:) = lat_chck_cast(1,:);
end
[I] = find(~isnan(p_chck_cast));
deltaSA_atlas = nan(45,12);
deltaSA_atlas(I) = gsw_deltaSA_atlas(p_chck_cast(I),long_chck_cast_temp(I),lat_chck_cast_temp(I));
deltaSA_atlas_error = nan(45,9);
for I = 1:9
    deltaSA_atlas_error(:,I) = abs(deltaSA_atlas(:,I) - deltaSA_atlas(:,I+3));
end
gsw_cv.deltaSA_atlas_ca = nanmax(nanmax(deltaSA_atlas_error));



%save library\gsw_chck_vals_errors_v3_0.mat gsw_cv
save gsw_chck_vals_errors_v3_06_12.mat gsw_cv