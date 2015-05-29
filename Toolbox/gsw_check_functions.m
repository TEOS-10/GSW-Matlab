if isempty(which('gsw_gibbs.html'))
    try
        cprintf('err','You need to add the GSW "html" subdirectory to your path.\n')
    catch
        fprintf(2,'You need to add the GSW "html" subdirectory to your path. \n');
    end
end
if isempty(which('gsw_gibbs.m'))
    try
        cprintf('err','You need to add the GSW "library" subdirectory to your path.\n')
    catch
        fprintf(2,'You need to add the GSW "library" subdirectory to your path. \n');
    end
end
if isempty(which('gibbs.pdf'))
    try
        cprintf('err','You need to add the GSW "pdf" subdirectory to your path.\n')
    catch
        fprintf(2,'You need to add the GSW "pdf" subdirectory to your path. \n');
    end
end
if isempty(which('gsw_gibbs.html')) | isempty(which('gsw_gibbs.m'))| isempty(which('gibbs.pdf'))
    error('You have not added the GSW subdirectories to you MATLAB Path')
end

gsw_data = 'gsw_data_v2_0.mat';
gsw_data_file = which(gsw_data);
load (gsw_data_file,'gsw_cv');

try
    cprintf('text',' \n');
    cprintf('strings','This function is running three stored vertical profiles through\n');
    cprintf('strings','all the functions in the GSW Oceanographic Toolbox, and then checks\n');
    cprintf('strings','that the outputs are all within pre-defined limits of the correct\n');
    cprintf('strings','answers.  These pre-defined limits are a factor of approximately\n');
    cprintf('strings','a hundred larger than the errors expected from numerical round-off.\n');
    cprintf('text',' \n');
    cprintf('text',' checking ');
catch
    fprintf(1,' \n');
    fprintf(1,'This function is running three stored vertical profiles through\n');
    fprintf(1,'all the functions in the GSW Oceanographic Toolbox, and then checks\n');
    fprintf(1,'that the outputs are all within pre-defined limits of the correct\n');
    fprintf(1,'answers.  These pre-defined limits are a factor of approximately\n');
    fprintf(1,'a hundred larger than the errors expected from numerical round-off.\n');
    fprintf(1,' \n');
    fprintf(1,' checking ');
end
        
gsw_chks = 1;

%% Absolute Salinity (SA) and Preformed Salinity (Sstar)

SA_chck_cast = gsw_SA_from_SP(gsw_cv.SP_chck_cast,gsw_cv.p_chck_cast,gsw_cv.long_chck_cast,gsw_cv.lat_chck_cast);
[ISA_from_SP] = find((gsw_cv.SA_from_SP - SA_chck_cast) >= gsw_cv.SA_from_SP_ca);
if ~isempty(ISA_from_SP)
    try
        cprintf('err','gsw_SA_from_SP:   Failed. Note that this will cause many other programmes in the GSW toolbox to fail.\n');
    catch
        fprintf(2,'gsw_SA_from_SP:   Failed. Note that this will cause many other programmes in the GSW toolbox to fail.\n');
    end
    gsw_chks = 0;
end

Sstar_from_SP = gsw_Sstar_from_SP(gsw_cv.SP_chck_cast,gsw_cv.p_chck_cast,gsw_cv.long_chck_cast,gsw_cv.lat_chck_cast);
[ISstar_from_SP] = find((gsw_cv.Sstar_from_SP - Sstar_from_SP) >= gsw_cv.Sstar_from_SP_ca);
if ~isempty(ISstar_from_SP)
    try
        cprintf('err','gsw_Sstar_from_SP:   Failed\n');
    catch
        fprintf(2,'gsw_Sstar_from_SP:   Failed\n');
    end
    gsw_chks = 0;
end

[SA_SA_Sstar_from_SP, Sstar_SA_Sstar_from_SP] = gsw_SA_Sstar_from_SP(gsw_cv.SP_chck_cast,gsw_cv.p_chck_cast,gsw_cv.long_chck_cast,gsw_cv.lat_chck_cast);
[ISA_Sstar_from_SP] = find((gsw_cv.SA_SA_Sstar_from_SP - SA_SA_Sstar_from_SP) >= gsw_cv.SA_SA_Sstar_from_SP_ca | ...
    (gsw_cv.Sstar_SA_Sstar_from_SP - Sstar_SA_Sstar_from_SP) >= gsw_cv.Sstar_SA_Sstar_from_SP_ca);
if ~isempty(ISA_Sstar_from_SP)
    try
        cprintf('err','gsw_SA_Sstar_from_SP:   Failed\n');
    catch
        fprintf(2,'gsw_SA_Sstar_from_SP:   Failed\n');
    end
    gsw_chks = 0;
end

%% Conservative Temperature (CT)

CT_chck_cast = gsw_CT_from_t(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[ICT_from_t] = find((gsw_cv.CT_from_t - CT_chck_cast) >= gsw_cv.CT_from_t_ca);
if ~isempty(ICT_from_t)
    try
        cprintf('err','gsw_CT_from_t:   Failed. Note that this will cause many other programmes in the GSW toolbox to fail.\n');
    catch
        fprintf(2,'gsw_CT_from_t:   Failed. Note that this will cause many other programmes in the GSW toolbox to fail.\n');
    end
    gsw_chks = 0;
end

%% other conversions between temperatures, salinities, pressure and height

t_from_CT =  gsw_t_from_CT(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[It_from_CT] = find((gsw_cv.t_chck_cast - t_from_CT) >= gsw_cv.t_from_CT_ca);
if ~isempty(It_from_CT)
    try
        cprintf('err','gsw_t_from_CT:   Failed.\n');
    catch
        fprintf(2,'gsw_t_from_CT:   Failed.\n');
    end
    gsw_chks = 0;
end

pt = gsw_pt_from_t(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast,gsw_cv.pr);
[Ipt_from_t] = find((gsw_cv.pt_from_t - pt) >= gsw_cv.pt_from_t_ca);
if ~isempty(Ipt_from_t)
    try
        cprintf('err','gsw_pt_from_t:   Failed\n');
    catch
        fprintf(2,'gsw_pt_from_t:   Failed\n');
    end
    gsw_chks = 0;
end

CT_from_pt = gsw_CT_from_pt(SA_chck_cast,pt);
[ICT_from_pt] = find((gsw_cv.CT_from_pt - CT_from_pt) >= gsw_cv.CT_from_pt_ca);
if ~isempty(ICT_from_pt)
    try
        cprintf('err','gsw_CT_from_pt:   Failed\n');
    catch
        fprintf(2,'gsw_CT_from_pt:   Failed\n');
    end
    gsw_chks = 0;
end

pot_enthalpy = gsw_pot_enthalpy_from_pt(SA_chck_cast,pt);
[Ipot_enthalpy] = find((gsw_cv.pot_enthalpy - pot_enthalpy) >= gsw_cv.pot_enthalpy_ca);
if ~isempty(Ipot_enthalpy)
    try
        cprintf('err','gsw_pot_enthalpy_from_pt:   Failed\n');
    catch
        fprintf(2,'gsw_pot_enthalpy_from_pt:   Failed\n');
    end
    gsw_chks = 0;
end

pt0 = gsw_pt0_from_t(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Ipt0] = find((gsw_cv.pt0 - pt0) >= gsw_cv.pt0_ca);
if ~isempty(Ipt0)
    try
        cprintf('err','gsw_pt0_from_t:   Failed\n');
    catch
        fprintf(2,'gsw_pt0_from_t:   Failed\n');
    end
    gsw_chks = 0;
end

if gsw_chks == 1 ;
    try
        cprintf('text','.');
    catch
        fprintf(1,'.');
    end
end

pt_from_CT = gsw_pt_from_CT(SA_chck_cast,CT_chck_cast);
[Ipt_from_CT] = find((gsw_cv.pt - pt_from_CT) >= gsw_cv.pt_ca);
if ~isempty(Ipt_from_CT)
    try
        cprintf('err','gsw_pt_from_CT:   Failed\n');
    catch
        fprintf(2,'gsw_pt_from_CT:   Failed\n');
    end
    gsw_chks = 0;
end

Sstar_from_SA = gsw_Sstar_from_SA(SA_chck_cast,gsw_cv.p_chck_cast,gsw_cv.long_chck_cast,gsw_cv.lat_chck_cast);
[ISstar_from_SA] = find((gsw_cv.Sstar_from_SA - Sstar_from_SA) >= gsw_cv.Sstar_from_SA_ca);
if ~isempty(ISstar_from_SA)
    try
        cprintf('err','gsw_Sstar_from_SA:   Failed\n');
    catch
        fprintf(2,'gsw_Sstar_from_SA:   Failed\n');
    end
    gsw_chks = 0;
end

SA_from_Sstar = gsw_SA_from_Sstar(SA_chck_cast,gsw_cv.p_chck_cast,gsw_cv.long_chck_cast,gsw_cv.lat_chck_cast);
[ISA_from_Sstar] = find((gsw_cv.SA_from_Sstar - SA_from_Sstar) >= gsw_cv.SA_from_Sstar_ca);
if ~isempty(ISA_from_Sstar)
    try
        cprintf('err','gsw_SA_from_Sstar:   Failed\n');
    catch
        fprintf(2,'gsw_SA_from_Sstar:   Failed\n');
    end
    gsw_chks = 0;
end

SP_from_SA = gsw_SP_from_SA(SA_chck_cast,gsw_cv.p_chck_cast,gsw_cv.long_chck_cast,gsw_cv.lat_chck_cast);
[ISP_from_SA] = find((gsw_cv.SP_chck_cast - SP_from_SA) >= gsw_cv.SP_from_SA_ca);
if ~isempty(ISP_from_SA)
    try
        cprintf('err','gsw_SP_from_SA:   Failed\n');
    catch
        fprintf(2,'gsw_SP_from_SA:   Failed\n');
    end
    gsw_chks = 0;
end

SP_from_Sstar = gsw_SP_from_SA(SA_chck_cast,gsw_cv.p_chck_cast,gsw_cv.long_chck_cast,gsw_cv.lat_chck_cast);
[ISP_from_Sstar] = find((gsw_cv.SP_from_Sstar - SP_from_Sstar) >= gsw_cv.SP_from_Sstar_ca);
if ~isempty(ISP_from_Sstar)
    try
        cprintf('err','gsw_SP_from_Sstar:   Failed\n');
    catch
        fprintf(2,'gsw_SP_from_Sstar:   Failed\n');
    end
    gsw_chks = 0;
end

z_from_p = gsw_z_from_p(gsw_cv.p_chck_cast,gsw_cv.lat_chck_cast);
[Iz_from_p] = find((gsw_cv.z_from_p - z_from_p) >= gsw_cv.z_from_p_ca);
if ~isempty(Iz_from_p)
    try
        cprintf('err','gsw_z_from_p:   Failed\n');
    catch
        fprintf(2,'gsw_z_from_p:   Failed\n');
    end
    gsw_chks = 0;
end

p_from_z = gsw_p_from_z(z_from_p,gsw_cv.lat_chck_cast);
[Ip_from_z] = find((gsw_cv.p_from_z - p_from_z) >= gsw_cv.p_from_z_ca);
if ~isempty(Ip_from_z)
    try
        cprintf('err','gsw_p_from_z:   Failed\n');
    catch
        fprintf(2,'gsw_p_from_z:   Failed\n');
    end
    gsw_chks = 0;
end

t90_from_t68 = gsw_t90_from_t68(gsw_cv.t_chck_cast);
[It90_from_t68] = find((gsw_cv.t90_from_t68 - t90_from_t68) >= gsw_cv.t90_from_t68_ca);
if ~isempty(It90_from_t68)
    try
        cprintf('err','gsw_t90_from_t68:   Failed\n');
    catch
        fprintf(2,'gsw_t90_from_t68:   Failed\n');
    end
    gsw_chks = 0;
end

t90_from_t48 = gsw_t90_from_t48(gsw_cv.t_chck_cast);
[It90_from_t48] = find((gsw_cv.t90_from_t48 - t90_from_t48) >= gsw_cv.t90_from_t48_ca);
if ~isempty(It90_from_t48)
    try
        cprintf('err','gsw_t90_from_t48:   Failed\n');
    catch
        fprintf(2,'gsw_t90_from_t48:   Failed\n');
    end
    gsw_chks = 0;
end

%% density and enthalpy, based on the 25-term expression for density
if gsw_chks == 1 ;
    try
        cprintf('text','.');
    catch
        fprintf(1,'.');
    end
end

rho_CT25 = gsw_rho_CT25(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Irho_CT25] = find((gsw_cv.rho_CT25 - rho_CT25) >= gsw_cv.rho_CT25_ca);
if ~isempty(Irho_CT25)
    try
        cprintf('err','gsw_rho_CT25:   Failed\n');
    catch
        fprintf(2,'gsw_rho_CT25:   Failed\n');
    end
    gsw_chks = 0;
end

[rho_CT25rab, alpha_CT25rab, beta_CT25rab] = gsw_rho_alpha_beta_CT25(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Irho_CT25rab] = find((gsw_cv.rho_CT25rab - rho_CT25rab) >= gsw_cv.rho_CT25rab_ca | ...
    (gsw_cv.alpha_CT25rab - alpha_CT25rab) >= gsw_cv.alpha_CT25rab_ca | ...
    (gsw_cv.beta_CT25rab - beta_CT25rab) >= gsw_cv.beta_CT25rab_ca);
if ~isempty(Irho_CT25rab)
    try
        cprintf('err','gsw_rho_alpha_beta_CT25:   Failed\n');
    catch
        fprintf(2,'gsw_rho_alpha_beta_CT25:   Failed\n');
    end
    gsw_chks = 0;
end

specvol_CT25 = gsw_specvol_CT25(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Ispecvol_CT25] = find((gsw_cv.specvol_CT25 - specvol_CT25) >= gsw_cv.specvol_CT25_ca);
if ~isempty(Ispecvol_CT25)
    try
        cprintf('err','gsw_specvol_CT25:   Failed\n');
    catch
        fprintf(2,'gsw_specvol_CT25:   Failed\n');
    end
    gsw_chks = 0;
end

specvol_anom_CT25 = gsw_specvol_anom_CT25(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Ispecvol_anom_CT25] = find((gsw_cv.specvol_anom_CT25 - specvol_anom_CT25) >= gsw_cv.specvol_anom_CT25_ca);
if ~isempty(Ispecvol_anom_CT25)
    try
        cprintf('err','gsw_specvol_anom_CT25:   Failed\n');
    catch
        fprintf(2,'gsw_specvol_anom_CT25:   Failed\n');
    end
    gsw_chks = 0;
end

enthalpy_CT25 =  gsw_enthalpy_CT25(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Ienthalpy_CT25] = find((gsw_cv.enthalpy_CT25 - enthalpy_CT25) >= gsw_cv.enthalpy_CT25_ca);
if ~isempty(Ienthalpy_CT25)
    try
        cprintf('err','gsw_enthalpy_CT25:   Failed\n');
    catch
        fprintf(2,'gsw_enthalpy_CT25:   Failed\n');
    end
    gsw_chks = 0;
end

enthalpy_diff_CT25 =  gsw_enthalpy_diff_CT25(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast_shallow,gsw_cv.p_chck_cast_deep);
[Ienthalpy_diff_CT25] = find((gsw_cv.enthalpy_diff_CT25 - enthalpy_diff_CT25) >= gsw_cv.enthalpy_diff_CT25_ca);
if ~isempty(Ienthalpy_diff_CT25)
    try
        cprintf('err','gsw_enthalpy_diff_CT25:   Failed\n');
    catch
        fprintf(2,'gsw_enthalpy_diff_CT25:   Failed\n');
    end
    gsw_chks = 0;
end

%% water column properties, based on the 25-term expression for density

[n2, p_mid_n2] = gsw_Nsquared_CT25(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast,gsw_cv.lat_chck_cast);
[INsquared] = find((gsw_cv.n2 - n2) >= gsw_cv.n2_ca | (gsw_cv.p_mid_n2 - p_mid_n2) >= gsw_cv.p_mid_n2_ca);
if ~isempty(INsquared)
    try
        cprintf('err','gsw_Nsquared_CT25:   Failed\n');
    catch
        fprintf(2,'gsw_Nsquared_CT25:   Failed\n');
    end
    gsw_chks = 0;
end

[Tu, Rsubrho, p_mid_TuRsr] = gsw_Turner_Rsubrho_CT25(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[ITurner] = find((gsw_cv.Tu - Tu) >= gsw_cv.Tu_ca | (gsw_cv.Rsubrho - Rsubrho) >= gsw_cv.Rsubrho_ca | ...
    (gsw_cv.p_mid_TuRsr - p_mid_TuRsr) >= gsw_cv.p_mid_TuRsr_ca);
if ~isempty(ITurner)
    try
        cprintf('err','gsw_Turner_Rsubrho_CT25:   Failed\n');
    catch
        fprintf(2,'gsw_Turner_Rsubrho_CT25:   Failed\n');
    end
    gsw_chks = 0;
end

[IPVfN2, p_mid_IPVfN2] = gsw_IPV_vs_fNsquared_ratio_CT25(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast,gsw_cv.pr);
[IIPVfN2] = find((gsw_cv.IPVfN2 - IPVfN2) >= gsw_cv.IPVfN2_ca | ...
    (gsw_cv.p_mid_IPVfN2 - p_mid_IPVfN2) >= gsw_cv.p_mid_IPVfN2_ca);
if ~isempty(IIPVfN2)
    try
        cprintf('err','gsw_IPV_vs_fNsquared_ratio_CT25:   Failed\n');
    catch
        fprintf(2,'gsw_IPV_vs_fNsquared_ratio_CT25:   Failed\n');
    end
    gsw_chks = 0;
end

if gsw_chks == 1 ;
    try
        cprintf('text','.');
    catch
        fprintf(1,'.');
    end
end

%% geostrophic streamfunctions, based on the 25-term expression for density

geo_strf_dyn_height = gsw_geo_strf_dyn_height(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Igeo_strf_dyn_height] = find((gsw_cv.geo_strf_dyn_height - geo_strf_dyn_height) >= gsw_cv.geo_strf_dyn_height_ca);
if ~isempty(Igeo_strf_dyn_height)
    try
        cprintf('err','gsw_geo_strf_dyn_height:   Failed\n');
    catch
        fprintf(2,'gsw_geo_strf_dyn_height:   Failed\n');
    end
    gsw_chks = 0;
end

[geo_strf_dyn_height_pc, dh_pmid] = gsw_geo_strf_dyn_height_pc(SA_chck_cast,CT_chck_cast,gsw_cv.delta_p_chck_cast);
[Igeo_strf_dyn_height_pc] = find((gsw_cv.geo_strf_dyn_height_pc - geo_strf_dyn_height_pc) >= gsw_cv.geo_strf_dyn_height_pc_ca | ...
    (gsw_cv.dh_pmid - dh_pmid) >= gsw_cv.dh_pmid_ca);
if ~isempty(Igeo_strf_dyn_height_pc)
    try
        cprintf('err','gsw_geo_strf_dyn_height_pc:   Failed\n');
    catch
        fprintf(2,'gsw_geo_strf_dyn_height_pc:   Failed\n');
    end
    gsw_chks = 0;
end

geo_strf_McD_Klocker = gsw_geo_strf_McD_Klocker(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast,gsw_cv.Neutral_Density,gsw_cv.p_Neutral_Density);
[Igeo_strf_McD_Klocker] = find((gsw_cv.geo_strf_McD_Klocker - geo_strf_McD_Klocker) >= gsw_cv.geo_strf_McD_Klocker_ca);
if ~isempty(Igeo_strf_McD_Klocker)
    try
        cprintf('err','gsw_geo_strf_McD_Klocker:   Failed\n');
    catch
        fprintf(2,'gsw_geo_strf_McD_Klocker:   Failed\n');
    end
    gsw_chks = 0;
end

[geo_strf_McD_Klocker_pc, mk_p_mid] = gsw_geo_strf_McD_Klocker_pc(SA_chck_cast,CT_chck_cast,gsw_cv.delta_p_chck_cast,gsw_cv.Neutral_Density(1),3);
[Igeo_strf_McD_Klocker_pc] = find((gsw_cv.geo_strf_McD_Klocker_pc - geo_strf_McD_Klocker_pc) >= gsw_cv.geo_strf_McD_Klocker_pc_ca |...
    (gsw_cv.mk_p_mid - mk_p_mid) >= gsw_cv.mk_p_mid_ca);
if ~isempty(Igeo_strf_McD_Klocker_pc)
    try
        cprintf('err','gsw_geo_strf_McD_Klocker_pc:   Failed\n');
    catch
        fprintf(2,'gsw_geo_strf_McD_Klocker_pc:   Failed\n');
    end
    gsw_chks = 0;
end

geo_strf_Montgomery = gsw_geo_strf_Montgomery(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Igeo_strf_Montgomery] = find((gsw_cv.geo_strf_Montgomery - geo_strf_Montgomery) >= gsw_cv.geo_strf_Montgomery_ca);
if ~isempty(Igeo_strf_Montgomery)
    try
        cprintf('err','gsw_geo_strf_Montgomery:   Failed\n');
    catch
        fprintf(2,'gsw_geo_strf_Montgomery:   Failed\n');
    end
    gsw_chks = 0;
end

geo_strf_Cunningham = gsw_geo_strf_Cunningham(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Igeo_strf_Cunningham] = find((gsw_cv.geo_strf_Cunningham - geo_strf_Cunningham) >= gsw_cv.geo_strf_Cunningham_ca);
if ~isempty(Igeo_strf_Cunningham)
    try
        cprintf('err','gsw_geo_strf_Cunningham:   Failed\n');
    catch
        fprintf(2,'gsw_geo_strf_Cunningham:   Failed\n');
    end
    gsw_chks = 0;
end

[geo_str_velocity, gv_mid_lat, gv_mid_long] = gsw_geostrophic_velocity(geo_strf_dyn_height,gsw_cv.long_chck_cast,gsw_cv.lat_chck_cast,gsw_cv.p_chck_cast);
[Igeostrophic_velo] = find((gsw_cv.geo_str_velocity - geo_str_velocity) >= gsw_cv.geo_str_velocity_ca | ...
    (gsw_cv.gv_mid_lat - gv_mid_lat) >= gsw_cv.gv_mid_lat_ca  | ...
    (gsw_cv.gv_mid_long - gv_mid_long) >= gsw_cv.gv_mid_long_ca);
if ~isempty(Igeostrophic_velo)
    try
        cprintf('err','gsw_geostrophic_velocity:   Failed\n');
    catch
        fprintf(2,'gsw_geostrophic_velocity:   Failed\n');
    end
    gsw_chks = 0;
end

%% neutral and non-linear properties, based on the 25-term expression for density

cabbeling_CT25 = gsw_cabbeling_CT25(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Icabbeling_CT25] = find((gsw_cv.cabbeling_CT25 - cabbeling_CT25) >= gsw_cv.cabbeling_CT25_ca);
if ~isempty(Icabbeling_CT25)
    try
        cprintf('err','gsw_cabbeling_CT25:   Failed\n');
    catch
        fprintf(2,'gsw_cabbeling_CT25:   Failed\n');
    end
    gsw_chks = 0;
end

thermobaric_CT25 = gsw_thermobaric_CT25(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Ithermobaric_CT25] = find((gsw_cv.thermobaric_CT25 - thermobaric_CT25) >= gsw_cv.thermobaric_CT25_ca);
if ~isempty(Ithermobaric_CT25)
    try
        cprintf('err','gsw_thermobaric_CT25:   Failed\n');
    catch
        fprintf(2,'gsw_thermobaric_CT25:   Failed\n');
    end
    gsw_chks = 0;
end

if gsw_chks == 1 ;
    try
        cprintf('text','.');
    catch
        fprintf(1,'.');
    end
end

isopycnal_slope_ratio_CT25 = gsw_isopycnal_slope_ratio_CT25(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast,gsw_cv.pr);
[Iisopycnal_slope_ratio_CT25] = find((gsw_cv.isopycnal_slope_ratio_CT25 - isopycnal_slope_ratio_CT25) >= gsw_cv.isopycnal_slope_ratio_CT25_ca);
if ~isempty(Iisopycnal_slope_ratio_CT25)
    try
        cprintf('err','gsw_isopycnal_slope_ratio_CT25:   Failed\n');
    catch
        fprintf(2,'gsw_isopycnal_slope_ratio_CT25:   Failed\n');
    end
    gsw_chks = 0;
end

[G_CT_CT25, p_mid_G_CT_CT25] = gsw_isopycnal_vs_ntp_CT_ratio_CT25(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast,gsw_cv.pr);
[IG_CT] = find((gsw_cv.G_CT_CT25 - G_CT_CT25) >= gsw_cv.G_CT_CT25_ca | ...
    (gsw_cv.p_mid_G_CT_CT25 - p_mid_G_CT_CT25) >= gsw_cv.p_mid_G_CT_CT25_ca);
if ~isempty(IG_CT)
    try
        cprintf('err','gsw_isopycnal_vs_ntp_CT_ratio_CT25:   Failed\n');
    catch
        fprintf(2,'gsw_isopycnal_vs_ntp_CT_ratio_CT25:   Failed\n');
    end
    gsw_chks = 0;
end

ntpptCT_CT25 = gsw_ntp_pt_vs_CT_ratio_CT25(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[IntpptCT_CT25] = find((gsw_cv.ntpptCT_CT25 - ntpptCT_CT25) >= gsw_cv.ntpptCT_CT25_ca);
if ~isempty(IntpptCT_CT25)
    try
        cprintf('err','gsw_ntp_pt_vs_CT_ratio_CT25:   Failed\n');
    catch
        fprintf(2,'gsw_ntp_pt_vs_CT_ratio_CT25:   Failed\n');
    end
    gsw_chks = 0;
end

%% basic thermodynamic properties

rho = gsw_rho(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Irho] = find((gsw_cv.rho - rho) >= gsw_cv.rho_ca);
if ~isempty(Irho)
    try
        cprintf('err','gsw_rho:   Failed\n');
    catch
        fprintf(2,'gsw_rho:   Failed\n');
    end
    gsw_chks = 0;
end

pot_rho = gsw_pot_rho(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast,gsw_cv.pr);
[Ipot_rho] = find((gsw_cv.pot_rho - pot_rho) >= gsw_cv.pot_rho_ca);
if ~isempty(Ipot_rho)
    try
        cprintf('err','gsw_pot_rho:   Failed\n');
    catch
        fprintf(2,'gsw_pot_rho:   Failed\n');
    end
    gsw_chks = 0;
end

specvol = gsw_specvol(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Ispecvol] = find((gsw_cv.specvol - specvol) >= gsw_cv.specvol_ca);
if ~isempty(Ispecvol)
    try
        cprintf('err','gsw_specvol:   Failed\n');
    catch
        fprintf(2,'gsw_specvol:   Failed\n');
    end
    gsw_chks = 0;
end

specvol_anom = gsw_specvol_anom(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Ispecvol_anom] = find((gsw_cv.specvol_anom - specvol_anom) >= gsw_cv.specvol_anom_ca);
if ~isempty(Ispecvol_anom)
    try
        cprintf('err','gsw_specvol_anom:   Failed\n');
    catch
        fprintf(2,'gsw_specvol_anom:   Failed\n');
    end
    gsw_chks = 0;
end

alpha_wrt_CT = gsw_alpha_wrt_CT(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Ialpha_wrt_CT] = find((gsw_cv.alpha_wrt_CT - alpha_wrt_CT) >= gsw_cv.alpha_wrt_CT_ca);
if ~isempty(Ialpha_wrt_CT)
    try
        cprintf('err','gsw_alpha_wrt_CT:   Failed\n');
    catch
        fprintf(2,'gsw_alpha_wrt_CT:   Failed\n');
    end
    gsw_chks = 0;
end

alpha_wrt_pt = gsw_alpha_wrt_pt(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Ialpha_wrt_pt] = find((gsw_cv.alpha_wrt_pt - alpha_wrt_pt) >= gsw_cv.alpha_wrt_pt_ca);
if ~isempty(Ialpha_wrt_pt)
    try
        cprintf('err','gsw_alpha_wrt_pt:   Failed\n');
    catch
        fprintf(2,'gsw_alpha_wrt_pt:   Failed\n');
    end
    gsw_chks = 0;
end

alpha_wrt_t = gsw_alpha_wrt_t(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Ialpha_wrt_t] = find((gsw_cv.alpha_wrt_t - alpha_wrt_t) >= gsw_cv.alpha_wrt_t_ca);
if ~isempty(Ialpha_wrt_t)
    try
        cprintf('err','gsw_alpha_wrt_t:   Failed\n');
    catch
        fprintf(2,'gsw_alpha_wrt_t:   Failed\n');
    end
    gsw_chks = 0;
end

if gsw_chks == 1 ;
    try
        cprintf('text','.');
    catch
        fprintf(1,'.');
    end
end

beta_const_CT = gsw_beta_const_CT(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Ibeta_const_CT] = find((gsw_cv.beta_const_CT - beta_const_CT) >= gsw_cv.beta_const_CT_ca);
if ~isempty(Ibeta_const_CT)
    try
        cprintf('err','gsw_beta_const_CT:   Failed\n');
    catch
        fprintf(2,'gsw_beta_const_CT:   Failed\n');
    end
    gsw_chks = 0;
end

beta_const_pt = gsw_beta_const_pt(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Ibeta_const_pt] = find((gsw_cv.beta_const_pt - beta_const_pt) >= gsw_cv.beta_const_pt_ca);
if ~isempty(Ibeta_const_pt)
    try
        cprintf('err','gsw_beta_const_pt:   Failed\n');
    catch
        fprintf(2,'gsw_beta_const_pt:   Failed\n');
    end
    gsw_chks = 0;
end

beta_const_t = gsw_beta_const_t(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Ibeta_const_t] = find((gsw_cv.beta_const_t - beta_const_t) >= gsw_cv.beta_const_t_ca);
if ~isempty(Ibeta_const_t)
    try
        cprintf('err','gsw_beta_const_t:   Failed\n');
    catch
        fprintf(2,'gsw_beta_const_t:   Failed\n');
    end
    gsw_chks = 0;
end

entropy = gsw_entropy(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Ientropy] = find((gsw_cv.entropy - entropy) >= gsw_cv.entropy_ca);
if ~isempty(Ientropy)
    try
        cprintf('err','gsw_entropy:   Failed\n');
    catch
        fprintf(2,'gsw_entropy:   Failed\n');
    end
    gsw_chks = 0;
end

internal_energy = gsw_internal_energy(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Iinternal_energy] = find((gsw_cv.internal_energy - internal_energy) >= gsw_cv.internal_energy_ca);
if ~isempty(Iinternal_energy)
    try
        cprintf('err','gsw_internal_energy:   Failed\n');
    catch
        fprintf(2,'gsw_internal_energy:   Failed\n');
    end
    gsw_chks = 0;
end

enthalpy = gsw_enthalpy(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Ienthalpy] = find((gsw_cv.enthalpy - enthalpy) >= gsw_cv.enthalpy_ca);
if ~isempty(Ienthalpy)
    try
        cprintf('err','gsw_enthalpy:   Failed\n');
    catch
        fprintf(2,'gsw_enthalpy:   Failed\n');
    end
    gsw_chks = 0;
end

cp = gsw_cp(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Icp] = find((gsw_cv.cp - cp) >= gsw_cv.cp_ca);
if ~isempty(Icp)
    try
        cprintf('err','gsw_cp:   Failed\n');
    catch
        fprintf(2,'gsw_cp:   Failed\n');
    end
    gsw_chks = 0;
end

isochoric_heat_cap = gsw_isochoric_heat_cap(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Iisochoric_heat_cap] = find((gsw_cv.isochoric_heat_cap - isochoric_heat_cap) >= gsw_cv.isochoric_heat_cap_ca);
if ~isempty(Iisochoric_heat_cap)
    try
        cprintf('err','gsw_isochoric_heat_cap:   Failed\n');
    catch
        fprintf(2,'gsw_isochoric_heat_cap:   Failed\n');
    end
    gsw_chks = 0;
end

chem_potential =  gsw_chem_potential_relative(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Ichem_potential] = find((gsw_cv.chem_potential - chem_potential) >= gsw_cv.chem_potential_ca);
if ~isempty(Ichem_potential)
    try
        cprintf('err','gsw_chem_potential_relative:   Failed\n');
    catch
        fprintf(2,'gsw_chem_potential_relative:   Failed\n');
    end
    gsw_chks = 0;
end

chem_potential_water =  gsw_chem_potential_water(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Ichem_potential_water] = find((gsw_cv.chem_potential_water - chem_potential_water) >= gsw_cv.chem_potential_water_ca);
if ~isempty(Ichem_potential_water)
    try
        cprintf('err','gsw_chem_potential_water:   Failed\n');
    catch
        fprintf(2,'gsw_chem_potential_water:   Failed\n');
    end
    gsw_chks = 0;
end

if gsw_chks == 1 ;
    try
        cprintf('text','.');
    catch
        fprintf(1,'.');
    end
end

chem_potential_salt =  gsw_chem_potential_salt(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Ichem_potential_salt] = find((gsw_cv.chem_potential_salt - chem_potential_salt) >= gsw_cv.chem_potential_salt_ca);
if ~isempty(Ichem_potential_salt)
    try
        cprintf('err','gsw_chem_potential_salt:   Failed\n');
    catch
        fprintf(2,'gsw_chem_potential_salt:   Failed\n');
    end
    gsw_chks = 0;
end

Helmholtz_energy = gsw_Helmholtz_energy(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[IHelmholtz_energy] = find((gsw_cv.Helmholtz_energy - Helmholtz_energy) >= gsw_cv.Helmholtz_energy_ca);
if ~isempty(IHelmholtz_energy)
    try
        cprintf('err','gsw_Helmholtz_energy:   Failed\n');
    catch
        fprintf(2,'gsw_Helmholtz_energy:   Failed\n');
    end
    gsw_chks = 0;
end

sound_speed = gsw_sound_speed(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Isound_speed] = find((gsw_cv.sound_speed - sound_speed) >= gsw_cv.sound_speed_ca);
if ~isempty(Isound_speed)
    try
        cprintf('err','gsw_sound_speed:   Failed\n');
    catch
        fprintf(2,'gsw_sound_speed:   Failed\n');
    end
    gsw_chks = 0;
end

kappa = gsw_kappa(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Ikappa] = find((gsw_cv.kappa - kappa) >= gsw_cv.kappa_ca);
if ~isempty(Ikappa)
    try
        cprintf('err','gsw_kappa:   Failed\n');
    catch
        fprintf(2,'gsw_kappa:   Failed\n');
    end
    gsw_chks = 0;
end

kappa_const_t = gsw_kappa_const_t(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Ikappa_const_t] = find((gsw_cv.kappa_const_t - kappa_const_t) >= gsw_cv.kappa_const_t_ca);
if ~isempty(Ikappa_const_t)
    try
        cprintf('err','gsw_kappa_const_t:   Failed\n');
    catch
        fprintf(2,'gsw_kappa_const_t:   Failed\n');
    end
    gsw_chks = 0;
end

adiabatic_lapse_rate = gsw_adiabatic_lapse_rate(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Iadiabatic_lapse_rate] = find((gsw_cv.adiabatic_lapse_rate - adiabatic_lapse_rate) >= gsw_cv.adiabatic_lapse_rate_ca);
if ~isempty(Iadiabatic_lapse_rate)
    try
        cprintf('err','gsw_adiabatic_lapse_rate:   Failed\n');
    catch
        fprintf(2,'gsw_adiabatic_lapse_rate:   Failed\n');
    end
    gsw_chks = 0;
end

molality = gsw_molality(SA_chck_cast);
[Imolality] = find((gsw_cv.molality - molality) >= gsw_cv.molality_ca);
if ~isempty(Imolality)
    try
        cprintf('err','gsw_molality:   Failed\n');
    catch
        fprintf(2,'gsw_molality:   Failed\n');
    end
    gsw_chks = 0;
end

ionic_strength = gsw_ionic_strength(SA_chck_cast);
[Iionic_strength] = find((gsw_cv.ionic_strength - ionic_strength) >= gsw_cv.ionic_strength_ca);
if ~isempty(Iionic_strength)
    try
        cprintf('err','gsw_ionic_strength:   Failed\n');
    catch
        fprintf(2,'gsw_ionic_strength:   Failed\n');
    end
    gsw_chks = 0;
end

osmotic_coefficient = gsw_osmotic_coefficient(SA_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Iosmotic_coefficient] = find((gsw_cv.osmotic_coefficient - osmotic_coefficient) >= gsw_cv.osmotic_coefficient_ca);
if ~isempty(Iosmotic_coefficient)
    try
        cprintf('err','gsw_osmotic_coefficient:   Failed\n');
    catch
        fprintf(2,'gsw_osmotic_coefficient:   Failed\n');
    end
    gsw_chks = 0;
end

[t_maxden, pt_maxden, CT_maxden] = gsw_temps_maxdensity(SA_chck_cast,gsw_cv.p_chck_cast);
[Itemps_maxd] = find((gsw_cv.t_maxden - t_maxden) >= gsw_cv.t_maxden_ca |  ...
    (gsw_cv.pt_maxden - pt_maxden) >= gsw_cv.pt_maxden_ca | ...
    (gsw_cv.CT_maxden - CT_maxden) >= gsw_cv.CT_maxden_ca);
if ~isempty(Itemps_maxd)
    try
        cprintf('err','gsw_temps_maxdensity:   Failed\n');
    catch
        fprintf(2,'gsw_temps_maxdensity:   Failed\n');
    end
    gsw_chks = 0;
end

if gsw_chks == 1 ;
    try
        cprintf('text','.');
    catch
        fprintf(1,'.');
    end
end

%% basic thermodynamic properties in terms of CT and pt

rho_CT = gsw_rho_CT(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Irho_CT] = find((gsw_cv.rho_CT - rho_CT) >= gsw_cv.rho_CT_ca);
if ~isempty(Irho_CT)
    try
        cprintf('err','gsw_rho_CT:   Failed\n');
    catch
        fprintf(2,'gsw_rho_CT:   Failed\n');
    end
    gsw_chks = 0;
end

[rho_CTrab, alpha_CTrab, beta_CTrab] = gsw_rho_alpha_beta_CT(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Irho_alpha_beta] = find((gsw_cv.rho_CTrab - rho_CTrab) >= gsw_cv.rho_CTrab_ca |...
    (gsw_cv.alpha_CTrab - alpha_CTrab) >= gsw_cv.alpha_CTrab_ca |...
    (gsw_cv.beta_CTrab - beta_CTrab) >= gsw_cv.beta_CTrab_ca);
if ~isempty(Irho_alpha_beta)
    try
        cprintf('err','gsw_rho_alpha_beta_CT:   Failed\n');
    catch
        fprintf(2,'gsw_rho_alpha_beta_CT:   Failed\n');
    end
    gsw_chks = 0;
end

specvol_CT = gsw_specvol_CT25(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Ispecvol_CT] = find((gsw_cv.specvol_CT - specvol_CT) >= gsw_cv.specvol_CT_ca);
if ~isempty(Ispecvol_CT)
    try
        cprintf('err','gsw_specvol_CT25:   Failed\n');
    catch
        fprintf(2,'gsw_specvol_CT25:   Failed\n');
    end
    gsw_chks = 0;
end

specvol_anom_CT = gsw_specvol_anom_CT(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Ispecvol_anom_CT] = find((gsw_cv.specvol_anom_CT - specvol_anom_CT) >= gsw_cv.specvol_anom_CT_ca);
if ~isempty(Ispecvol_anom_CT)
    try
        cprintf('err','gsw_specvol_anom_CT:   Failed\n');
    catch
        fprintf(2,'gsw_specvol_anom_CT:   Failed\n');
    end
    gsw_chks = 0;
end

sigma0_CT = gsw_sigma0_CT(SA_chck_cast,CT_chck_cast);
[Isigma0_CT] = find((gsw_cv.sigma0_CT - sigma0_CT) >= gsw_cv.sigma0_CT_ca);
if ~isempty(Isigma0_CT)
    try
        cprintf('err','gsw_sigma0_CT:   Failed\n');
    catch
        fprintf(2,'gsw_sigma0_CT:   Failed\n');
    end
    gsw_chks = 0;
end

sigma1_CT = gsw_sigma1_CT(SA_chck_cast,CT_chck_cast);
[Isigma1_CT] = find((gsw_cv.sigma1_CT - sigma1_CT) >= gsw_cv.sigma1_CT_ca);
if ~isempty(Isigma1_CT)
    try
        cprintf('err','gsw_sigma1_CT:   Failed\n');
    catch
        fprintf(2,'gsw_sigma1_CT:   Failed\n');
    end
    gsw_chks = 0;
end

sigma2_CT = gsw_sigma2_CT(SA_chck_cast,CT_chck_cast);
[Isigma2_CT] = find((gsw_cv.sigma2_CT - sigma2_CT) >= gsw_cv.sigma2_CT_ca);
if ~isempty(Isigma2_CT)
    try
        cprintf('err','gsw_sigma2_CT:   Failed\n');
    catch
        fprintf(2,'gsw_sigma2_CT:   Failed\n');
    end
    gsw_chks = 0;
end

sigma3_CT = gsw_sigma3_CT(SA_chck_cast,CT_chck_cast);
[Isigma3_CT] = find((gsw_cv.sigma3_CT - sigma3_CT) >= gsw_cv.sigma3_CT_ca);
if ~isempty(Isigma3_CT)
    try
        cprintf('err','gsw_sigma3_CT:   Failed\n');
    catch
        fprintf(2,'gsw_sigma3_CT:   Failed\n');
    end
    gsw_chks = 0;
end

sigma4_CT = gsw_sigma4_CT(SA_chck_cast,CT_chck_cast);
[Isigma4_CT] = find((gsw_cv.sigma4_CT - sigma4_CT) >= gsw_cv.sigma4_CT_ca);
if ~isempty(Isigma4_CT)
    try
        cprintf('err','gsw_sigma4_CT:   Failed\n');
    catch
        fprintf(2,'gsw_sigma4_CT:   Failed\n');
    end
    gsw_chks = 0;
end

enthalpy_CT =  gsw_enthalpy_CT(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Ienthalpy_CT] = find((gsw_cv.enthalpy_CT - enthalpy_CT) >= gsw_cv.enthalpy_CT_ca);
if ~isempty(Ienthalpy_CT)
    try
        cprintf('err','gsw_enthalpy_CT:   Failed\n');
    catch
        fprintf(2,'gsw_enthalpy_CT:   Failed\n');
    end
    gsw_chks = 0;
end

if gsw_chks == 1 ;
    try
        cprintf('text','.');
    catch
        fprintf(1,'.');
    end
end

enthalpy_diff_CT =  gsw_enthalpy_diff_CT(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast_shallow,gsw_cv.p_chck_cast_deep);
[Ienthalpy_diff_CT] = find((gsw_cv.enthalpy_diff_CT - enthalpy_diff_CT) >= gsw_cv.enthalpy_diff_CT_ca);
if ~isempty(Ienthalpy_diff_CT)
    try
        cprintf('err','gsw_enthalpy_diff_CT:   Failed\n');
    catch
        fprintf(2,'gsw_enthalpy_diff_CT:   Failed\n');
    end
    gsw_chks = 0;
end
 
entropy_from_pt =  gsw_entropy_from_pt(SA_chck_cast,pt);
[Ientropy_from_pt] = find((gsw_cv.entropy_from_pt - entropy_from_pt) >= gsw_cv.entropy_from_pt_ca);
if ~isempty(Ientropy_from_pt)
    try
        cprintf('err','gsw_entropy_from_pt:   Failed\n');
    catch
        fprintf(2,'gsw_entropy_from_pt:   Failed\n');
    end
    gsw_chks = 0;
end

entropy_from_CT =  gsw_entropy_from_CT(SA_chck_cast,CT_chck_cast);
[Ientropy_from_CT] = find((gsw_cv.entropy_from_CT - entropy_from_CT) >= gsw_cv.entropy_from_CT_ca);
if ~isempty(Ientropy_from_CT)
    try
        cprintf('err','gsw_entropy_from_CT:   Failed\n');
    catch
        fprintf(2,'gsw_entropy_from_CT:   Failed\n');
    end
    gsw_chks = 0;
end

CT_from_entropy =  gsw_CT_from_entropy(SA_chck_cast,entropy);
[ICT_from_entropy] = find((gsw_cv.CT_from_entropy - CT_from_entropy) >= gsw_cv.CT_from_entropy_ca);
if ~isempty(ICT_from_entropy)
    try
        cprintf('err','gsw_CT_from_entropy:   Failed\n');
    catch
        fprintf(2,'gsw_CT_from_entropy:   Failed\n');
    end
    gsw_chks = 0;
end

pt_from_entropy =  gsw_pt_from_entropy(SA_chck_cast,entropy);
[Ipt_from_entropy] = find((gsw_cv.pt_from_entropy - pt_from_entropy) >= gsw_cv.pt_from_entropy_ca);
if ~isempty(Ipt_from_entropy)
    try
        cprintf('err','gsw_pt_from_entropy:   Failed\n');
    catch
        fprintf(2,'gsw_pt_from_entropy:   Failed\n');
    end
    gsw_chks = 0;
end

%% derivatives of enthalpy, entropy, CT and pt

[CT_SA, CT_pt] = gsw_CT_first_derivatives(SA_chck_cast,pt);
[ICT_first_deriv] = find((gsw_cv.CT_SA - CT_SA) >= gsw_cv.CT_SA_ca | ...
    (gsw_cv.CT_pt - CT_pt) >= gsw_cv.CT_pt_ca);
if ~isempty(ICT_first_deriv)
    try
        cprintf('err','gsw_CT_first_derivatives:   Failed\n');
    catch
        fprintf(2,'gsw_CT_first_derivatives:   Failed\n');
    end
    gsw_chks = 0;
end

[CT_SA_SA, CT_SA_pt, CT_pt_pt] = gsw_CT_second_derivatives(SA_chck_cast,pt);
[ICT_second_deriv] = find((gsw_cv.CT_SA_SA - CT_SA_SA) >= gsw_cv.CT_SA_SA_ca | ...
    (gsw_cv.CT_SA_pt - CT_SA_pt) >= gsw_cv.CT_SA_pt_ca | ...
    (gsw_cv.CT_pt_pt - CT_pt_pt) >= gsw_cv.CT_pt_pt_ca);
if ~isempty(ICT_second_deriv)
    try
        cprintf('err','gsw_CT_second_derivatives:   Failed\n');
    catch
        fprintf(2,'gsw_CT_second_derivatives:   Failed\n');
    end
    gsw_chks = 0;
end

[h_SA, h_CT, h_P] = gsw_enthalpy_first_derivatives(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Ienthalpy_first_deriv] = find((gsw_cv.h_SA - h_SA) >= gsw_cv.h_SA_ca | ...
    (gsw_cv.h_CT - h_CT) >= gsw_cv.h_CT_ca | ...
    (gsw_cv.h_P - h_P) >= gsw_cv.h_P_ca);
if ~isempty(Ienthalpy_first_deriv)
    try
        cprintf('err','gsw_enthalpy_first_derivatives:   Failed\n');
    catch
        fprintf(2,'gsw_enthalpy_first_derivatives:   Failed\n');
    end
    gsw_chks = 0;
end

[h_SA_SA, h_SA_CT, h_CT_CT] = gsw_enthalpy_second_derivatives(SA_chck_cast,CT_chck_cast,gsw_cv.p_chck_cast);
[Ienthalpy_second_deriv] = find((gsw_cv.h_SA_SA - h_SA_SA) >= gsw_cv.h_SA_SA_ca  | ...
    (gsw_cv.h_SA_CT - h_SA_CT) >= gsw_cv.h_SA_CT_ca | ...
    (gsw_cv.h_CT_CT - h_CT_CT) >= gsw_cv.h_CT_CT_ca);
if ~isempty(Ienthalpy_second_deriv)
    try
        cprintf('err','gsw_enthalpy_second_derivatives:   Failed\n');
    catch
        fprintf(2,'gsw_enthalpy_second_derivatives:   Failed\n');
    end
    gsw_chks = 0;
end

[eta_SA, eta_CT] = gsw_entropy_first_derivatives(SA_chck_cast,CT_chck_cast);
[Ientropy_first_deriv] = find((gsw_cv.eta_SA - eta_SA) >= gsw_cv.eta_SA_ca | ...
    (gsw_cv.eta_CT - eta_CT) >= gsw_cv.eta_CT_ca);
if ~isempty(Ientropy_first_deriv)
    try
        cprintf('err','Ienthalpy_second_deriv:   Failed\n');
    catch
        fprintf(2,'Ienthalpy_second_deriv:   Failed\n');
    end
    gsw_chks = 0;
end

if gsw_chks == 1 ;
    try
        cprintf('text','.');
    catch
        fprintf(1,'.');
    end
end

[eta_SA_SA, eta_SA_CT, eta_CT_CT] = gsw_entropy_second_derivatives(SA_chck_cast,CT_chck_cast);
[Ientropy_second_deriv] = find(((gsw_cv.eta_SA_SA - eta_SA_SA)) >= gsw_cv.eta_SA_SA_ca |...
    (gsw_cv.eta_SA_CT - eta_SA_CT) >= gsw_cv.eta_SA_CT_ca |...
    (gsw_cv.eta_CT_CT - eta_CT_CT) >= gsw_cv.eta_CT_CT_ca);
if ~isempty(Ientropy_second_deriv)
    try
        cprintf('err','gsw_entropy_second_derivatives:   Failed\n');
    catch
        fprintf(2,'gsw_entropy_second_derivatives:   Failed\n');
    end
    gsw_chks = 0;
end

[pt_SA, pt_CT] = gsw_pt_first_derivatives(SA_chck_cast,CT_chck_cast);
[Ipt_first_deriv] = find((gsw_cv.pt_SA - pt_SA) >= gsw_cv.pt_SA_ca |...
    (gsw_cv.pt_CT - pt_CT) >= gsw_cv.pt_CT_ca);
if ~isempty(Ipt_first_deriv)
    try
        cprintf('err','gsw_pt_first_derivatives:   Failed\n');
    catch
        fprintf(2,'gsw_pt_first_derivatives:   Failed\n');
    end
    gsw_chks = 0;
end

[pt_SA_SA, pt_SA_CT, pt_CT_CT] = gsw_pt_second_derivatives(SA_chck_cast,CT_chck_cast);
[Ipt_second_deriv] = find((gsw_cv.pt_SA_SA - pt_SA_SA) >= gsw_cv.pt_SA_SA_ca  | ...
    (gsw_cv.pt_SA_CT - pt_SA_CT) >= gsw_cv.pt_SA_CT_ca | ...
    (gsw_cv.pt_CT_CT - pt_CT_CT) >= gsw_cv.pt_CT_CT_ca);
if ~isempty(Ipt_second_deriv)
    try
        cprintf('err','gsw_pt_second_derivatives:   Failed\n');
    catch
        fprintf(2,'gsw_pt_second_derivatives:   Failed\n');
    end
    gsw_chks = 0;
end

%% planet earth properties

f = gsw_f(gsw_cv.lat_chck_cast);
[If] = find((gsw_cv.f - f) >= gsw_cv.f_ca);
if ~isempty(If)
    try
        cprintf('err','gsw_f:   Failed\n');
    catch
        fprintf(2,'gsw_f:   Failed\n');
    end
    gsw_chks = 0;
end

grav = gsw_grav(gsw_cv.lat_chck_cast,gsw_cv.p_chck_cast);
[Igrav] = find((gsw_cv.grav - grav) >= gsw_cv.grav_ca);
if ~isempty(Igrav)
    try
        cprintf('err','gsw_grav:   Failed\n');
    catch
        fprintf(2,'gsw_grav:   Failed\n');
    end
    gsw_chks = 0;
end

distance = gsw_distance(gsw_cv.long_chck_cast,gsw_cv.lat_chck_cast,gsw_cv.p_chck_cast);
[Idistance] = find((gsw_cv.distance - distance) >= gsw_cv.distance_ca);
if ~isempty(Idistance)
    try
        cprintf('err','gsw_distance:   Failed\n');
    catch
        fprintf(2,'gsw_distance:   Failed\n');
    end
    gsw_chks = 0;
end

%% Absolute Salinity from direct density measurements:- a laboratory function

SA_from_rho = gsw_SA_from_rho(rho,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[ISA_from_rho] = find((gsw_cv.SA_from_rho - SA_from_rho) >= gsw_cv.SA_from_rho_ca);
if ~isempty(ISA_from_rho)
    try
        cprintf('err','gsw_SA_from_rho:   Failed\n');
    catch
        fprintf(2,'gsw_SA_from_rho:   Failed\n');
    end
    gsw_chks = 0;
end

sigma0_pt = gsw_sigma0_pt(SA_chck_cast,pt0);
[Isigma0_pt] = find((gsw_cv.sigma0_pt - sigma0_pt) >= gsw_cv.sigma0_pt_ca);
if ~isempty(Isigma0_pt)
    try
        cprintf('err','gsw_sigma0_pt:   Failed\n');
    catch
        fprintf(2,'gsw_sigma0_pt:   Failed\n');
    end
    gsw_chks = 0;
end

%% Practical Salinity (SP):- PSS-78

cndr = gsw_cndr_from_SP(gsw_cv.SP_chck_cast,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[Icndr] = find((gsw_cv.cndr - cndr) >= gsw_cv.cndr_ca);
if ~isempty(Icndr)
    try
        cprintf('err','gsw_cndr_from_SP:   Failed\n');
    catch
        fprintf(2,'gsw_cndr_from_SP:   Failed\n');
    end
    gsw_chks = 0;
end

SP_from_cndr = gsw_SP_from_cndr(cndr,gsw_cv.t_chck_cast,gsw_cv.p_chck_cast);
[ISP_from_cndr] = find((gsw_cv.SP_from_cndr - SP_from_cndr) >= gsw_cv.SP_from_cndr_ca);
if ~isempty(ISP_from_cndr)
    try
        cprintf('err','gsw_SP_from_cndr:   Failed\n');
    catch
        fprintf(2,'gsw_SP_from_cndr:   Failed\n');
    end
    gsw_chks = 0;
end

if gsw_chks == 1 ;
    try
        cprintf('text',' Finished.\n');
        cprintf('text',' \n');
    catch
        fprintf(1,' Finished.\n');
        fprintf(1,'\n');
    end
end

if gsw_chks == 0
    try
        cprintf('err','Your installation of the Gibbs SeaWater (GSW) Oceanographic Toolbox has errors !\n');
    catch
        fprintf(2,'Your installation of the Gibbs SeaWater (GSW) Oceanographic Toolbox has errors !\n');
    end
else
    try
        cprintf('comment','Well done! The gsw_check_values function confirms that the \n')
        cprintf('comment','Gibbs SeaWater (GSW) Oceanographic Toolbox is installed correctly.\n');
        cprintf('text',' \n')
        cprintf('text','Press enter to continue. \n')
        pause
        cprintf('text',' \n')
        cprintf('keywords','A short demonstration of the GSW Oceanographic Toolbox now follows.\n');
        pause(3)
        cprintf('comment','The following vertical profile, from the North Pacific, is of\n');
        cprintf('comment','Practical Salinity, SP, and in situ temperature, t, as a function\n');
        cprintf('comment','of pressure, p,\n');
        pause(6)
        cprintf('text','%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %s \n' ,'SP = [',gsw_cv.SP_chck_cast([1,22,29:4:45],1)',']');
        cprintf('text','%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %s \n' ,'t  = [',gsw_cv.t_chck_cast([1,22,29:4:45],1)',']');
        cprintf('text','%s %7.0f  %7.0f  %7.0f  %7.0f  %7.0f  %7.0f  %7.0f %s \n' ,'p  = [',gsw_cv.p_chck_cast([1,22,29:4:45],1)',']');
        cprintf('comment','We have shown only seven bottles from the full vertical profile.\n');
        cprintf('text',' \n')
        pause(6)
        cprintf('comment','We now convert Practical Salinity, SP, into Absolute Salinity, SA,\n');
        cprintf('comment','using the function "gsw_SA_from_SP",\n');
        pause(6)
        cprintf('text','SA = gsw_SA_from_SP(SP,p,long,lat)\n');
        cprintf('text','%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %s \n' ,'SA = [',SA_chck_cast([1,22,29:4:45],1)',']');
        cprintf('text',' \n')
        pause(6)
        cprintf('comment','We now convert in situ temperature, t, into Conservative Temperature, CT\n');
        cprintf('comment','using the function "gsw_CT_from_t",\n');
        pause(6)
        cprintf('text','CT = gsw_CT_from_t(SA,t,p)\n');
        cprintf('text','%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'CT = [',CT_chck_cast([1,22,29:4:45],1)',']');
        cprintf('text',' \n')
        pause(6)
        cprintf('keywords','We now plot the profile on the Absolute Salinity - Conservative Temperature diagram\n');
        cprintf('comment','Potential density anomaly contours are shown for two different reference\n'); 
        cprintf('comment','pressures; 0 dbar and 2000 dbar. These values are obtained by using the function\n')
        cprintf('comment','"gsw_rho_CT" - 1000 kg/m3, as follows,\n')
        cprintf('text','sigma_Theta = gsw_rho_CT(SA,CT,0) - 1000\n');
        cprintf('text','sigma_2 = gsw_rho_CT(SA,CT,2000) - 1000\n');
        pause(8)
    catch
        disp('Well done! The gsw_check_values function confirms that the')
        disp('Gibbs SeaWater (GSW) Oceanographic Toolbox is installed correctly.');
        disp(' ')
        disp('Press enter to continue.')
        pause
        disp(' ')
        disp('A short demonstation of the GSW Oceanographic Toolbox now follows.');
        pause(3)
        disp('The following vertical profile, from the North Pacific, is of');
        disp('Practical Salinity, SP, and in situ temperature, t, as a function');
        disp('of pressure, p,');
        fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %s \n' ,'SP = [',gsw_cv.SP_chck_cast([1,22,29:4:45],1)',']');
        fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %s \n' ,'t  = [',gsw_cv.t_chck_cast([1,22,29:4:45],1)',']');
        fprintf(1,'%s %7.0f  %7.0f  %7.0f  %7.0f  %7.0f  %7.0f  %7.0f %s \n' ,'p  = [',gsw_cv.p_chck_cast([1,22,29:4:45],1)',']');
        disp('We have shown only seven bottles from the full vertical profile.');
        disp(' ')
        pause(6)
        disp('We now convert Practical Salinity, SP, into Absolute Salinity, SA,');
        disp('using the function "gsw_SA_from_SP",');
        disp('SA = gsw_SA_from_SP(SP,p,long,lat)');
        fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %s \n' ,'SA = [',SA_chck_cast([1,22,29:4:45],1)',']');
        disp(' ')
        pause(6)
        disp('We now convert in situ temperature, t, into Conservative Temperature, CT');
        disp('using the function "gsw_CT_from_t",');
        disp('CT = gsw_CT_from_t(SA,t,p)');
        fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %s \n' ,'CT = [',CT_chck_cast([1,22,29:4:45],1)',']');
        disp(' ')
        pause(6)
        disp('We now plot the profile on the Absolute Salinity - Conservative Temperature diagram');
        disp('Potential density anomaly contours are shown for two different reference'); 
        disp('pressures; 0 dbar and 2000 dbar. These values are obtained by using the function')
        disp('"gsw_rho_CT" - 1000 kg/m3, as follows,')
        disp('sigma_Theta = gsw_rho_CT(SA,CT,0) - 1000');
        disp('sigma_2 = gsw_rho_CT(SA,CT,2000) - 1000');
        pause(8)
    end
    pause(4)
    try
        SA_gridded = meshgrid([34:0.01:35.6],1:301);
        CT_gridded = meshgrid([0:0.1:30],1:161)';
        isopycs_0 = gsw_rho_CT(SA_gridded,CT_gridded,0)-1000;
        isopycs_2 = gsw_rho_CT(SA_gridded,CT_gridded,2000)-1000;
        figure
        [c0,h0]=contour(SA_gridded,CT_gridded,isopycs_0,[21:0.5:29],'k','linewidth',0.5);
        hold on
        [c2,h2]=contour(SA_gridded,CT_gridded,isopycs_2,[21:0.5:40],'r','linewidth',0.5);     
        clabel(c0,h0,'labelspacing',360);
        clabel(c2,h2,'labelspacing',240,'color','r');
        plot(SA_chck_cast(:,1),CT_chck_cast(:,1),'b','linewidth',2)
        plot(SA_chck_cast([1,22,29:4:45],1),CT_chck_cast([1,22,29:4:45],1),'co','linewidth',2,'Markersize',8)
        xlabel('Absolute Salinity, S_A [g/kg]','fontsize',15)
        ylabel('Conservative Temperature,  {\Theta} [{\circ}C]','fontsize',15)
        title(['S_A - {\Theta} for a profile at (',num2str(gsw_cv.long_chck_cast(1)),' {\circ}E, ',num2str(gsw_cv.lat_chck_cast(1)),' {\circ}N)'],'fontsize',15)
        lh = plot(1:2,1:2,'k',1:2,1:2,'r');
        lh2 = legend(lh,'{\sigma^\Theta}','{\sigma_2}','Location','southEast');
        set(lh2,'FontSize',15);
    catch
        disp('It appears that you are running MATLAB without the desktop or display,')
        disp('so we can not show the resulting figure.')
    end
    clear all
end

