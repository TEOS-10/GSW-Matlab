% GSW Oceanographic Toolbox 
% Version 3.03 (R2011a) 23-May-2013
%
% documentation set
%  gsw_front_page              - front page to the GSW Oceanographic Toolbox
%  gsw_check_functions         - checks that all the GSW functions work correctly
%  gsw_demo                    - demonstrates many GSW functions and features
%
% Practical Salinity (SP), PSS-78  
%  gsw_SP_from_C               - Practical Salinity from conductivity, C (inc. for SP < 2)
%  gsw_C_from_SP               - conductivity, C, from Practical Salinity (inc. for SP < 2)
%  gsw_SP_from_R               - Practical Salinity from conductivity ratio, R (inc. for SP < 2)
%  gsw_R_from_SP               - conductivity ratio, R, from Practical Salinity (inc. for SP < 2)
%  gsw_SP_salinometer          - Practical Salinity from a laboratory salinometer (inc. for SP < 2)
%  gsw_SP_from_SK              - Practical Salinity from Knudsen Salinity
%
% Absolute Salinity (SA), Preformed Salinity (Sstar) and Conservative Temperature (CT) 
%  gsw_SA_from_SP              - Absolute Salinity from Practical Salinity
%  gsw_Sstar_from_SP           - Preformed Salinity from Practical Salinity
%  gsw_CT_from_t               - Conservative Temperature from in-situ temperature
%
% Absolute Salinity - Conservative Temperature plotting function 
%  gsw_SA_CT_plot              - function to plot Absolute Salinity - Conservative Temperature
%                                profiles on the SA-CT diagram, including the freezing line
%                                and selected potential density contours
%
% other conversions between temperatures, salinities, pressure and height
%  gsw_deltaSA_from_SP               - Absolute Salinity Anomaly from Practical Salinity
%  gsw_SA_Sstar_from_SP              - Absolute Salinity & Preformed Salinity from Practical Salinity
%  gsw_SR_from_SP                    - Reference Salinity from Practical Salinity
%  gsw_SP_from_SR                    - Practical Salinity from Reference Salinity
%  gsw_SP_from_SA                    - Practical Salinity from Absolute Salinity
%  gsw_Sstar_from_SA                 - Preformed Salinity from Absolute Salinity
%  gsw_SA_from_Sstar                 - Absolute Salinity from Preformed Salinity
%  gsw_SP_from_Sstar                 - Practical Salinity from Preformed Salinity
%  gsw_pt_from_CT                    - potential temperature from Conservative Temperature
%  gsw_t_from_CT                     - in-situ temperature from Conservative Temperature
%  gsw_CT_from_pt                    - Conservative Temperature from potential temperature
%  gsw_pot_enthalpy_from_pt          - potential enthalpy from potential temperature
%  gsw_pt_from_t                     - potential temperature
%  gsw_pt0_from_t                    - potential temperature with a reference pressure of zero dbar
%  gsw_t_from_pt0                    - in-situ temperature from potential temperature with p_ref = 0 dbar
%  gsw_t90_from_t48                  - ITS-90 temperature from IPTS-48 temperature
%  gsw_t90_from_t68                  - ITS-90 temperature from IPTS-68 temperature
%  gsw_z_from_p                      - height from pressure
%  gsw_p_from_z                      - pressure from height
%  gsw_z_from_depth                  - height from depth
%  gsw_depth_from_z                  - depth from height
%  gsw_Abs_Pressure_from_p           - Absolute Pressure,P, from pressure, p
%  gsw_p_from_Abs_Pressure           - pressure, p, from Absolute Pressure, P
%  gsw_entropy_from_CT               - entropy from Conservative Temperature
%  gsw_CT_from_entropy               - Conservative Temperature from entropy 
%  gsw_entropy_from_pt               - entropy from potential temperature
%  gsw_pt_from_entropy               - potential temperature from entropy
%  gsw_entropy_from_t                - entropy from in-situ temperature
%  gsw_t_from_entropy                - in-situ temperature from entropy
%  gsw_adiabatic_lapse_rate_from_CT  - adiabatic lapse rate from Conservative Temperature
%  gsw_adiabatic_lapse_rate_from_t   - adiabatic lapse rate from in-situ temperature
%  gsw_molality_from_SA              - molality of seawater
%  gsw_ionic_strength_from_SA        - ionic strength of seawater
%
% density and enthalpy, based on the 48-term expression for density 
%  gsw_rho                         - in-situ density and potential density
%  gsw_alpha                       - thermal expansion coefficient with respect to CT
%  gsw_beta                        - saline contraction coefficient at constant CT
%  gsw_rho_alpha_beta              - in-situ density, thermal expansion & saline contraction coefficients
%  gsw_alpha_on_beta               - alpha divied by beta
%  gsw_rho_first_derivaties        - first derivaties of rho
%  gsw_specvol                     - specific volume
%  gsw_specvol_anom                - specific volume anomaly
%  gsw_sigma0                      - sigma0 with reference pressure of 0 dbar
%  gsw_sigma1                      - sigma1 with reference pressure of 1000 dbar
%  gsw_sigma2                      - sigma2 with reference pressure of 2000 dbar
%  gsw_sigma3                      - sigma3 with reference pressure of 3000 dbar
%  gsw_sigma4                      - sigma4 with reference pressure of 4000 dbar
%  gsw_sound_speed                 - sound speed (approximate, with r.m.s. error of 0.067 m/s)
%  gsw_kappa                       - isentropic compressibility
%  gsw_cabbeling                   - cabbeling coefficient
%  gsw_thermobaric                 - thermobaric coefficient
%  gsw_SA_from_rho                 - Absolute Salinity from density
%  gsw_CT_from_rho                 - Conservative Temperature from density
%  gsw_CT_maxdensity               - Conservative Temperature of maximum density of seawater
%  gsw_internal_energy             - internal energy
%  gsw_enthalpy                    - enthalpy
%  gsw_enthalpy_diff               - difference of enthalpy between two pressures
%  gsw_dynamic_enthalpy            - dynamic enthalpy
%  gsw_enthalpy_first_derivaties   - first derivaties of enthalpy
%  gsw_enthalpy_second_derivaties  - second derivaties of enthalpy
%
% water column properties, based on the 48-term expression for density  
%  gsw_Nsquared                    - buoyancy (Brunt-Vaisala) frequency squared (N^2)
%  gsw_Turner_Rsubrho              - Turner angle & Rsubrho
%  gsw_IPV_vs_fNsquared_ratio      - ratio of the vertical gradient of potential density
%                                    (with reference pressure, p_ref), to the vertical 
%                                    gradient of locally-referenced potential density
%
% neutral properties, based on the 48-term expression for density
%  gsw_isopycnal_slope_ratio       - ratio of the slopes of isopycnals on the SA-CT diagram 
%                                    for p & p_ref
%  gsw_isopycnal_vs_ntp_CT_ratio   - ratio of the gradient of Conservative Temperature
%                                    in a potential density surface to that in the neutral 
%                                    tangent plane
%  gsw_ntp_pt_vs_CT_ratio          - ratio of gradients of potential temperature &
%                                    Conservative Temperature in a neutral tangent plane
%                                    (i.e. in a locally-referenced potential density surface)
%
% geostrophic streamfunctions, based on the 48-term expression for density
%  gsw_geo_strf_dyn_height          - dynamic height anomaly
%  gsw_geo_strf_dyn_height_pc       - dynamic height anomaly for piecewise constant profiles
%  gsw_geo_strf_isopycnal           - approximate isopycnal geostrophic streamfunction
%  gsw_geof_str_isopycnal_pc        - approximate isopycnal geostrophic streamfunction for
%                                     piecewise constant profiles
%  gsw_geo_strf_Montgomery          - Montgomery geostrophic streamfunction
%  gsw_geo_strf_Cunningham          - Cunningham geostrophic streamfunction
%
% geostrophic velocity 
%  gsw_geostrophic_velocity         - geostrophic velocity
%
% derivatives of entropy, CT and pt 
%  gsw_CT_first_derivatives         - first derivatives of Conservative Temperature
%  gsw_CT_second_derivatives        - second derivatives of Conservative Temperature
%  gsw_entropy_first_derivatives    - first derivatives of entropy
%  gsw_entropy_second_derivatives   - second derivatives of entropy
%  gsw_pt_first_derivatives         - first derivatives of potential temperature
%  gsw_pt_second_derivatives        - second derivatives of potential temperature
%
% freezing temperatures
%  gsw_CT_freezing         - Conservative Temperature freezing temperature of seawater
%  gsw_t_freezing          - in-situ freezing temperature of seawater
%  gsw_brineSA_CT          - Absolute Salinity of seawater at the freezing point (for given CT)
%  gsw_brineSA_t           - Absolute Salinity of seawater at the freezing point (for given t)
% 
% isobaric melting enthalpy and isobaric evaporation enthalpy
%  gsw_latentheat_melting  - latent heat of melting of ice into seawater (isobaric melting enthalpy)
%  gsw_latentheat_evap_CT  - latent heat of evaporation of water from seawater (isobaric
%                            evaporation enthalpy) with CT as input temperature
%  gsw_latentheat_evap_t   - latent heat of evaporation of water from seawater (isobaric
%                            evaporation enthalpy) with in-situ temperature, t, as input
%
% Planet Earth properties
%  gsw_f                   - Coriolis parameter
%  gsw_grav                - gravitational acceleration
%  gsw_distance            - spherical earth distance between points in the ocean
%
% steric height 
%  gsw_steric_height       - dynamic height anomaly divided by 9.7963 m s^-2
%
% TEOS-10 constants
%  gsw_T0                  - Celcius zero point; 273.15 K
%  gsw_P0                  - one standard atmosphere; 101 325 Pa
%  gsw_SS0                 - Standard Ocean Reference Salinity; 35.165 04 g/kg
%  gsw_uPS                 - unit conversion factor for salinities; (35.165 04/35) g/kg
%  gsw_cp0                 - the "specific heat" for use with CT; 3991.867 957 119 63 (J/kg)/K
%  gsw_C3515               - conductivity of SSW at SP=35, t_68=15, p=0; 42.9140 mS/cm
%  gsw_SonCl               - ratio of SP to Chlorinity; 1.80655 (g/kg)^-1
%  gsw_valence_factor      - valence factor of sea salt; 1.2452898
%  gsw_atomic_weight       - mole-weighted atomic weight of sea salt; 31.4038218... g/mol
%  
% density and enthalpy in terms of CT, based on the exact Gibbs function
%  gsw_rho_CT_exact                         - in-situ density and potential density
%  gsw_alpha_CT_exact                       - thermal expansion coefficient with respect to CT
%  gsw_beta_CT_exact                        - saline contraction coefficientat constant CT
%  gsw_rho_alpha_beta_CT_exact              - in-situ density, thermal expansion & saline contraction coefficient
%  gsw_alpha_on_beta_CT_exact               - alpha divied by beta
%  gsw_rho_first_derivaties_CT_exact        - first derivaties of rho 
%  gsw_specvol_CT_exact                     - specific volume
%  gsw_specvol_anom_CT_exact                - specific volume anomaly
%  gsw_sigma0_CT_exact                      - sigma0 with reference pressure of 0 dbar
%  gsw_sigma1_CT_exact                      - sigma1 with reference pressure of 1000 dbar
%  gsw_sigma2_CT_exact                      - sigma2 with reference pressure of 2000 dbar
%  gsw_sigma3_CT_exact                      - sigma3 with reference pressure of 3000 dbar
%  gsw_sigma4_CT_exact                      - sigma4 with reference pressure of 4000 dbar
%  gsw_sound_speed_CT_exact                 - sound speed
%  gsw_kappa_CT_exact                       - isentropic compressibility
%  gsw_cabbeling_CT_exact                   - cabbeling coefficient
%  gsw_thermobaric_CT_exact                 - thermobaric coefficient
%  gsw_SA_from_rho_CT_exact                 - Absolute Salinity from density
%  gsw_CT_from_rho_exact                    - Conservative Temperature from density
%  gsw_CT_maxdensity_exact                  - Conservative Temperature of maximum density of seawater
%  gsw_internal_energy_CT_exact             - internal energy
%  gsw_enthalpy_CT_exact                    - enthalpy
%  gsw_enthalpy_diff_CT_exact               - difference of enthalpy between two pressures
%  gsw_dynamic_enthalpy_CT_exact            - dynamic enthalpy
%  gsw_enthalpy_first_derivaties_CT_exact   - first derivaties of enthalpy
%  gsw_enthalpy_second_derivaties_CT_exact  - second derivaties of enthalpy
%
% Labroratory functions, for use with a densimeter measuremants
%  gsw_SA_from_rho_t_exact              - Absolute Salinity from density
%  gsw_deltaSA_from_rho_t_exact         - Absolute Salinity Anomaly from density
%  gsw_rho_t_exact                      - in-situ density
%
% basic thermodynamic properties in terms of in-situ t, based on the exact Gibbs function 
%  gsw_rho_t_exact                      - in-situ density
%  gsw_pot_rho_t_exact                  - potential density
%  gsw_sigma0_pt0_exact                 - sigma0 from pt0 with reference pressure of 0 dbar
%  gsw_alpha_wrt_CT_t_exact             - thermal expansion coefficient with respect to 
%                                         Conservative Temperature.
%  gsw_alpha_wrt_pt_t_exact             - thermal expansion coefficient with respect to 
%                                         potential temperature
%  gsw_alpha_wrt_t_exact                - thermal expansion coefficient with respect to 
%                                         in-situ temperature
%  gsw_beta_const_CT_t_exact            - saline contraction coefficient at constant 
%                                         Conservative Temperature.
%  gsw_beta_const_pt_t_exact            - saline contraction coefficient at constant 
%                                         potential temperature
%  gsw_beta_const_t_exact               - saline contraction coefficient at constant 
%                                         in-situ temperature
%  gsw_specvol_t_exact                  - specific volume
%  gsw_specvol_anom_t_exact             - specific volume anomaly
%  gsw_sound_speed_t_exact              - sound speed
%  gsw_kappa_t_exact                    - isentropic compressibility
%  gsw_kappa_const_t_exact              - isothermal compressibility
%  gsw_internal_energy_t_exact          - internal energy
%  gsw_enthalpy_t_exact                 - enthalpy
%  gsw_dynamic_enthalpy_t_exact         - dynamic enthalpy
%  gsw_SA_from_rho_t_exact              - Absolute Salinity from density
%  gsw_t_from_rho_exact                 - in-situ temperature from density
%  gsw_t_maxdensity_exact               - in-situ temperature of maximum density of seawater
%  gsw_cp_t_exact                       - isobaric heat capacity
%  gsw_isochoric_heat_cap_t_exact       - isochoric heat capacity
%  gsw_chem_potential_relative_t_exact  - relative chemical potential
%  gsw_chem_potential_water_t_exact     - chemical potential of water in seawater
%  gsw_chem_potential_salt_t_exact      - chemical potential of salt in seawater
%  gsw_Helmholtz_energy_t_exact         - Helmholtz energy
%  gsw_osmotic_coefficient_t_exact      - osmotic coefficient of seawater
%  gsw_osmotic_pressure_t_exact         - osmotic pressure of seawater
%
% Library functions of the GSW toolbox (internal functions; not intended to be called by users) 
%  (The GSW functions above call the following library functions.)
%  gsw_gibbs                 - the TEOS-10 Gibbs function and its derivatives
%  gsw_SAAR                  - Absolute Salinity Anomaly Ratio (excluding the Baltic Sea)
%  gsw_Fdelta                - ratio of Absolute to Preformed Salinity, minus 1
%  gsw_deltaSA_atlas         - Absolute Salinity Anomaly atlas value (excluding the Baltic Sea)
%  gsw_SA_from_SP_Baltic     - Calculates Absolute Salinity in the Baltic Sea
%  gsw_SP_from_SA_Baltic     - Calculates Practical Salinity in the Baltic Sea
%  gsw_infunnel              - "oceanographic funnel" check for the 48-term equation
%  gsw_entropy_part          - entropy minus the terms that are a function of only SA
%  gsw_entropy_part_zerop    - entropy_part evaluated at 0 dbar
%  gsw_interp_ref_cast       - linearly interpolates the reference cast
%  gsw_interp_SA_CT          - linearly interpolates (SA,CT,p) to the desired p
%  gsw_gibbs_pt0_pt0         - gibbs(0,2,0,SA,t,0)
%  gsw_specvol_SSO_0_p       - specvol_CT(35.16504,0,p)
%  gsw_enthalpy_SSO_0_p      - enthalpy_CT(35.16504,0,p)
%  gsw_Hill_ratio_at_SP2     - Hill ratio at a Practical Salinity of 2
% 
%  The GSW data set.
%  gsw_data_v3_0        - contains
%                          (1) the global data set of Absolute Salinity Anomaly Ratio,
%                          (2) the global data set of Absolute Salinity Anomaly atlas,                                    
%                          (3) a reference cast (for the isopycnal streamfunction), 
%                          (4) two reference casts that are used by gsw_demo 
%                          (5) three vertical profiles of (SP, t, p) at known long & lat, plus
%                              the outputs of all the GSW functions for these 3 profiles, and
%                              the required accuracy of all these outputs.
%