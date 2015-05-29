gsw_data = 'gsw_data_v3_0.mat';
gsw_data_file = which(gsw_data);
load (gsw_data_file,'gsw_demo_data');

clear gsw_data gsw_data_file

%test if Java Virtual Machine is running
try
    JavaVirtMach = system_dependent('useJava','jvm');
catch
%    assume no Java Virtual Machine
    JavaVirtMach = 0;
end

try
    cprintf('keywords','Welcome the Gibbs Seawater (GSW) Oceanographic Toolbox (version 3). \n');
    pause(3)
    cprintf('comment','This is a short demonstration of some of the features of the \n');
    cprintf('comment','GSW Oceanographic Toolbox. \n');
    cprintf('text',' \n');
    cprintf('keywords','The most important functions are the first two functions. \n');
    cprintf('text',' \n');
    cprintf('comment','The following vertical profiles, from the North Pacific, are of \n');
    cprintf('comment','Practical Salinity, SP, and in-situ temperature, t, as a function \n');
    cprintf('comment','of pressure, p, \n');
    pause(6)
    cprintf('text','%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %s \n' ,'SP = [',gsw_demo_data.SP([1,22,29:4:45],1)',']');
    cprintf('text','%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %s \n' ,'t  = [',gsw_demo_data.t([1,22,29:4:45],1)',']');
    cprintf('text','%s %7.0f  %7.0f  %7.0f  %7.0f  %7.0f  %7.0f  %7.0f %s \n' ,'p  = [',gsw_demo_data.p([1,22,29:4:45],1)',']');
    cprintf('comment','Note that, we have shown only seven bottles from the full vertical profile. \n');
    cprintf('text',' \n');
    pause(6)
    cprintf('keywords','The first step under TEOS-10 is to convert Practical Salinity, SP, \n');
    cprintf('keywords','into Absolute Salinity, SA. This is done with the function "gsw_SA_from_SP" \n');
    pause(6)
    cprintf('text','SA = gsw_SA_from_SP(SP,p,long,lat) \n');
    gsw_demo_data.SA = gsw_SA_from_SP(gsw_demo_data.SP,gsw_demo_data.p,gsw_demo_data.long,gsw_demo_data.lat);
    cprintf('text','%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %s \n' ,'SA = [',gsw_demo_data.SA([1,22,29:4:45],1)',']');
    cprintf('text',' \n');
    pause(6)
    cprintf('keywords','The second step is to convert in-situ temperature, t, into \n');
    cprintf('keywords','Conservative Temperature, CT, using the function \n');
    cprintf('keywords','"gsw_CT_from_t", \n');
    pause(6)
    cprintf('text','CT = gsw_CT_from_t(SA,t,p) \n');
    gsw_demo_data.CT = gsw_CT_from_t(gsw_demo_data.SA,gsw_demo_data.t,gsw_demo_data.p);
    cprintf('text','%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'CT = [',gsw_demo_data.CT([1,22,29:4:45],1)',']');
    cprintf('text',' \n');
    cprintf('comment','At this point the data has been converted into SA and CT, which are \n');
    cprintf('comment','the TEOS-10 salinity and temperature variables.  With these variables it \n');
    cprintf('comment','is possible to compute the complete range of water column properties. \n');
    cprintf('text',' \n');
    pause(6)
    cprintf('comment','The first property to be demonstrated is density (rho) as a function \n');
    cprintf('comment','of SA and CT.  This is computed by using the function "gsw_rho_CT". \n');
    cprintf('comment','The use of a single algorithm for seawater density (the 48-term computationally \n');
    cprintf('comment','efficient expression) ensures consistency between ocean modelling, observational \n');
    cprintf('comment','oceanography, and  theoretical studies.  Note that this is not been the case to \n');
    cprintf('comment','date under EOS-80. \n');
    cprintf('text','rho_CT = gsw_rho_CT(SA,CT,p) \n');
    gsw_demo_data.rho_CT = gsw_rho_CT(gsw_demo_data.SA,gsw_demo_data.CT,gsw_demo_data.p);
    cprintf('text','%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'rho_CT = [',gsw_demo_data.rho_CT([1,22,29:4:45],1)',']');
    cprintf('text',' \n');
    pause(6)
    cprintf('comment','Using this same programme, gsw_rho_CT, it is possible to compute potential \n');
    cprintf('comment','density by replacing the in-situ pressure, p with the reference pressure, \n');
    cprintf('comment','p_ref. \n');
    cprintf('text',' \n');
    pause(2)
    cprintf('comment','An example. We have set p_ref to be 2000 dbar, thus we have the potential \n');
    cprintf('comment','density referenced to 2000 dbars. \n');
    cprintf('text','pot_rho_CT_2 = gsw_rho_CT(SA,CT,p_ref) \n');
    gsw_demo_data.pot_rho_CT_2 = gsw_rho_CT(gsw_demo_data.SA,gsw_demo_data.CT,gsw_demo_data.p_ref);
    cprintf('text','%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'pot_rho_CT = [',gsw_demo_data.pot_rho_CT_2([1,22,29:4:45],1)',']');
    cprintf('text',' \n');
    pause(6)
    cprintf('comment','The potential density anomaly can be obtained by using the function \n');
    cprintf('comment','"gsw_rho_CT" - 1000 kg/m^3. \n');
    cprintf('comment','Two examples of this are sigma_0 and sigma_2 which can be calculated \n');
    cprintf('comment','as follows \n');
    cprintf('text','sigma_0 = gsw_rho_CT(SA,CT,0) - 1000 \n');
    gsw_demo_data.sigma_0 = gsw_rho_CT(gsw_demo_data.SA,gsw_demo_data.CT,0) -1000;
    cprintf('text','%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'sigma_0 = [',gsw_demo_data.sigma_0([1,22,29:4:45],1)',']');
    cprintf('text',' \n');
    pause(6)
    cprintf('text','sigma_2 = gsw_rho_CT(SA,CT,2000) - 1000 \n');
    gsw_demo_data.sigma_2 = gsw_rho_CT(gsw_demo_data.SA,gsw_demo_data.CT,2000) - 1000;
    cprintf('text','%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'sigma_2 = [',gsw_demo_data.sigma_2([1,22,29:4:45],1)',']');
    cprintf('text',' \n');
    pause(6)
    cprintf('comment','However, there are alternatives to the last two calls, we have provided  \n');
    cprintf('comment','some short-cuts for the standard oceaongraphic variables as functions of \n');
    cprintf('comment','SA and CT, the alternative short-cuts to the above two calls are: \n');
    cprintf('text','sigma_0 = gsw_sigma0_CT(SA,CT) \n');
    cprintf('comment',' and \n');
    cprintf('text','sigma_2 = gsw_sigma2_CT(SA,CT) \n');
    cprintf('text',' \n');
    pause(6)
    cprintf('comment','Calculating the Conservative Temperature at which seawater freezes is \n');
    cprintf('comment','done with the function \n');
    cprintf('text','"gsw_CT_freezing" \n');
    cprintf('comment','This programme allows the user to choose the amount of air which the water \n');
    cprintf('comment','contains, at zero the water is unsaturated and at 1 it is completely \n');
    cprintf('comment','saturated, we have opted to set the default saturation level at maximum \n');
    cprintf('text','CT_freezing = gsw_CT_freezing(SA,p) \n');
    gsw_demo_data.CT_freezing = gsw_CT_freezing(gsw_demo_data.SA,gsw_demo_data.p);
    cprintf('text','%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'CT_freezing = [',gsw_demo_data.CT_freezing([1,22,29:4:45],1)',']');
    cprintf('text',' \n');
    cprintf('comment','Press enter to continue. \n');
    pause
    cprintf('keywords','We now plot the profile on the Absolute Salinity - Conservative Temperature diagram\n');
    cprintf('comment','This can be done by calling "gsw_SA_CT_plot".  This function plots the \n');
    cprintf('comment','Absolute Salinity and Conservative Temperature profile data on a SA-CT diagram \n');
    cprintf('comment','with user definied potential density contours and the Conservative Temperature \n');
    cprintf('comment','freezing line at p of 0 dbar.  The potential density anomaly contours are \n');
    cprintf('comment','referenced to user supplied depth are also included.  In this example we have \n');
    cprintf('comment','set the reference pressure to be 2000 dbar. \n');
    cprintf('comment','note that this plotting function relies on the functions \n');
    cprintf('comment','"gsw_rho_CT" and "gsw_CT_freezing" \n');
    cprintf('text',' \n');
    cprintf('text','p_ref = 2000 \n');
    cprintf('text','gsw_SA_CT_plot(SA,CT,p_ref,''\\itS_A - \\Theta plot'') \n');
    pause(6)
    if JavaVirtMach == 1
        try
            gsw_SA_CT_plot(gsw_demo_data.SA,gsw_demo_data.CT,gsw_demo_data.p_ref,[33:0.2:38],'\itS\rm_A - \Theta plot')
            pause(6)
            close all
        catch
            fprintf(1,' \n');
            fprintf(1,'It appears that you are running MATLAB without the Java Virtual Machine, \n');
            fprintf(1,'so we can not show the resulting figure. \n');
        end
    else
        fprintf(1,' \n');
        fprintf(1,'It appears that you are running MATLAB without the Java Virtual Machine, \n');
        fprintf(1,'so we can not show the resulting figure. \n');
    end
    cprintf('text',' \n');
    pause(2)
    cprintf('comment','The bouyancy (Brunt Vasaila) frequency squared (N^2) at the mid point \n');
    cprintf('comment','pressure (p_mid) between the "bottles" can be obtained by using the \n');
    cprintf('comment','function "gsw_Nsquared" \n');
    cprintf('text','[N2, p_mid] = gsw_Nsquared(SA,CT,p) \n');
    [gsw_demo_data.N2, gsw_demo_data.p_mid] = gsw_Nsquared(gsw_demo_data.SA,gsw_demo_data.CT,gsw_demo_data.p);
    cprintf('text','%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %s \n' ,'N2 = [',1e5*gsw_demo_data.N2([1,22,29:4:44],1)','] (*1e-5)');
    cprintf('text',' \n');
    cprintf('text','%s %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f %s \n' ,'p_mid = [',gsw_demo_data.p_mid([1,22,29:4:44],1)',']');
    cprintf('text',' \n');
    pause(6)
    cprintf('comment','The dynamic height anomaly, commmonly shortened to "dynamic height", can be \n');
    cprintf('comment','calculated with the function "gsw_geo_strf_dyn_height".  In this function \n');
    cprintf('comment','the user defines the the reference pressure that they want the dymanic height \n');
    cprintf('comment','relative to.  In this example we set p_ref to be 2000 dbar. \n');
    gsw_demo_data.geo_strf_dyn_height = gsw_geo_strf_dyn_height(gsw_demo_data.SA,gsw_demo_data.CT,gsw_demo_data.p,gsw_demo_data.p_ref);
    cprintf('text','geo_strf_dyn_height = gsw_geo_strf_dyn_height(SA,CT,p,p_ref) \n');
    cprintf('text','%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'geo_strf_dyn_height  = [',gsw_demo_data.geo_strf_dyn_height([1,22,29:4:45],1)',']');
    pause(4)
    cprintf('comment','The end. \n');
catch
    fprintf(1,'Welcome the Gibbs Seawater (GSW) Oceanographic Toolbox (version 3). \n');
    pause(3)
    fprintf(1,'This is a short demonstration of some of the features of the \n');
    fprintf(1,'GSW Oceanographic toolbox. \n');
    fprintf(1,' \n');
    fprintf(1,'The most important functions are the first two functions. \n');
    fprintf(1,' \n');
    fprintf(1,'The following vertical profiles, from the North Pacific, are of \n');
    fprintf(1,'Practical Salinity, SP, and in-situ temperature, t, as a function \n');
    fprintf(1,'of pressure, p, \n');
    pause(6)
    fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %s \n' ,'SP = [',gsw_demo_data.SP([1,22,29:4:45],1)',']');
    fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %s \n' ,'t  = [',gsw_demo_data.t([1,22,29:4:45],1)',']');
    fprintf(1,'%s %7.0f  %7.0f  %7.0f  %7.0f  %7.0f  %7.0f  %7.0f %s \n' ,'p  = [',gsw_demo_data.p([1,22,29:4:45],1)',']');
    fprintf(1,'Note that, we have shown only seven bottles from the full vertical profile. \n');
    fprintf(1,' \n');
    pause(6)
    fprintf(1,'The first step under TEOS-10 is to convert Practical Salinity, SP, \n');
    fprintf(1,'into Absolute Salinity, SA. This is done with the function "gsw_SA_from_SP" \n');
    pause(6)
    fprintf(1,'SA = gsw_SA_from_SP(SP,p,long,lat) \n');
    gsw_demo_data.SA = gsw_SA_from_SP(gsw_demo_data.SP,gsw_demo_data.p,gsw_demo_data.long,gsw_demo_data.lat);
    fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %s \n' ,'SA = [',gsw_demo_data.SA([1,22,29:4:45],1)',']');
    fprintf(1,' \n');
    pause(6)
    fprintf(1,'The second step is to convert in-situ temperature, t, into \n');
    fprintf(1,'Conservative Temperature, CT, using the function \n');
    fprintf(1,'"gsw_CT_from_t", \n');
    pause(6)
    fprintf(1,'CT = gsw_CT_from_t(SA,t,p) ');
    gsw_demo_data.CT = gsw_CT_from_t(gsw_demo_data.SA,gsw_demo_data.t,gsw_demo_data.p);
    fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'CT = [',gsw_demo_data.CT([1,22,29:4:45],1)',']');
    fprintf(1,' \n');
    fprintf(1,'At this point the data has been converted into SA and CT, which are  \n');
    fprintf(1,'the TEOS-10 salinity and temperature variables.  With these variables it \n');
    fprintf(1,'is possible to compute the complete range of water column properties. \n');
    fprintf(1,' \n');
    pause(6)
    fprintf(1,'The first property to be demonstrated is density (rho) as a function \n');
    fprintf(1,'of SA and CT.  This is computed by using the function "gsw_rho_CT". \n');
    fprintf(1,'The use of a single algorithm for seawater density (the 48-term computationally \n');
    fprintf(1,'efficient expression) ensures consistency between ocean modelling, observational \n');
    fprintf(1,'oceanography, and  theoretical studies.  Note that this is not been the case to \n');
    fprintf(1,'date under EOS-80. \n');
    fprintf(1,'rho_CT = gsw_rho_CT(SA,CT,p) \n');
    gsw_demo_data.rho_CT = gsw_rho_CT(gsw_demo_data.SA,gsw_demo_data.CT,gsw_demo_data.p);
    fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'rho_CT = [',gsw_demo_data.rho_CT([1,22,29:4:45],1)',']');
    fprintf(1,' \n');
    pause(6)
    fprintf(1,'Using this same programme, gsw_rho_CT, it is possible to compute potential \n');
    fprintf(1,'density by replacing the in-situ pressure, p with the reference pressure, \n');
    fprintf(1,'p_ref. \n');
    fprintf(1,' \n');
    pause(2)
    fprintf(1,'An example. We have set p_ref to be 2000 dbar, thus we have the potential \n');
    fprintf(1,'density referenced to 2000 dbars. \n');
    fprintf(1,'pot_rho_CT_2 = gsw_rho_CT(SA,CT,p_ref) \n');
    gsw_demo_data.pot_rho_CT_2 = gsw_rho_CT(gsw_demo_data.SA,gsw_demo_data.CT,gsw_demo_data.p_ref);
    fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'pot_rho_CT = [',gsw_demo_data.pot_rho_CT_2([1,22,29:4:45],1)',']');
    fprintf(1,' \n');
    pause(6)
    fprintf(1,'The potential density anomaly can be obtained by using the function \n');
    fprintf(1,'"gsw_rho_CT" - 1000 kg/m^3. \n');
    fprintf(1,'Two examples of this are sigma_Theta and sigma_2 which can be calculated \n');
    fprintf(1,'as follows \n');
    fprintf(1,'sigma_0 = gsw_rho_CT(SA,CT,0) - 1000 \n');
    gsw_demo_data.sigma_0 = gsw_rho_CT(gsw_demo_data.SA,gsw_demo_data.CT,0) -1000;
    fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'sigma_0 = [',gsw_demo_data.sigma_0([1,22,29:4:45],1)',']');
    fprintf(1,' \n');
    pause(6)
    fprintf(1,'sigma_2 = gsw_rho_CT(SA,CT,2000) - 1000 \n');
    gsw_demo_data.sigma_2 = gsw_rho_CT(gsw_demo_data.SA,gsw_demo_data.CT,2000) - 1000;
    fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'sigma_2 = [',gsw_demo_data.sigma_2([1,22,29:4:45],1)',']');
    fprintf(1,' \n');
    pause(6)
    fprintf(1,'However, there are alternatives to the last two calls, we have provided \n');
    fprintf(1,'some short-cuts for the standard oceaongraphic variables as functions of \n');
    fprintf(1,'SA and CT, the alternative short-cuts to the above two calls are: \n');
    fprintf(1,'sigma_0 = gsw_sigma0_CT(SA,CT) \n');
    fprintf(1,' and  \n');
    fprintf(1,'sigma_2 = gsw_sigma2_CT(SA,CT) \n');
    fprintf(1,' \n');
    pause(6)
    fprintf(1,'Calculating the Conservative Temperature at which seawater freezes is \n');
    fprintf(1,'done with the function \n');
    fprintf(1,'"gsw_CT_freezing" \n');
    fprintf(1,'This programme allows the user to choose the amount of air which the water \n');
    fprintf(1,'contains, at zero the water is unsaturated and at 1 it is completely \n');
    fprintf(1,'saturated, we have opted to set the default saturation level at maximum \n');
    fprintf(1,'CT_freezing = gsw_CT_freezing(SA,p) \n');
    gsw_demo_data.CT_freezing = gsw_CT_freezing(gsw_demo_data.SA,gsw_demo_data.p);
    fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'CT_freezing = [',gsw_demo_data.CT_freezing([1,22,29:4:45],1)',']');
    fprintf(1,' \n');
    fprintf(1,'Press enter to continue. \n');
    pause
    fprintf(1,'We now plot the profile on the Absolute Salinity - Conservative Temperature diagram \n');
    fprintf(1,'This can be done by calling "gsw_SA_CT_plot".  This function plots the \n');
    fprintf(1,'Absolute Salinity and Conservative Temperature profile data on a SA-CT diagram \n');
    fprintf(1,'with user definied potential density contours and the Conservative Temperature \n');
    fprintf(1,'freezing line at p of 0 dbar.  The potential density anomaly contours are \n');
    fprintf(1,'referenced to user supplied depth are also included.  In this example we have \n');
    fprintf(1,'set the reference pressure to be 2000 dbar. \n');
    fprintf(1,'note that this plotting function relies on the functions \n');
    fprintf(1,'"gsw_rho_CT" and "gsw_CT_freezing" \n');
    fprintf(1,' \n');
    fprintf(1,'p_ref = 2000 \n');
    fprintf(1,'gsw_SA_CT_plot(SA,CT,p_ref,''\\itS_A - \\Theta plot'') \n');
    pause(6)
    if JavaVirtMach == 1
        try
            gsw_SA_CT_plot(gsw_demo_data.SA,gsw_demo_data.CT,gsw_demo_data.p_ref,[33:0.2:38],'\itS\rm_A - \Theta plot')
            pause(6)
            close all
        catch
            fprintf(1,' \n');
            fprintf(1,'It appears that you are running MATLAB without the Java Virtual Machine, \n');
            fprintf(1,'so we can not show the resulting figure. \n');
        end
    else
        fprintf(1,' \n');
        fprintf(1,'It appears that you are running MATLAB without the Java Virtual Machine, \n');
        fprintf(1,'so we can not show the resulting figure. \n');
    end        
    fprintf(1,' \n');
    pause(2)
    fprintf(1,'The bouyancy (Brunt Vasaila) frequency squared (N^2) at the mid point \n');
    fprintf(1,'pressure (p_mid) between the "bottles" can be obtained by using the \n');
    fprintf(1,'function "gsw_Nsquared" \n');
    fprintf(1,'[N2, p_mid] = gsw_Nsquared(SA,CT,p) \n');
    [gsw_demo_data.N2, gsw_demo_data.p_mid] = gsw_Nsquared(gsw_demo_data.SA,gsw_demo_data.CT,gsw_demo_data.p);
    fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %s \n' ,'N2 = [',1e5*gsw_demo_data.N2([1,22,29:4:44],1)','] (*1e-5)');
    fprintf(1,' \n');
    fprintf(1,'%s %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f %s \n' ,'p_mid = [',gsw_demo_data.p_mid([1,22,29:4:44],1)',']');
    fprintf(1,' \n');
    pause(6)
    fprintf(1,'The dynamic height anomaly, commmonly shortened to "dynamic height", can be \n');
    fprintf(1,'calculated with the function "gsw_geo_strf_dyn_height".  In this function \n');
    fprintf(1,'the user defines the the reference pressure that they want the dymanic height \n');
    fprintf(1,'relative to. In this example we set p_ref to be 2000 dbar. \n');
    gsw_demo_data.geo_strf_dyn_height = gsw_geo_strf_dyn_height(gsw_demo_data.SA,gsw_demo_data.CT,gsw_demo_data.p,gsw_demo_data.p_ref);
    fprintf(1,'geo_strf_dyn_height = gsw_geo_strf_dyn_height(SA,CT,p,p_ref) \n');
    fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'geo_strf_dyn_height  = [',gsw_demo_data.geo_strf_dyn_height([1,22,29:4:45],1)',']');
    pause(4)
    fprintf(1,'The end. \n');
end

clear gsw_demo_data JavaVirtMach
