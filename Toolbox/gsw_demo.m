gsw_data = 'gsw_data_v3_0.mat';
gsw_data_file = which(gsw_data);
load (gsw_data_file,'gsw_demo_data','version_number');
clear gsw_data gsw_data_file

%test if Java Virtual Machine is running
try
    JavaVirtMach = system_dependent('useJava','jvm');
catch
    %   assume no Java Virtual Machine
    JavaVirtMach = 0;
end

fprintf(1,'%s %s %s \n','Welcome the Gibbs Seawater (GSW) Oceanographic Toolbox ( Version ',version_number,').');
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
fprintf(1,'CT = gsw_CT_from_t(SA,t,p) \n');
gsw_demo_data.CT = gsw_CT_from_t(gsw_demo_data.SA,gsw_demo_data.t,gsw_demo_data.p);
fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'CT = [',gsw_demo_data.CT([1,22,29:4:45],1)',']');
fprintf(1,' \n');
fprintf(1,'At this point the data has been converted into SA and CT, which are  \n');
fprintf(1,'the TEOS-10 salinity and temperature variables.  With these variables it \n');
fprintf(1,'is possible to compute the complete range of water column properties. \n');
fprintf(1,' \n');
pause(6)
fprintf(1,'The first property to be demonstrated is density (rho) as a function \n');
fprintf(1,'of SA and CT.  This is computed by using the function "gsw_rho". \n');
fprintf(1,'The use of a single algorithm for seawater density (the 75-term computationally \n');
fprintf(1,'efficient expression) ensures consistency between ocean modelling, observational \n');
fprintf(1,'oceanography, and  theoretical studies.  Note that this is not been the case to \n');
fprintf(1,'date under EOS-80. \n');
fprintf(1,'rho = gsw_rho(SA,CT,p) \n');
gsw_demo_data.rho = gsw_rho(gsw_demo_data.SA,gsw_demo_data.CT,gsw_demo_data.p);
fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'rho = [',gsw_demo_data.rho([1,22,29:4:45],1)',']');
fprintf(1,' \n');
pause(6)
fprintf(1,'Using this same programme, gsw_rho, it is possible to compute potential \n');
fprintf(1,'density by replacing the in-situ pressure, p with the reference pressure, \n');
fprintf(1,'p_ref. \n');
fprintf(1,' \n');
pause(2)
fprintf(1,'An example. We have set p_ref to be 2000 dbar, thus we have the potential \n');
fprintf(1,'density referenced to 2000 dbars. \n');
fprintf(1,'p_ref = 2000 \n');
fprintf(1,'pot_rho_2 = gsw_rho(SA,CT,p_ref) \n');
gsw_demo_data.pot_rho_2 = gsw_rho(gsw_demo_data.SA,gsw_demo_data.CT,gsw_demo_data.p_ref);
fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'pot_rho_2 = [',gsw_demo_data.pot_rho_2([1,22,29:4:45],1)',']');
fprintf(1,' \n');
pause(6)
fprintf(1,'The potential density anomaly can be obtained by using the function \n');
fprintf(1,'"gsw_rho" - 1000 kg/m^3. \n');
fprintf(1,'Two examples of this are sigma_0 and sigma_2 which can be calculated \n');
fprintf(1,'as follows \n');
fprintf(1,'sigma_0 = gsw_rho(SA,CT,0) - 1000 \n');
gsw_demo_data.sigma_0 = gsw_rho(gsw_demo_data.SA,gsw_demo_data.CT,0) -1000;
fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'sigma_0 = [',gsw_demo_data.sigma_0([1,22,29:4:45],1)',']');
fprintf(1,' \n');
pause(6)
fprintf(1,'sigma_2 = gsw_rho(SA,CT,2000) - 1000 \n');
gsw_demo_data.sigma_2 = gsw_rho(gsw_demo_data.SA,gsw_demo_data.CT,2000) - 1000;
fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'sigma_2 = [',gsw_demo_data.sigma_2([1,22,29:4:45],1)',']');
fprintf(1,' \n');
pause(6)
fprintf(1,'However, there are alternatives to the last two calls, we have provided \n');
fprintf(1,'some short-cuts for the standard oceaongraphic variables as functions of \n');
fprintf(1,'SA and CT, the alternative short-cuts to the above two calls are: \n');
fprintf(1,'sigma_0 = gsw_sigma0(SA,CT) \n');
fprintf(1,' and  \n');
fprintf(1,'sigma_2 = gsw_sigma2(SA,CT) \n');
fprintf(1,' \n');
pause(6)
fprintf(1,'Calculating the Conservative Temperature at which seawater freezes is \n');
fprintf(1,'done with the function \n');
fprintf(1,'"gsw_CT_freezing_poly" \n');
fprintf(1,'This programme allows the user to choose the amount of air which the water \n');
fprintf(1,'contains. When saturation_fraction is 0 the seawater contains no air, and \n'); 
fprintf(1,'when saturation_fraction is 1 the seawater is completely saturated with air. \n');
fprintf(1,'The default setting is to have the seawater air free. \n');
fprintf(1,'CT_freezing = gsw_CT_freezing_poly(SA,p) \n');
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
fprintf(1,'We have interpolated the raw observations to a 1 dbar regular spaced data using the \n');
fprintf(1,'function "gsw_SA_CT_interp". \n');
fprintf(1,'The raw data are the black dots and the interpolated data is the blue line \n');
fprintf(1,'Note that this plotting function relies on the functions \n');
fprintf(1,'"gsw_rho" and "gsw_CT_freezing_poly" \n');
fprintf(1,' \n');
fprintf(1,'p_i = [min(p):max(p)]; \n');
fprintf(1,'[SA_i, CT_i] = gsw_SA_CT_interp(SA,CT,p,p_i); \n');
fprintf(1,'p_ref = 2000 \n');
fprintf(1,'gsw_SA_CT_plot(SA,CT,p_ref,''\\itS\\rm_A - \\Theta plot'') \n');
fprintf(1,'hold on \n');
fprintf(1,'plot(SA_i,CT_i,''b-'') \n');
pause(6)
if JavaVirtMach == 1
    try
        fig1 = figure;
        set(fig1,'pos','default','menubar','none','numbertitle','off', ...
            'name','GSW demo')
        hax1 = axes('pos',[0 0 1 1],'parent',fig1);
        axis(hax1,'off');
        str1(1) = {'An \it{S}\rm_A - {\Theta} diagram will soon appear.'};
        str1(2) = {' '};
        str1(3) = {'When you have finished studying the plot,'};
        str1(4) = {'minimise or close the window '};
        str1(5) = {'and press enter to continue.'};
        text(.5,.6,str1,'parent',hax1,'horizontalalignment','center','fontsize',18)
        pause(6)
        clf(fig1)
        set(fig1,'pos','default','menubar','none','numbertitle','off', ...
            'name','GSW demo - An example SA_CT diagram')
        gsw_demo_data.p_i = [min(gsw_demo_data.p):max(gsw_demo_data.p)];
        [gsw_demo_data.SA_i,gsw_demo_data.CT_i] = gsw_SA_CT_interp(gsw_demo_data.SA,gsw_demo_data.CT,gsw_demo_data.p,gsw_demo_data.p_i);
        gsw_SA_CT_plot(gsw_demo_data.SA(:,1),gsw_demo_data.CT(:,1),gsw_demo_data.p_ref,[32.5:0.5:38],'\itS\rm_A - \Theta  diagram');
        hold on
        plot(gsw_demo_data.SA_i(:,1),gsw_demo_data.CT_i(:,1),'b-')
        fprintf(1,' \n');
        fprintf(1,'When you have finished studying the plot press enter to continue\n');
        pause
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
fprintf(1,'p_ref = 2000 \n');
fprintf(1,'geo_strf_dyn_height = gsw_geo_strf_dyn_height(SA,CT,p,p_ref) \n');
fprintf(1,'%s %7.4f  %7.4f  %7.4f  %7.4f  %7.4f  %7.4f %7.4f %s \n' ,'geo_strf_dyn_height  = [',gsw_demo_data.geo_strf_dyn_height([1,22,29:4:45],1)',']');
pause(4)
fprintf(1,'The end. \n');

clear gsw_demo_data version_number JavaVirtMach
