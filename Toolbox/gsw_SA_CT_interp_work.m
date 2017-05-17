
 load SA_CT_interp_data_2
 
 SA_rr = NaN(length(p_i),length(lon_raw));
 CT_rr = SA_rr;
 SA_linear = SA_rr;
 CT_linear = SA_rr;
 SA_i = SA_rr;
 CT_i = SA_rr;

 for I = 1:length(lon_raw)
     [Idata] = find(~isnan(SA_SAonly(:,I) + CT_SAonly(:,I) + p_raw(:,I)));
     if length(Idata) > 3
         [SA_rr(:,I),CT_rr(:,I)] = gsw_rr68_interp_SA_CT(SA_SAonly(Idata,I),CT_SAonly(Idata,I),p_raw(Idata,I),p_i);
         [SA_linear(:,I),CT_linear(:,I)] = gsw_linear_interp_SA_CT(SA_SAonly(Idata,I),CT_SAonly(Idata,I),p_raw(Idata,I),p_i);
         [SA_i(:,I), CT_i(:,I)] = gsw_SA_CT_interp(SA_SAonly(Idata,I),CT_SAonly(Idata,I),p_raw(Idata,I),p_i);
         
     end
 end
 
 
  
%---------------------------------------------------------------------------

p_all = unique(sort([p; p_i]));
[I1,I2,I3] = intersect(p_i,p_all);
intergral_N2_all = interp1q(p,intergral_N2,p_all);
[SA_i_all, CT_i_all] = gsw_rr68_interp_SA_CT(SA,CT,intergral_N2,intergral_N2_all);
SA_i = SA_i_all(I2);
CT_i = CT_i_all(I3);



figure,plot(SA,CT,'o',SA_i, CT_i)
figure, plot(SA,p,'o', SA_i,p_i), set(gca,'ydir','reverse')
figure, plot(CT,p,'o', CT_i,p_i), set(gca,'ydir','reverse')

figure, plot(SA,p,'o', SA_i_all,p_all), set(gca,'ydir','reverse')
figure, plot(CT,p,'o', CT_i_all,p_all), set(gca,'ydir','reverse')


%---------------------------------------------------------------------------

dp = diff(p);
if any(dp > 60)  
    p50 = [min(p_i):50:max(p_i)];
    p50 = p50(:);
    p_plus50 = unique(sort([p; p50]));    
    intergral_N2_i_p_plus50 = interp1q(p,intergral_N2,p_plus50);
    [SA_i_p_plus50, CT_i_p_plus50] = gsw_rr68_interp_SA_CT(SA,CT,intergral_N2, intergral_N2_i_p_plus50);
    [SA_i50, CT_i50] = gsw_rr68_interp_SA_CT(SA_i_p_plus50,CT_i_p_plus50,p_plus50, p_i);
else
    [SA_i50, CT_i50] = gsw_rr68_interp_SA_CT(SA, CT,p,p_i);
end
figure,plot(SA,CT,'o',SA_i50, CT_i50)
figure, plot(SA,p,'o', SA_i50,p_i), set(gca,'ydir','reverse')
figure, plot(CT,p,'o', CT_i50,p_i), set(gca,'ydir','reverse')

%---------------------------------------------------------------------------

dp = diff(p);
if any(dp > 100)        
    p100 = [0:100:8000];
    p100 = p100(:);
    p_plus100 = unique(sort([p; p100]));    
   % intergral_N2_i_p_plus100 = interp1q(p,intergral_N2,p_plus100);    % linear
    intergral_N2_i_p_plus100 = interp1(p,intergral_N2,p_plus100,'pchip');
    [SA_i_p_plus100, CT_i_p_plus100] = gsw_rr68_interp_SA_CT(SA,CT,intergral_N2,intergral_N2_i_p_plus100);
    [SA_i100, CT_i100] = gsw_rr68_interp_SA_CT(SA_i_p_plus100, CT_i_p_plus100,p_plus100,p_i);
else
    [SA_i100, CT_i100] = gsw_rr68_interp_SA_CT(SA, CT,p,p_i);
end
figure,plot(SA,CT,'o',SA_i100, CT_i100)
figure, plot(SA,p,'o', SA_i100,p_i), set(gca,'ydir','reverse')
figure, plot(CT,p,'o', CT_i100,p_i), set(gca,'ydir','reverse')

%--------------------------------------------------------------------------

%     max_dp = 50;
%     dp = diff(p_tmp);
%     if any(dp > max_dp)
%         Ibg = find(dp > max_dp);
%         lbg = length(Ibg);
%         p_gt_max_dp = NaN(ceil(max(p_tmp)/max_dp),lbg);
%         for I = 1:lbg
%             dummy = [p_tmp(Ibg(I)): dp(Ibg(I))/ceil((p_tmp(Ibg(I)+1)-p_tmp(Ibg(I)))/max_dp)  :p_tmp(Ibg(I)+1)];
%             p_gt_max_dp(1:length(dummy),I) = dummy;
%         end
%         p_plus = unique(sort([p_tmp; p_gt_max_dp(~isnan(p_gt_max_dp))]));
%         intergral_N2_plus = interp1q(p_tmp,intergral_N2,p_plus);
%         [SA_plus, CT_plus] = gsw_rr68_interp_SA_CT(SA_tmp,CT_tmp,intergral_N2,intergral_N2_plus);
%         [SA_i(:,Iprofile), CT_i(:,Iprofile)] = gsw_rr68_interp_SA_CT(SA_plus,CT_plus,p_plus,p_i_tmp); 
%     else
%         [SA_i(:,Iprofile), CT_i(:,Iprofile)] = gsw_rr68_interp_SA_CT(SA_tmp,CT_tmp,p_tmp,p_i_tmp);
%     end
