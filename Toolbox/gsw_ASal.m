function result = gsw_ASal(SP,p0,longs0,lats0)

%%
% result = gsw_ASal(SP,p0,longs0,lats0)
%
% converts Practical Salinity to Absolute Salinity 
%
% SP                  : Practical Salinity                 [psu]
% p0                  : sea (gauge) pressure               [dbar]
% longs0              : longitude                          [deg E]     
% lats0               : latitude                           [deg N]
%
% result              : Absolute Salinity                  [g/kg]

%%

if gsw_check_arrays(SP,longs0,lats0,p0)
    error('****    input array dimensions in gsw_ASal do not agree    ****')
end

inds = find(isfinite(SP)); result = nan*ones(size(SP));

result(inds) = (35.16504/35)*SP(inds) + gsw_delta_SA(p0(inds),longs0(inds),lats0(inds));

flag = 1;

result(inds) = adjust_Baltic(result(inds),SP(inds),longs0(inds),lats0(inds),flag);

end

function result = adjust_Baltic(SA,SP,longs,lats,flag)

%%
% result = adjust_Baltic(SA,SP,longs,lats)
%
% for the Baltic Sea, overwrite Absolute Salinity with a value
% computed analytically from Practical Salinity, or vice versa
%
% SA                  : Absolute Salinity                  [g/kg]
% SP                  : Practical Salinity                 [psu]
% longs               : longitude                          [deg E]     
% lats                : latitude                           [deg N]
% flag                : flag - 1 or 2 
%
% result              : Absolute Salinity                  [g/kg]
%                         when flag = 1
%                     : Practical Salinity                 [psu]
%                         when flag = 2

%%

xb1 = 12.6; xb2 = 7; xb3 = 26; xb1a = 45; xb3a = 26;

yb1 = 50; yb2 = 59; yb3 = 69;

inds = find(xb2<longs & longs<xb1a & yb1<lats & lats<yb3);

if flag==1
    result = SA;
elseif flag==2
    result = SP;
end 

if ~isempty(inds)
    
    xx_left = interp1([yb1,yb2,yb3],[xb1,xb2,xb3],lats(inds));
    
    xx_right = interp1([yb1,yb3],[xb1a,xb3a],lats(inds));
    
    inds1 = find(xx_left<=longs(inds) & longs(inds)<=xx_right);
    
    if ~isempty(inds)
      if flag==1
        result(inds(inds1)) = (35.04104d0/35.d0)*SP(inds(inds1)) + 0.124d0;
      elseif flag==2
        result(inds(inds1)) = (35.d0/35.04104d0)*(SA(inds(inds1)) - 0.124d0);
      end
    end
    
    result = reshape(result,size(longs));
    
end

end

    
        

    

