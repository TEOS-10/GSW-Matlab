function beta_CT = gsw_beta_CT(SA,CT,p)

% gsw_beta_CT                    saline contraction coefficient at constant
%                               Conservative Temperature (48-term equation)
%==========================================================================
% This function has changed name to "gsw_beta".
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_beta_CT"
%==========================================================================

warning('This function has changed name to "gsw_beta".  This is the final version of GSW that will support "gsw_beta_CT"')

beta_CT = gsw_beta(SA,CT,p);

end
