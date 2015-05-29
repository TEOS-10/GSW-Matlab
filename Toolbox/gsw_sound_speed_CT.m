function sound_speed = gsw_sound_speed_CT(SA,CT,p)

% gsw_sound_speed                            sound speed (48-term equation)
%                                (approximate with a r.m.s. of 6.7 cm s^-1)
%==========================================================================
% This function has changed name to "gsw_sound_speed"
% Note that this is the final version of GSW (Version 3.04) that will 
% contain the function "gsw_sound_speed_CT"
%==========================================================================

warning('This function has changed name to "gsw_sound_speed".  This is the final version of GSW that will support "gsw_sound_speed_CT" ')

sound_speed = gsw_sound_speed(SA,CT,p);

end