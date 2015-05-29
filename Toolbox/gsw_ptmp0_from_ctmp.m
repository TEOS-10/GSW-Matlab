function result = gsw_ptmp0_from_ctmp(SA,ct)

%%   
% result = gsw_ptmp0_from_ctmp(SA,ct)
%
% potential temperature of seawater from conservative temperature
%
% SA                  : Absolute Salinity                  [g/kg]
% ct                  : conservative temperature           [deg C]
%
% result              : potential temperature              [deg C]

%%

cp0 = 3991.86795711963;     % in concrete 02/12/08;

n0 = 0; n2 = 2;

if gsw_check_arrays(SA,ct)
    error('****    input array dimesnions in gsw_ptmp0_from_ctmp do not agree    ****')
end

s1 = SA*35./35.16504; ct1 = ct; p0 = zeros(size(ct));

a0 = -1.446013646344788d-2;     b0 =  1.000000000000000d+0;
a1 = -3.305308995852924d-3;     b1 =  6.506097115635800d-4;
a2 =  1.062415929128982d-4;     b2 =  3.830289486850898d-3;
a3 =  9.477566673794488d-1;     b3 =  1.247811760368034d-6;
a4 =  2.166591947736613d-3;
a5 =  3.828842955039902d-3;

a5ct = a5*ct1; b3ct = b3*ct1;

ct_factor = (a3+a4*s1+a5ct);

th0_num = a0+s1.*(a1+a2*s1)+ct1.*ct_factor;

rec_th0_den = ones(size(s1))./(b0+b1*s1+ct1.*(b2+b3ct));

th0 = th0_num.*rec_th0_den;

ct0 = gsw_ctmp(SA,th0); 

dth_dct = (ct_factor+a5ct-(b2+b3ct+b3ct).*th0).*rec_th0_den; 

theta = th0-(ct0-ct).*dth_dct; 

nloops = 1;                 % default

%%    NOTE: nloops = 1 gives theta with a maximum error of ... 

n = 0; 

while n<=nloops
    dct = gsw_ctmp(SA,theta)-ct;
    dct_dth = -(theta+273.15).*gsw_g(n0,n2,n0,SA,theta,p0)./cp0;
    theta = theta - dct./dct_dth;
    n = n+1;
end

result = theta;

end


