function result = gsw_ptmp(SA,t,p,pr)

%%   
% result = gsw_ptmp(SA,t,p,pr)
%
% potential temperature of seawater from specific entropy 
%
% SA                  : Absolute Salinity                  [g/kg]
% t                   : temperature                        [deg C]
% p                   : sea (gauge) pressure               [dbar]
% pr                  : reference (gauge) pressure         [dbar]
%
% result              : potential temperature              [deg C]

%%

if gsw_check_arrays(SA,t,p,pr)
    error('****    input array dimesnions in gsw_ptmp do not agree    ****')
end


n0 = 0; n2 = 2;

s1 = SA*35./35.16504; t1 = t; p1 = p; pr1 = pr;

theta = t1+(p1-pr1).*( 8.65483913395442d-6  - ...
               s1  .*  1.41636299744881d-6  - ...
          (p1+pr1) .*  7.38286467135737d-9  + ...
               t1  .*(-8.38241357039698d-6  + ...
               s1  .*  2.83933368585534d-8  + ...
               t1  .*  1.77803965218656d-8  + ...
          (p1+pr1) .*  1.71155619208233d-10));

de_dt = 13.6;

nloops = 1;             % default

%    NOTE: nloops = 1 gives theta with a maximum error of ...
%          nloops = 2 gives theta with a maximum error of ...

n = 0; 

while n<=nloops
    n = n+1; theta_old = theta;
    dentropy = gsw_entropy(SA,theta,pr) - gsw_entropy(SA,t,p);
    theta = theta - dentropy./de_dt;
    theta = 0.5d0.*(theta+theta_old); 
    dentropy_dt = -gsw_g(n0,n2,n0,SA,theta,pr);
    theta = theta_old - dentropy./dentropy_dt;
end

result = theta;

end
