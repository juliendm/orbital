%-----------------------------------------------%
% Begin Function:  get_pdyn_hr                  %
%-----------------------------------------------%
function [pdyn,hr] = get_pdyn_hr(in,auxdat,name)

re = auxdat.re;
f_ell = auxdat.f_ell;
date0_doy = auxdat.date0_doy;
date0_sec = auxdat.date0_sec;
atm_model = auxdat.atm_model;
rho0 = auxdat.rho0;

r = in.state(:,1);
lon = in.state(:,2);
lat = in.state(:,3);
v = in.state(:,4);
aoa = in.state(:,8);

tt = in.time;
n = length(r);

[h,lon,glat] = latgeo(int32(n),double(re),double(f_ell),...
    double(r),double(lon),double(lat));

h = reshape(h,n,1);
lon = reshape(lon,n,1);
glat = reshape(glat,n,1);

[rho,p,ss] = get_rho_p_ss(int32(n),double(h),double(lon),double(glat),...
    double(tt),double(re),int32(date0_doy),double(date0_sec),int32(atm_model));

rho = reshape(rho,n,1);
p = reshape(p,n,1);
ss = reshape(ss,n,1);

if(strcmp(name,'soar'))
    rnose = get_rnose(int32(n),double(v./ss),double(aoa*180/pi));
elseif(strcmp(name,'stg'))
    rnose = ones(n,1);
end
rnose = reshape(rnose,n,1);

qc=110.3E6./sqrt(rnose).*sqrt(rho/rho0).*(v/7925E0).^(3.15);
qf=rho.*v.^3/2;
hr=qc.*(1-exp(-qf./qc));
hr(isnan(hr)) = 0;

pdyn = 0.5*rho.*v.^2;

end
%-----------------------------------------------%
% End Function:  get_pdyn_hr                    %
%-----------------------------------------------%
