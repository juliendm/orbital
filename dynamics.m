
%--------------------------------------%
% BEGIN: function dynamics             %
%--------------------------------------%

function [dr,dlon,dlat,dv,dgam,dal,dm,da,db,dt,pdyn,hr,nx,ny,nz,thu,cd,cl,Fdrag,Flift,rho,p,Tenv,mach,rey1m,trim_fwd,trim_aft,ka_fwd,ka_aft] = dynamics(in,param,auxdat,k)

  engine_type={'booster' 'spaceplane' 'spaceplane' 'spaceplane' 'spaceplane' 'spaceplane' 'spaceplane'};
  engine_state = {'on' 'off' 'on' 'off' 'off' 'off' 'off'};
  pha = {'a' 'a' 'a' 'r' 'r' 'r' 'r'};

  re = auxdat.re;
  gm = auxdat.gm;
  omega = auxdat.omega;
  j2 = auxdat.j2;
  f_ell = auxdat.f_ell;
  ispg0 = auxdat.ispg0;
  atm_model = auxdat.atm_model;
  date0_doy = auxdat.date0_doy;
  date0_sec = auxdat.date0_sec;

  tt = in.time;

  r = in.state(:,1);
  lon = in.state(:,2);
  lat = in.state(:,3);
  v = in.state(:,4);
  gam = in.state(:,5);
  al = in.state(:,6);
  m = in.state(:,7);
  aoa = in.state(:,8);
  bank = in.state(:,9);
  tva = in.state(:,10);

  dv1 = param(:,1);
  dv2 = param(:,2);
  dv3 = param(:,3);
  dv4 = param(:,4);

  daoa = in.control(:,1);
  dbank = in.control(:,2);
  dtva = in.control(:,3);

  n = length(r);

  [h,lon,glat] = latgeo(int32(n),double(re),double(f_ell),...
      double(r),double(lon),double(lat));
  h = reshape(h,n,1);
  lon = reshape(lon,n,1);
  glat = reshape(glat,n,1);


  if(strcmp(pha(k),'a')) ar_flag=1;end
  if(strcmp(pha(k),'r')) ar_flag=2;end

  % [cd,cl,rho,p] = get_cd_cl_rho_p(int32(n),double(h),double(lon),...
  %     double(glat),double(aoa*180/pi),double(v),double(tt),double(re),...
  %     int32(date0_doy),double(date0_sec),int32(atm_model),int32(ar_flag));


  % if(strcmp(engine_state(k),'on'))
  %   aoa = 0.0*r;   %%%%%%%%%%%%%%%%%%%%%%%%%%% WHEN PROPU IS ON: AOA = TVA, SO NO AERODYNAMICS EFFECTS
  %   tva = control_angle;
  % else
  %   aoa = control_angle;
  %   tva = 0.0*r;
  % end




  [cd,cl,rho,p,Tenv,mach,rey1m,trim_fwd,trim_aft,ka_fwd,ka_aft] = aerodynamics(k,n,h,lon,...
      glat,aoa*180/pi,v,tt,re,...
      date0_doy,date0_sec,atm_model,ar_flag,...
      dv1,dv2,dv3,dv4);

  cd = reshape(cd,n,1);
  cl = reshape(cl,n,1);
  rho = reshape(rho,n,1);
  p = reshape(p,n,1);
  Tenv = reshape(Tenv,n,1);
  mach = reshape(mach,n,1);
  rey1m = reshape(rey1m,n,1);
  trim_fwd = reshape(trim_fwd,n,1);
  trim_aft = reshape(trim_aft,n,1);
  ka_fwd = reshape(ka_fwd,n,1);
  ka_aft = reshape(ka_aft,n,1);


  if(strcmp(engine_type(k),'booster'))
     isp = auxdat.booster_isp;
     thu0 = param(:,end-3);
     % thu0 = auxdat.booster_thu0;
     sref = auxdat.booster_sref;
     na = auxdat.booster_na;
  elseif(strcmp(engine_type(k),'spaceplane'))
     isp = auxdat.spaceplane_isp;
     thu0 = param(:,end-1);
     % thu0 = auxdat.spaceplane_thu0;
     sref = auxdat.spaceplane_sref;
     na = auxdat.spaceplane_na;
  end

  if(strcmp(engine_state(k),'on'))
      thu = thu0-na*p;
      C = cos(aoa);
      S = sin(aoa);
      ca = -cl.*S+cd.*C;
      cn = cl.*C+cd.*S;
      ca = ca*0.85;
      cd = cn.*S+ca.*C;
      cl = cn.*C-ca.*S;
  else
      thu0 = 0.0;
      thu = zeros(size(p));
  end


  Fi = 0.0;
  Fj = 0.0;
  Fk = 0.0;

  j2coef = gm*re^2*j2;
  Fx = m.*(j2coef./r.^4).*cos(lat).*cos(lon).*(6.0-15.0/2.0.*cos(lat).^2);
  Fy = m.*(j2coef./r.^4).*cos(lat).*sin(lon).*(6.0-15.0/2.0.*cos(lat).^2);
  Fz = m.*(j2coef./r.^4).*sin(lat).*(3.0-15.0/2.0.*cos(lat).^2);
  Fi = Fi+(cos(lon).*cos(lat).*Fx+cos(lat).*sin(lon).*Fy+sin(lat).*Fz);
  Fj = Fj+(-sin(lon).*Fx+cos(lon).*Fy);
  Fk = Fk+(-sin(lat).*cos(lon).*Fx-sin(lat).*sin(lon).*Fy+cos(lat).*Fz);

  dr = v.*sin(gam);
  dlon = v.*sin(al).*cos(gam)./r./cos(lat);
  dlat = v./r.*cos(al).*cos(gam);
  dv =  -gm./r.^2.*sin(gam) ...
        -rho.*v.^2*sref.*cd./m/2.0 ...
        +thu.*cos(aoa+tva)./m ...
        +omega^2*r.*cos(lat).*(sin(gam).*cos(lat)-cos(gam).*sin(lat).*cos(al)) ...
        +Fi.*sin(gam)./m+Fj.*sin(al).*cos(gam)./m+Fk.*cos(al).*cos(gam)./m;
  dgam = -gm./r.^2.*cos(gam)./v ...
         +rho.*v*sref.*cl.*cos(bank)./m/2.0 ...
         +thu.*sin(aoa+tva).*cos(bank)./m./v ...
         +v./r.*cos(gam) ...
         +2.0*omega*sin(al).*cos(lat) ...
         +omega^2*r.*cos(lat).*(cos(gam).*cos(lat)+sin(gam).*sin(lat).*cos(al))./v ...
         +Fi./m./v.*cos(gam)-Fj./m./v.*sin(al).*sin(gam)-Fk./m./v.*sin(al).*sin(gam);
  dal = rho.*v*sref.*cl.*sin(bank)./m./cos(gam)/2.0 ...
        +thu.*sin(aoa+tva).*sin(bank)./m./v./cos(gam) ...
        +v./r.*cos(gam).*sin(al).*sin(lat)./cos(lat) ...
        +2.0*omega*(sin(lat)-cos(lat).*cos(al).*tan(gam)) ...
        +omega^2*r.*sin(lat).*cos(lat).*sin(al)./cos(gam)./v ...
        -Fk./m./v.*sin(al)./cos(gam)+Fj./m./v.*cos(al)./cos(gam);
  dm = r*0.0-thu0/isp/ispg0;

  da = daoa;
  db = dbank;
  dt = dtva;

  %
  % load factor in the body frame of the vehicle
  %
  Flift=cl.*rho.*v.^2*sref/2;
  Fdrag=cd.*rho.*v.^2*sref/2;

  z=tt*0;z=[z z z];


  ct=cos(tva);st=sin(tva);

  FT=z;
  FT(:,1)= ct.*thu;      % body frame 
  FT(:,3)=-st.*thu;      % body frame    % + or - ????????????

  FD=z;FD(:,1)=-Fdrag;   % aerodynamic frame
  FL=z;FL(:,3)=-Flift;   % aerodynamic frame

  ca=cos(aoa);sa=sin(aoa);

  x=ca.*FD(:,1)-sa.*FD(:,3);
  y=FD(:,2);
  z=sa.*FD(:,1)+ca.*FD(:,3);
  FD=[x y z];

  x=ca.*FL(:,1)-sa.*FL(:,3);
  y=FL(:,2);
  z=sa.*FL(:,1)+ca.*FL(:,3);
  FL=[x y z];
  % FD and Fl are now in body frame

  g0=gm/re^2;kn=length(tt);
  Fn=zeros(n,3);
  for k=1:kn
  Fn(k,1:3)=(FT(k,1:3)+FD(k,1:3)+FL(k,1:3))/(m(k)*g0);
  end

  % absolute value of the load factor in 
  % the body frame z axis

  nx = Fn(:,1);
  ny = Fn(:,2);
  nz = Fn(:,3);

  %nz = abs(Fn(:,3));
  %nz = sqrt(Fn(:,1).^2+Fn(:,2).^2+Fn(:,3).^2);

  [pdyn,hr] = get_pdyn_hr(in,auxdat,'soar');

end

%------------------------------------%
% END: function dynamics             %
%------------------------------------%

