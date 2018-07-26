
%--------------------------------------%
% BEGIN: function continuous           %
%--------------------------------------%

function phaseout = s3toContinuous(input)

  global load_cases

  %disp('Call Continuous ...')

  time_cumul = [];
  mach_cumul = [];
  rey_cumul = [];
  aoa_cumul = [];
  nx_cumul = [];
  nz_cumul = [];
  pdyn_cumul = [];
  thu_cumul = [];
  mass_cumul = [];

  Nph = length(input.phase);

  for k=1:Nph

    [dr,dlon,dlat,dv,dgam,dal,dm,da,db,dt,pdyn,hr,nx,ny,nz,thu,cd,cl,Fdrag,Flift,rho,p,Tenv,mach,rey1m,el_def,bf_def,trim_fwd,trim_aft,ka_fwd,ka_aft] = dynamics(input.phase(k),input.phase(k).parameter,input.auxdata,k);

    phaseout(k).dynamics = [dr,dlon,dlat,dv,dgam,dal,dm,da,db,dt];

    % phaseout(k).path = [pdyn,hr,abs(nz)];
    % phaseout(k).path = [pdyn,hr,abs(nz),trim_fwd,trim_aft];
    phaseout(k).path = [pdyn,hr,abs(nz),trim_fwd,trim_aft,ka_fwd,ka_aft];




    mass = input.phase(k).state(:,7);
    aoa = input.phase(k).state(:,8);

    time_cumul = [time_cumul;input.phase(k).time];
    mach_cumul = [mach_cumul;mach];
    rey_cumul = [rey_cumul;rey1m];
    aoa_cumul = [aoa_cumul;aoa];
    nx_cumul = [nx_cumul;nx];
    nz_cumul = [nz_cumul;nz];
    pdyn_cumul = [pdyn_cumul;pdyn];
    thu_cumul = [thu_cumul;thu];
    mass_cumul = [mass_cumul;mass];

  end

  % for k=4:6

  %   daoa=input.phase(k).control(:,1);
  %   dbnk=input.phase(k).control(:,2);

  %   phaseout(k).integrand = abs(daoa)+abs(dbnk);

  % end









  max_fuel_mass = input.phase(2).state(1,7) - mass_cumul(end);
  fuel_mass_cumul = min(max_fuel_mass*ones(size(mass_cumul)), mass_cumul - mass_cumul(end));

  % NEED TO ADD A MIN DIST IN BETWEEN PEAKS BECAUSE AT THE BEGINING WHEN SOLUTION ARE VERY UNSMOOTH, IT WOULD IDENTIFY TOO MUCH PEAKS

  % Find Peaks
  [pks_nx,locs_nx] = findpeaks(nx_cumul,'MinPeakHeight',1.0,'MinPeakProminence',1.0);
  [pks_nz,locs_nz] = findpeaks(-nz_cumul,'MinPeakHeight',1.0,'MinPeakProminence',1.0);
  [pks_pdyn,locs_pdyn] = findpeaks(pdyn_cumul,'MinPeakHeight',5.0e3,'MinPeakProminence',4.0e3);

  nx   = [nx_cumul(locs_nx);nx_cumul(locs_nz);nx_cumul(locs_pdyn)];
  nz   = [nz_cumul(locs_nx);nz_cumul(locs_nz);nz_cumul(locs_pdyn)];
  pdyn = [pdyn_cumul(locs_nx);pdyn_cumul(locs_nz);pdyn_cumul(locs_pdyn)];
  mach = [mach_cumul(locs_nx);mach_cumul(locs_nz);mach_cumul(locs_pdyn)];
  rey  = [rey_cumul(locs_nx);rey_cumul(locs_nz);rey_cumul(locs_pdyn)];
  aoa  = [aoa_cumul(locs_nx);aoa_cumul(locs_nz);aoa_cumul(locs_pdyn)];
  thu  = [thu_cumul(locs_nx);thu_cumul(locs_nz);thu_cumul(locs_pdyn)];
  fuel_mass  = [fuel_mass_cumul(locs_nx);fuel_mass_cumul(locs_nz);fuel_mass_cumul(locs_pdyn)];

  load_cases = [mach,rey,aoa,nx,nz,thu,pdyn,fuel_mass];



  % dv_geo1 = input.phase(Nph).parameter(end,1);
  % dv_geo2 = input.phase(Nph).parameter(end,2);
  % dv_geo3 = input.phase(Nph).parameter(end,3);
  % dv_geo4 = input.phase(Nph).parameter(end,4);
  % dv_geo5 = input.phase(Nph).parameter(end,5);
  % dv_geo6 = input.phase(Nph).parameter(end,6);



%  disp(['Dry mass: ',num2str(spaceplane_dry_mass_global)]);

  % % Visualize Peaks

  % figure
  % plot(time_cumul,nx_cumul,time_cumul(locs_nx),pks_nx,'or')
  % hold on;
  % plot(time_cumul,-nz_cumul,time_cumul(locs_nz),pks_nz,'or')
  % saveas(gcf,'load_factor.png')

  % figure
  % plot(time_cumul,pdyn_cumul,time_cumul(locs_pdyn),pks_pdyn,'or')
  % saveas(gcf,'pdyn.png')



end

%------------------------------------%
% END: function continuous           %
%------------------------------------%





