
%-----------------------------------------------%
% Begin Function:  endpoint                     %
%-----------------------------------------------%

function output = s3toEndpoint(input)

	%
	% Variables at Start and End of each phase
	%
	Nph=length(input.phase);
	for k=1:Nph
		t0{k} = input.phase(k).initialtime;
		tf{k} = input.phase(k).finaltime;
		x0{k} = input.phase(k).initialstate;
		xf{k} = input.phase(k).finalstate;
	end

	s = input.parameter;

	dv_geo1 = s(1);
	dv_geo2 = s(2);
	dv_geo3 = s(3);
	dv_geo4 = s(4);

	spaceplane_dry_mass = compute_dry_mass(dv_geo1,dv_geo2,dv_geo3,dv_geo4);

	% booster_fuel_mass = s(end-2);
	% spaceplane_fuel_mass = s(end);
	% booster_fuel_mass = input.auxdata.booster_fuel_mass;
	% spaceplane_fuel_mass = input.auxdata.spaceplane_fuel_mass;

	%
	% Objective Function:
	%
	% output.objective = -xf{2}(4)/1000; % Max v_MECO

	%output.objective = xf{6}(7)/1000; % Min Dry_Mass

	%output.objective = x0{1}(7)/1000; % Min Mass
	%output.objective = xf{2}(7)/1000; % Min Mass
	
	%output.objective = dry_mass/1000; % Min Mass

	% output.objective = booster_fuel_mass + spaceplane_fuel_mass;

	output.objective = (s(end-3)+s(end-1))/1000;

	%
	% phase links
	% event groups 1-6
	%
	output.eventgroup(1).event = [tf{1}-t0{2},xf{1}-x0{2}];
	output.eventgroup(2).event = [tf{2}-t0{3},xf{2}-x0{3}];
	output.eventgroup(3).event = [tf{3}-t0{4},xf{3}-x0{4}];
	output.eventgroup(4).event = [tf{4}-t0{5},xf{4}-x0{5}];
	output.eventgroup(5).event = [tf{5}-t0{6},xf{5}-x0{6}];
	output.eventgroup(6).event = [tf{6}-t0{7},xf{6}-x0{7}];

	%
	% mass losses at separation
	%

	output.eventgroup(1).event(1+7) = xf{1}(7) - x0{2}(7) - input.auxdata.booster_dry_mass;

	%
	% phase durations
	% event group 7
	%
	output.eventgroup(7).event = [...
	tf{1}-t0{1},...
	tf{2}-t0{2},...
	tf{3}-t0{3},...
	tf{4}-t0{4},...
	tf{5}-t0{5},...
	tf{6}-t0{6},...
	tf{7}-t0{7}];

	% %
	% % upper stage injection conditions
	% % event group 7
	% %
	% r   = xf{2}(1);
	% lon = xf{2}(2);
	% lat = xf{2}(3);
	% v   = xf{2}(4);
	% gam = xf{2}(5);
	% al  = xf{2}(6);
	% [alt,lon,glat] = latgeo(int32(1),input.auxdata.re,input.auxdata.f_ell,r,lon,lat);
	% output.eventgroup(7).event = [gam, v, alt];

	% %
	% % regularization for re-entry
	% % event group 8
	% %
	% output.eventgroup(7).event = [input.phase(5).integral+input.phase(6).integral];

	%
	% final dry mass
	% event group 8
	%

	output.eventgroup(8).event = [xf{2}(7) - spaceplane_dry_mass - input.auxdata.us_mass]; % US IS STILL THERE AT LANDING !!!!!!!!!!!!!!!!!

	% output.eventgroup(9).event = [x0{1}(7) - spaceplane_dry_mass - input.auxdata.us_mass - booster_fuel_mass - spaceplane_fuel_mass - input.auxdata.booster_dry_mass]; % US IS STILL THERE AT LANDING !!!!!!!!!!!!!!!!!

	%
	% fuel limitations
	% event group 9
	%


	output.eventgroup(9).event = [x0{1}(7) - xf{1}(7), x0{2}(7) - xf{2}(7)];

	%
	% landing heading
	% event group 10
	%

	al  = xf{6}(6);

	output.eventgroup(10).event = [al]; %[al - 2*pi*floor(al/2/pi)];



end

%-----------------------------------------------%
% End Function:  endpoint                       %
%-----------------------------------------------%

function dry_mass = compute_dry_mass(dv_geo1,dv_geo2,dv_geo3,dv_geo4)

  % strake
  % x = 4.44779 y + 2.6669425 # - 5.11669
  a_strake = 4.44779;
  b_strake = 2.6669425;

  % lead
  % x = 0.551852 y + 8.672531 # 7.70679
  a_lead = 0.551852;
  b_lead_ini = 8.672531;
  b_lead = b_lead_ini + dv_geo1;

  % trail
  % x = -0.0549317 y + 13.984169525 # 14.0803 
  a_trail = -0.0549317;
  b_trail_ini = 13.984169525;
  b_trail = b_trail_ini + dv_geo2 + dv_geo4;

  wing_width_section_1_ini = 1.54;
  wing_width_section_2_ini = 2.41;
  wing_width_ini = wing_width_section_1_ini + wing_width_section_2_ini;
  wing_width = wing_width_ini + dv_geo3;
  wing_width_section_1 = (b_lead-b_strake)/(a_strake-a_lead);
  wing_width_section_2 = wing_width - wing_width_section_1;

  % Section 1
  % \int _ fuse_radius                          ^ (wing_width_section_1 + fuse_radius) [ (a_trail * y + b_trail) - (a_strake * y + b_strake)] dy

  x1 = 0.0;
  x2 = wing_width_section_1;
  a = a_trail;
  b = b_trail;
  c = a_strake;
  d = b_strake;

  area_1 = -1.0/2.0 * (x1 - x2) * (2.0*b - 2.0*d + (a - c) * (x1 + x2));

  % Section 2
  % \int _ (wing_width_section_1 + fuse_radius) ^ (wing_width + fuse_radius)           [ (a_trail * y + b_trail) - (a_lead * y + b_lead) ] dy

  x1 = wing_width_section_1;
  x2 = wing_width;
  a = a_trail;
  b = b_trail;
  c = a_lead;
  d = b_lead;

  area_2 = -1.0/2.0 * (x1 - x2) * (2.0*b - 2.0*d + (a - c) * (x1 + x2));

  % Mass

  area = area_1 + area_2;

  projected_area = 80.0 + 2.0*area;

%  dry_mass = 150.0 * projected_area - 10000.0;

  dry_mass = 90.0 * projected_area - 2900.0;

  % IS TOO HEAVY, WINGS DO NOT HAVE ENAUGH LIFT FOR LANDING PHASE
  % IS AT AOA 10 ALL THE TIME AND WING AREAS ARE GETTING AS BIG AS THEY CAN
  % REDUCING MAX SPEED ACHIEVABLE
  % dry_mass = 105.927 * projected_area - 697.025;

end



