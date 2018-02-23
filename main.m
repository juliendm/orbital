
clear all;
close all;
% clc;

addpath('/home/juliendm/gpopsii/GPOPS-II/lib/gpopsRPMIntegration/gpopsSnoptRPMI/');
addpath('/home/juliendm/gpopsii/GPOPS-II/lib/gpopsRPMIntegration/gpopsIpoptRPMI/');
addpath('/home/juliendm/gpopsii/GPOPS-II/lib/gpopsRPMIntegration/');
addpath('/home/juliendm/gpopsii/GPOPS-II/lib/gpopsRPMDifferentiation/gpopsSnoptRPMD/');
addpath('/home/juliendm/gpopsii/GPOPS-II/lib/gpopsRPMDifferentiation/gpopsIpoptRPMD/');
addpath('/home/juliendm/gpopsii/GPOPS-II/lib/gpopsRPMDifferentiation/');
addpath('/home/juliendm/gpopsii/GPOPS-II/lib/gpopsFiniteDifference/');
addpath('/home/juliendm/gpopsii/GPOPS-II/lib/gpopsADiGator/');
addpath('/home/juliendm/gpopsii/GPOPS-II/lib/gpopsAutomaticScaling/');
addpath('/home/juliendm/gpopsii/GPOPS-II/lib/gpopsMeshRefinement/');
addpath('/home/juliendm/gpopsii/GPOPS-II/lib/gpopsCommon/');
addpath('/home/juliendm/gpopsii/GPOPS-II/gpopsUtilities/');
addpath('/home/juliendm/gpopsii/GPOPS-II/license/');
addpath('/home/juliendm/gpopsii/GPOPS-II/nlp/snopt/');

global auxdata ph N bounds guess mesh ...
       t_0 t_f t_min t_max r_0 r_f r_min r_max h_0 h_f h_min h_max ...
       lon_0 lon_f lon_min lon_max glat_0 glat_f glat_min glat_max ...
       lat_0 lat_f lat_min lat_max v_0 v_f v_min v_max gam_0 gam_f gam_min ...
       gam_max al_0 al_f al_min al_max m_0 m_f m_min m_max aoa_0 aoa_f ...
       aoa_min aoa_max bank_0 bank_f bank_min bank_max daoa_0 daoa_f ...
       daoa_min daoa_max dbank_0 dbank_f dbank_min dbank_max ...
       pdyn_min pdyn_max hr_min hr_max nz_min nz_max ...
       % dv1_0 dv1_f dv2_0 dv2_f dv3_0 dv3_f dv4_0 dv4_f ...
       % dv1_min dv1_max dv2_min dv2_max dv3_min dv3_max dv4_min dv4_max ...


addpath('./externals');
addpath('./utils');

%initialitaion (linterp)
% get_cd_cl_rho_p();
get_rnose();


% REF_LENGTH_MOMENT= 17.0
% REF_AREA= 60.0

% surfpack_load('lift_sub','externals/auxfiles/aero/build_points_lift_sub.dat',11,1,'externals/auxfiles/aero/model_lift_sub.sps');
% surfpack_load('lift_sup','externals/auxfiles/aero/build_points_lift_sup.dat',11,1,'externals/auxfiles/aero/model_lift_sup.sps');
% surfpack_load('drag_sub','externals/auxfiles/aero/build_points_drag_sub.dat',11,1,'externals/auxfiles/aero/model_drag_sub.sps');
% surfpack_load('drag_sup','externals/auxfiles/aero/build_points_drag_sup.dat',11,1,'externals/auxfiles/aero/model_drag_sup.sps');

% surfpack_load('trim_fwd','externals/auxfiles/perfo/build_points_TRIM_FWD.dat',7,1,'externals/auxfiles/perfo/model_TRIM_FWD.sps');
% surfpack_load('trim_aft','externals/auxfiles/perfo/build_points_TRIM_AFT.dat',7,1,'externals/auxfiles/perfo/model_TRIM_AFT.sps');
% surfpack_load('k_alpha_fwd','externals/auxfiles/perfo/build_points_K_ALPHA_FWD.dat',7,1,'externals/auxfiles/perfo/model_K_ALPHA_FWD.sps');
% surfpack_load('k_alpha_aft','externals/auxfiles/perfo/build_points_K_ALPHA_AFT.dat',7,1,'externals/auxfiles/perfo/model_K_ALPHA_AFT.sps');




surfpack_load_sps('lift_sub','externals/auxfiles/aero/model_lift_sub.sps');
surfpack_load_sps('lift_sup','externals/auxfiles/aero/model_lift_sup.sps');
surfpack_load_sps('drag_sub','externals/auxfiles/aero/model_drag_sub.sps');
surfpack_load_sps('drag_sup','externals/auxfiles/aero/model_drag_sup.sps');

surfpack_load_sps('trim_fwd','externals/auxfiles/perfo/model_TRIM_FWD.sps');
surfpack_load_sps('trim_aft','externals/auxfiles/perfo/model_TRIM_AFT.sps');
surfpack_load_sps('k_alpha_fwd','externals/auxfiles/perfo/model_K_ALPHA_FWD.sps');
surfpack_load_sps('k_alpha_aft','externals/auxfiles/perfo/model_K_ALPHA_AFT.sps');



d2r = pi/180;
r2d = 180/pi;

auxdata.re = 6378137.0;
auxdata.gm = 3986004.418e+08;
auxdata.omega = 7292115e-11;
auxdata.j2 = 1.0826298213133048e-03;
auxdata.f_ell = 3.3528106647474805e-03;
auxdata.ispg0 = 9.80665;

auxdata.atm_model = 1;
auxdata.rho0 = 1.3;
auxdata.date0 = '2010 MAY 01 11:00:00.000';
[auxdata.date0_sec,auxdata.date0_doy] = utc2doydec(auxdata.date0);



%
% Data
%

auxdata.booster_dry_mass = 25000.0;
auxdata.booster_isp = 330.0;
auxdata.booster_sref = 120.0;
auxdata.booster_na = 2.55;

auxdata.spaceplane_isp = 375.0;
auxdata.spaceplane_sref = 120.0;
auxdata.spaceplane_na = 2.55;

auxdata.us_mass = 6500.0;

%%
%% phase definitions:
%% setup of bounds,guesses,meshing
%%

%
% Configuration
%

% auxdata.fuel_mass = 25000.0;

% dv1_min = -0.5   ; dv1_max = 0.5;
% dv2_min = -0.5   ; dv2_max = 0.5;
% dv3_min = -0.5   ; dv3_max = 0.5;
% dv4_min = -0.2   ; dv4_max = 0.5;



spaceplane_dry_mass_min = 5000.0;
spaceplane_dry_mass_max = 20000.0;

booster_thu0_min = 1200000.0; booster_thu0_max = 3000000.0;
booster_fuel_mass_min = 30000.0; booster_fuel_mass_max = 300000.0; 

spaceplane_thu0_min = 390000.0; spaceplane_thu0_max = 1000000.0;
spaceplane_fuel_mass_min = 50000.0; spaceplane_fuel_mass_max = 50000.0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%% VOLUME INSIDE SPACEPLANE !!!!!!!!!!!!!!!!!


dv1_min = -0.5   ; dv1_max = 0.5;
dv2_min = -0.5   ; dv2_max = 0.5;
dv3_min = -0.5   ; dv3_max = 0.5;
dv4_min = -0.2   ; dv4_max = 0.5;


bounds.parameter.lower = [dv1_min,dv2_min,dv3_min,dv4_min,booster_thu0_min,booster_fuel_mass_min,spaceplane_thu0_min,spaceplane_fuel_mass_min];
bounds.parameter.upper = [dv1_max,dv2_max,dv3_max,dv4_max,booster_thu0_max,booster_fuel_mass_max,spaceplane_thu0_max,spaceplane_fuel_mass_max];


dv1_0 = 0.0       ; dv1_f = 0.0;
dv2_0 = 0.0       ; dv2_f = 0.0;
dv3_0 = 0.0       ; dv3_f = 0.0;
dv4_0 = 0.0       ; dv4_f = 0.0;





%
% VEHICLE ascent
%
t_min    = 0.0       ; t_max    = 10000.0;
h_min    = 0e3       ; h_max    = 1000e3;

lon_min  =   0.0*d2r ; lon_max  = 360.0*d2r;
glat_min = -90.0*d2r ; glat_max =  90.0*d2r;

v_min    = 0.0       ; v_max    = 15000.0;
gam_min  = 0.0*d2r ; gam_max  = 90.0*d2r;
al_min   = -2000*d2r ; al_max   = 2000*d2r;
m_min    = spaceplane_dry_mass_min+auxdata.us_mass     ; m_max    = spaceplane_dry_mass_max+auxdata.booster_dry_mass+auxdata.us_mass+spaceplane_fuel_mass_max+booster_fuel_mass_max;

aoa_min  = -5.0*d2r  ; aoa_max  = 20.0*d2r;
bank_min = 0.0*d2r   ; bank_max = 0.0*d2r;
daoa_min  = -0.5*d2r ; daoa_max = 0.5*d2r;
dbank_min = 0.0*d2r  ; dbank_max= 0.0*d2r;



%------- Phase 1 -------%
ph=1;N=10;

t_0    = 0.0           ; t_f    = 5.0;

h_0    = 0.01e3           ; h_f    = 9e3;
v_0    = 0.1         ; v_f    = 210.0;
gam_0  = 89.999*d2r       ; gam_f  = -4.0*d2r;


% h_0    = 1e3           ; h_f    = 9e3;
% v_0    = 10.0         ; v_f    = 210.0;
% gam_0  = 89.0*d2r       ; gam_f  = -4.0*d2r;



m_0    = m_max         ; m_f    = m_max;




% lon_0  = (360.0-74.0059728)*d2r  ; lon_f  = (360.0-122.1697189)*d2r;
% glat_0 = 40.7127753*d2r     ; glat_f = 37.4274745*d2r;


lon_0  = (360.0-118.2436849)*d2r  ; lon_f  = (360.0-122.1697189)*d2r;
glat_0 = 34.0522342*d2r     ; glat_f = 37.4274745*d2r;


al_0   = 160.0*d2r     ; al_f   = 160.0*d2r;

aoa_0  = 0.0*d2r       ; aoa_f  = 0.0*d2r;
bank_0 = 0.0*d2r       ; bank_f = 0.0*d2r;

daoa_0  = 0.0*d2r      ; daoa_f  = 0.0*d2r;
dbank_0 = 0.0*d2r      ; dbank_f = 0.0*d2r;


autofillPhase();

bounds.phase(ph).initialtime.lower = t_0;                                                   % initial conditions
bounds.phase(ph).initialtime.upper = t_0;                                                   % (heading is free)
%bounds.phase(ph).initialstate.lower = [r_0,lon_0,lat_0,v_0,gam_0,al_min,spaceplane_dry_mass_min+auxdata.booster_dry_mass+auxdata.us_mass+spaceplane_fuel_mass_min+booster_fuel_mass_min,aoa_0,bank_0];  %
bounds.phase(ph).initialstate.lower = [r_0,lon_0,lat_0,v_0,gam_0,al_min,170000,aoa_0,bank_0];  %
bounds.phase(ph).initialstate.upper = [r_0,lon_0,lat_0,v_0,gam_0,al_max,spaceplane_dry_mass_max+auxdata.booster_dry_mass+auxdata.us_mass+spaceplane_fuel_mass_max+booster_fuel_mass_max,aoa_0,bank_0];  %

bounds.phase(ph).control.lower(1) = 0;                                                      % force daoa=0
bounds.phase(ph).control.upper(1) = 0;                                                      %

%------- Phase 2 -------%
ph=2;N=50;

t_0    = 5.0           ; t_f    = 210.0;

h_0    = 9e3           ; h_f    = 95e3;
v_0    = 210.0         ; v_f    = 2000.0;
gam_0  = -4.0*d2r      ; gam_f  = 40.0*d2r;
m_0    = m_max         ; m_f    = m_min;

lon_0  = (360.0-122.1697189)*d2r  ; lon_f  = (360.0-122.1697189)*d2r;
glat_0 = 37.4274745*d2r     ; glat_f = 37.4274745*d2r;
al_0   = 160.0*d2r     ; al_f   = 160.0*d2r;

aoa_0  = 10.0*d2r      ; aoa_f  = 15.0*d2r;
bank_0 = 0.0*d2r       ; bank_f = 0.0*d2r;

daoa_0  = 0.0*d2r      ; daoa_f  = 0.0*d2r;
dbank_0 = 0.0*d2r      ; dbank_f = 0.0*d2r;

autofillPhase();

% bounds.phase(ph).finalstate.lower(7) = m_min;  % force to reach min mass
% bounds.phase(ph).finalstate.upper(7) = m_min;  %

bounds.phase(ph).control.lower(1) = 0;                                                      % force daoa=0   !!!!!!!!!!!!!! MAYBE NOT IF STILL IN THE ATMOSPHERE !!!!!!!!!
bounds.phase(ph).control.upper(1) = 0;                                                      %





%
% VEHICLE Gliding
%

t_min    = 0.0       ; t_max    = 10000.0;
h_min    = -40E3     ; h_max    = 1000e3;

lon_min  =   0.0*d2r ; lon_max  = 360.0*d2r;
glat_min = -90.0*d2r ; glat_max =  90.0*d2r;

v_min    = 1.0       ; v_max    = 15000.0;
gam_min  = -90.0*d2r ; gam_max  = 90.0*d2r;
al_min   = -2000*d2r ; al_max   = 2000*d2r;
m_min    = spaceplane_dry_mass_min+auxdata.us_mass  ; m_max    = spaceplane_dry_mass_max+auxdata.us_mass;

aoa_min  = -10*d2r   ; aoa_max  = 55.0*d2r;
bank_min = 0*d2r     ; bank_max = 90*d2r;
daoa_min  = -1*d2r   ; daoa_max = 0.1*d2r;
dbank_min = -5*d2r   ; dbank_max= 5*d2r;

%------- Phase 3 -------%
ph=3;N=10;

t_0    = 210.0           ; t_f    = 220.0;

h_0    = 95e3           ; h_f    = 100e3;
v_0    = 1800.0         ; v_f    = 1950.0;
gam_0  = 40.0*d2r      ; gam_f  = 40.0*d2r;
m_0    = m_min         ; m_f    = m_min;

lon_0  = (360.0-122.1697189)*d2r  ; lon_f  = (360.0-122.1697189)*d2r;
glat_0 = 37.4274745*d2r     ; glat_f = 37.4274745*d2r;
al_0   = 160.0*d2r     ; al_f   = 160.0*d2r;

aoa_0  = 15.0*d2r      ; aoa_f  = 15.0*d2r;
bank_0 = 0.0*d2r       ; bank_f = 0.0*d2r;

daoa_0  = 0.0*d2r      ; daoa_f  = 0.0*d2r;
dbank_0 = 0.0*d2r      ; dbank_f = 0.0*d2r;

autofillPhase();

bounds.phase(ph).control.lower(1) = 0;  % force daoa=0
bounds.phase(ph).control.upper(1) = 0;  %


%------- Phase 4 -------%
ph=4;N=50;

t_0    = 220.0           ; t_f    = 467.0;

h_0    = 100e3           ; h_f    = 20e3;
v_0    = 1950.0         ; v_f    = 2200.0;
gam_0  = 40.0*d2r      ; gam_f  = -40.0*d2r;
m_0    = m_min         ; m_f    = m_min;

lon_0  = (360.0-122.1697189)*d2r  ; lon_f  = (360.0-122.1697189)*d2r;
glat_0 = 37.4274745*d2r     ; glat_f = 37.4274745*d2r;
al_0   = 160.0*d2r     ; al_f   = 160.0*d2r;

aoa_0  = 15.0*d2r      ; aoa_f  = 2.0*d2r;
bank_0 = 0.0*d2r       ; bank_f = 0.0*d2r;

daoa_0  = 0.0*d2r      ; daoa_f  = 0.0*d2r;
dbank_0 = 0.0*d2r      ; dbank_f = 0.0*d2r;

autofillPhase();

% bounds.phase(ph).control.lower(1) = 0*d2r;                                           % 0 < daoa < 1
% bounds.phase(ph).control.upper(1) = 1*d2r;                                           %
% bounds.phase(ph).control.lower(2) = 0*d2r;                                           % dbank = 0
% bounds.phase(ph).control.upper(2) = 0*d2r;                                           %


bounds.phase(ph).control.lower(1) = 0*d2r;                                           % 0 < daoa < 1
bounds.phase(ph).control.upper(1) = 4*d2r;                                           %
bounds.phase(ph).control.lower(2) = 0*d2r;                                           % dbank = 0
bounds.phase(ph).control.upper(2) = 0*d2r;                                           %



%------- Phase 5 -------%
ph=5;N=50;

t_0    = 467.0          ; t_f    = 530.0;

h_0    = 20e3           ; h_f    = 20e3;
v_0    = 2200.0         ; v_f    = 250.0;
gam_0  = -4.0*d2r      ; gam_f  = 40.0*d2r;
m_0    = m_min         ; m_f    = m_min;

lon_0  = (360.0-122.1697189)*d2r  ; lon_f  = (360.0-122.1697189)*d2r;
glat_0 = 37.4274745*d2r     ; glat_f = 37.4274745*d2r;
al_0   = 160.0*d2r     ; al_f   = 160.0*d2r;

aoa_0  = 2.0*d2r       ; aoa_f  = 10.0*d2r;
bank_0 = 0.0*d2r       ; bank_f = 0.0*d2r;

daoa_0  = 0.0*d2r      ; daoa_f  = 0.0*d2r;
dbank_0 = 0.0*d2r      ; dbank_f = 0.0*d2r;

autofillPhase();




% bounds.phase(ph).control.lower(1) = -4*d2r;                                          % -4 < daoa < 0
% bounds.phase(ph).control.upper(1) = 0*d2r;                                           %
% bounds.phase(ph).control.lower(2) = -5*d2r;                                          % -5 < dbank < 5
% bounds.phase(ph).control.upper(2) = 5*d2r;                                           %

% bounds.phase(ph).state.lower(9) = 0*d2r;                                             % 0 < bank < 90
% bounds.phase(ph).state.upper(9) = 90*d2r;                                            %

% bounds.phase(ph).state.lower(8) = 10*d2r;                                            % aoa > 10deg for M > 1

% bounds.phase(ph).finalstate.lower(4) = 250;                                          % velocity
% bounds.phase(ph).finalstate.upper(4) = 300;                                          %


bounds.phase(ph).control.lower(1) = -4*d2r;                                          % -4 < daoa < 0
bounds.phase(ph).control.upper(1) = 0*d2r;                                           %
bounds.phase(ph).control.lower(2) = 0*d2r;                                           % 0 < dbank < 0
bounds.phase(ph).control.upper(2) = 0*d2r;                                           %


%------- Phase 6 -------%
ph=6;N=50;

t_0    = 530.0           ; t_f    = 1170.0;

h_0    = 20e3           ; h_f    = 0e3;
v_0    = 250.0         ; v_f    = 1.0;
gam_0  = 40.0*d2r      ; gam_f  = -10.0*d2r;
m_0    = m_min         ; m_f    = m_min;

lon_0  = (360.0-122.1697189)*d2r  ; lon_f  = (360.0-122.1697189)*d2r;
glat_0 = 37.4274745*d2r     ; glat_f = 37.4274745*d2r;
al_0   = 160.0*d2r     ; al_f   = 160.0*d2r;

aoa_0  = 10.0*d2r      ; aoa_f  = 2.0*d2r;
bank_0 = 0.0*d2r       ; bank_f = 0.0*d2r;

daoa_0  = 0.0*d2r      ; daoa_f  = 0.0*d2r;
dbank_0 = 0.0*d2r      ; dbank_f = 0.0*d2r;

autofillPhase();


% bounds.phase(ph).control.lower(1) = -0.1*d2r;                                        % -0.1 < daoa < 0.1
% bounds.phase(ph).control.upper(1) = 0.1*d2r;                                         %
% bounds.phase(ph).control.lower(2) = -1*d2r;                                          % -1 < dbank < 1
% bounds.phase(ph).control.upper(2) = 1*d2r;                                           %

% bounds.phase(ph).state.lower(8) = 0*d2r;                                             % 0 < aoa < 15
% bounds.phase(ph).state.upper(8) = 15*d2r;                                            %
% bounds.phase(ph).state.lower(9) = -45*d2r;                                           % 0 < bank < 45
% bounds.phase(ph).state.upper(9) = 45*d2r;                                            %

% bounds.phase(ph).state.lower(5) = -20*d2r;                                           % -20 < gam



bounds.phase(ph).control.lower(1) = -1*d2r;                                          % -1 < daoa < 1
bounds.phase(ph).control.upper(1) = 1*d2r;                                           %
bounds.phase(ph).control.lower(2) = -1*d2r;                                          % -1 < dbank < 1
bounds.phase(ph).control.upper(2) = 1*d2r;                                           %

bounds.phase(ph).state.lower(8) = 0*d2r;                                             % 0 < aoa < 20
bounds.phase(ph).state.upper(8) = 20*d2r;                                            %
bounds.phase(ph).state.lower(9) = -45*d2r;                                           % -45 < bank < 45
bounds.phase(ph).state.upper(9) = 45*d2r;                                            %

bounds.phase(ph).state.lower(5) = -20*d2r;                                           % -20 < gam






%
% Guess and Bounds
%


guess.parameter = [0.0,0.0,0.0,0.0,25000.0,400000.0];



h = 300; %699;
lon = 151.171605*d2r; % (360.0-117.85648998651124)*d2r;
glat = -33.929492*d2r; % 34.919316233549665*d2r;




[r_lnd,lon_lnd,lat_lnd] = geolat(int32(1),auxdata.re,auxdata.f_ell,h,lon,glat);

bounds.phase(ph).finalstate.lower(1) = r_lnd+0;                                   % altitude (landing conditions)
bounds.phase(ph).finalstate.upper(1) = r_lnd+200;                                   %

bounds.phase(ph).finalstate.lower(2) = lon_lnd;                                     % longitude
bounds.phase(ph).finalstate.upper(2) = lon_lnd;                                     %
bounds.phase(ph).finalstate.lower(3) = lat_lnd;                                     % latitude
bounds.phase(ph).finalstate.upper(3) = lat_lnd;                                     %

%force_loop_number = 2;
%bounds.phase(ph).finalstate.lower(6) = (238.0+force_loop_number*360.0)*d2r;                          % heading     238
%bounds.phase(ph).finalstate.upper(6) = (238.0+force_loop_number*360.0)*d2r;                          %

bounds.phase(ph).finalstate.lower(5) = -20*d2r;                                     % fpa
bounds.phase(ph).finalstate.upper(5) = -10*d2r;                                     %
bounds.phase(ph).finalstate.lower(9) = 0*d2r;                                       % bank
bounds.phase(ph).finalstate.upper(9) = 0*d2r;                                       %

%
% integral bounds
%

for k=4:6
	guess.phase(k).integral = 0;
	bounds.phase(k).integral.lower = 0;
	bounds.phase(k).integral.upper = 10000;
end

%
% constraints on pdyn,hr,nz
%

% trim_fwd,trim_aft,sm_fwd,sm_aft

pdyn_min  = 0.0      ; pdyn_max = 100E3;
hr_min    = 0.0      ; hr_max   = 100E3;
nz_min  = 0.0        ; nz_max = 10;

pdyn_min  = 0.0      ; pdyn_max = 100E3;
hr_min    = 0.0      ; hr_max   = 100E3;
nz_min  = 0.0        ; nz_max = 50;

for k=1:3

       bounds.phase(k).path.lower = [pdyn_min,hr_min,nz_min];
       bounds.phase(k).path.upper = [pdyn_max,hr_max,nz_max];

       % bounds.phase(k).path.lower = [pdyn_min,hr_min,nz_min,-4.0,-4.0];
       % bounds.phase(k).path.upper = [pdyn_max,hr_max,nz_max, 6.0, 6.0];

	% bounds.phase(k).path.lower = [pdyn_min,hr_min,nz_min,-4.0,-4.0,-10.0,-10.0];
	% bounds.phase(k).path.upper = [pdyn_max,hr_max,nz_max, 6.0, 6.0,  4.0,  4.0];

end

pdyn_min  = 0.0      ; pdyn_max = 10E3;       % 10E3;
hr_min    = 0.0      ; hr_max   = 80E3;       % 65E3;
nz_min    = 0.0      ; nz_max   = 5;          % 5.0;


pdyn_min  = 0.0      ; pdyn_max = 100E3;       % 10E3;
hr_min    = 0.0      ; hr_max   = 10E3;       % 65E3;
nz_min    = 0.0      ; nz_max   = 50;          % 5.0;

for k=4:6

       bounds.phase(k).path.lower = [pdyn_min,hr_min,nz_min];
       bounds.phase(k).path.upper = [pdyn_max,hr_max,nz_max];

       % bounds.phase(k).path.lower = [pdyn_min,hr_min,nz_min,-4.0,-4.0];
       % bounds.phase(k).path.upper = [pdyn_max,hr_max,nz_max, 6.0, 6.0];

	% bounds.phase(k).path.lower = [pdyn_min,hr_min,nz_min,-4.0,-4.0,-10.0,-10.0];
	% bounds.phase(k).path.upper = [pdyn_max,hr_max,nz_max, 6.0, 6.0,  4.0,  4.0];

end

%%
%% bounds for event constraints
%%
for k=1:5                                                % phase links
	bounds.eventgroup(k).lower = zeros(1,10);            %
	bounds.eventgroup(k).upper = zeros(1,10);            %
end                                                      %

% bounds.eventgroup(3).lower = [zeros(1,7),0,0,0];         % mass discontinuity
% bounds.eventgroup(3).upper = [zeros(1,7),0,0,0];         % at SOAR/UST separation

bounds.eventgroup(6).lower = [100, 100,  1200,10  ,10  ,10];   % phase durations
bounds.eventgroup(6).upper = [3000,3000,10000,5000,5000,5000]; % (fixed time between MECO and US release)

bounds.eventgroup(7).lower = [-20.0*d2r,  3000.0,  100000.0];         % gam, speed, altitude
bounds.eventgroup(7).upper = [ 50.0*d2r, 90000.0, 1000000.0];         %

bounds.eventgroup(8).lower = [0];             % constraint on |daoa| & |dbnk|
bounds.eventgroup(8).upper = [6];             % DIFFERENT w.r. to dsc

bounds.eventgroup(9).lower = [0];             % dry mass
bounds.eventgroup(9).upper = [0];             %

bounds.eventgroup(10).lower = [0];             % fuel limitations
bounds.eventgroup(10).upper = [0];             %

bounds.eventgroup(11).lower = [0];             % fuel limitations
bounds.eventgroup(11).upper = [0];             %

bounds.eventgroup(12).lower = [165.0*d2r]; %[238.0*d2r];             % landing heading
bounds.eventgroup(12).upper = [165.0*d2r]; %[238.0*d2r];             %







%%
%% optimizer setup
%%
setup.name = 'trajectory-optimization';
setup.functions.continuous = @continuous;
setup.functions.endpoint = @endpoint;
setup.nlp.solver = 'snopt';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-4;
setup.mesh.colpointsmin = 4;
setup.mesh.colpointsmax = 50;
setup.derivatives.supplier = 'sparseFD';
setup.derivatives.derivativelevel = 'first';
setup.derivatives.dependencies = 'sparse';
setup.scales.method = 'automatic-bounds';

setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
setup.mesh.phase = mesh;


%%
%% solution
%%

%
% using precomputed guesses and meshing
%

restart = 2;

if restart == 1

	init = load('../../ORIGINAL/20150730_s3to_asc_g/output.mat');
	for k=1:3
		setup.guess.phase(k).time=init.output.result.solution.phase(k).time;
		setup.guess.phase(k).state=init.output.result.solution.phase(k).state;
		setup.guess.phase(k).control=init.output.result.solution.phase(k).control;
		setup.mesh.phase(k)=init.output.result.setup.mesh.phase(k);
	end
	tstg = init.output.result.solution.phase(3).time(end);

	init = load('../../ORIGINAL/20150730_s3to_dsc_g/output.mat');
	for k=4:6
		setup.guess.phase(k).time=init.output.result.solution.phase(k-3).time+tstg;
		setup.guess.phase(k).state=init.output.result.solution.phase(k-3).state;
		setup.guess.phase(k).control=init.output.result.solution.phase(k-3).control;
		setup.mesh.phase(k)=init.output.result.setup.mesh.phase(k-3);
	end

elseif restart == 2

	init = load('output_solution.mat');

       setup.guess = init.output.result.solution;
       setup.mesh = init.output.result.setup.mesh;

elseif restart == 3

	init = load('output_solution.mat');

	% for k=1:6
	% 	setup.guess.phase(k).time = init.output.result.solution.phase(k).time;
	% 	state = init.output.result.solution.phase(k).state;
	% 	setup.guess.phase(k).state = state(:,1:end-4);
	% 	control = init.output.result.solution.phase(k).control;
	% 	setup.guess.phase(k).control = control(:,1:end-4);
	% 	setup.mesh.phase(k) = init.output.result.setup.mesh.phase(k);
	% end

       init.output.result.solution.parameter = [0.5,-0.27,-0.5,0.1,0.0,0.0,400000.0,25000.0]; % SHOULD PUT THE CONVERGED GUESS HERE

       setup.guess = init.output.result.solution;
       setup.mesh = init.output.result.setup.mesh;

end

setup.nlp.snoptoptions.maxiterations = 20;
setup.mesh.maxiterations = 1;


compute_first_run = 0;

if compute_first_run == 1
	output = gpops2(setup);

       %%
       %% save & plot
       %%

       disp('Output')
       write_outputs(output.result.interpsolution,[1:6],'soar',auxdata);
       disp('klm')
       make_kml(output.result);
       disp('Solution')
       save('output','output');
       disp('Done')

else
	output = init.output;
end

% clear mex; % clear mex memory

%%
%% optional smoothing and filtering
%%

smooth = 0;
filter = 0;

if smooth == 1

	sp = [0.001 0.01 0.001 0.00001 0.01 0.000001];
       N  = [10 50 10 50 50 50];

	for k=1:6

              disp(['Smoothing Phase ',num2str(k)])

%		tt = output.result.solution.phase(k).time;
              tt = output.result.nextsetup.guess.phase(k).time;

		if filter == 1
                     tt_filtered = linspace(tt(1),tt(end),N(k)+1)';
		else
			tt_filtered = tt;
		end

		%
		% Solution
		%

		for s = 1:9

%			xx = output.result.solution.phase(k).state(:,s);
                     xx = output.result.nextsetup.guess.phase(k).state(:,s);
			fxx = fit(tt,xx,'smoothingspline',fitoptions('Method','Smooth','SmoothingParam',sp(k)));

			xx_filtered = fxx(tt_filtered);
			new_phase(k).state(:,s) = xx_filtered;

			if s > 7
				dxx_filtered = differentiate(fxx,tt_filtered);
				new_phase(k).control(:,s-7) = dxx_filtered;
			end

		end

%		output.result.solution.phase(k).time = tt_filtered;
%		output.result.solution.phase(k).state = new_phase(k).state;
%		output.result.solution.phase(k).control = new_phase(k).control;

              output.result.nextsetup.guess.phase(k).time = tt_filtered;
              output.result.nextsetup.guess.phase(k).state = new_phase(k).state;
              output.result.nextsetup.guess.phase(k).control = new_phase(k).control;


		%
		% Mesh
		%

		if filter == 1
			mesh_phase(k).colpoints = N(k);
			mesh_phase(k).fraction = 1;
		end

	end

	if filter == 1
%		output.result.setup.mesh.phase = mesh_phase;
              output.result.nextsetup.mesh.phase = mesh_phase;
	end

       disp('Solution Smoothed')
       save('output','output');
       disp('Done')

end


%%
%% 2nd run
%%


% initialitaion (linterp)
% get_cd_cl_rho_p();
get_rnose();

compute_second_run = 1;

if compute_second_run == 1

       % setup.guess = output.result.solution;     % using previous solution and meshing
       % setup.mesh = output.result.setup.mesh;    %

       for iteration_index = 1:5

              setup.guess = output.result.nextsetup.guess;     % using previous solution and meshing
              setup.mesh = output.result.nextsetup.mesh;    %

              setup.nlp.snoptoptions.maxiterations = 20;
              setup.mesh.maxiterations = 0;

              output = gpops2(setup);

              %%
              %% save & plot
              %%

              disp('Output')
              %write_outputs(output.result.interpsolution,[1:6],'soar',auxdata);
              write_outputs(output.result.solution,[1:6],'soar',auxdata);
              disp('klm')
              make_kml(output.result);
              disp('Solution')
              save(['output_',num2str(iteration_index)],'output');
              disp('Done')

%              setup.guess = output.result.nextsetup.guess;     % using previous solution and meshing
%              setup.mesh = output.result.nextsetup.mesh;    %

       end

end







% clear mex; % clear mex memory










