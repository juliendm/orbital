
%-----------------------------------------------%
% Begin Function: compute_dry_mass              %
%-----------------------------------------------%

function dry_mass = compute_dry_mass(dvs, dry_mass_current_step, auxdata)

    % Structure Mass

    sizes = size(dvs);
    n_load_cases = sizes(1);
    n_members = 132;
    material_density = 2780.0;

    % Area

    areas = zeros(n_members); % AREA IS ONLY DEPENDENT ON GEOMETRIC DVS !!!!!!!!!!!!!
    dvs_geom = dvs(1,:);
    % USE DEFAULT VALUE FOR ALL DVS OTHER THAN GEOMETRIC
    dvs_geom(1) = 3.0; dvs_geom(2) = 0.0; dvs_geom(3) = 0.0; dvs_geom(4) = 3.0; dvs_geom(5) = -3.0; dvs_geom(6) = 1.5e6; dvs_geom(7) = 20.0e3; dvs_geom(8) = 20.0e3; 
    for member_index = 1:n_members
        if member_index == 29
            areas(member_index) = 3.328734e-01;
        elseif member_index == 87
            areas(member_index) = 3.190315e-01;
        elseif member_index == 88
            areas(member_index) = 3.467154e-01;
        elseif member_index == 119
            areas(member_index) = 1.047058e+00;
        elseif member_index == 122
            areas(member_index) = 1.386107e+00;
        else
            areas(member_index) = surfpack_eval(['area_',num2str(member_index,'%03d')],dvs_geom);
        end
    end

    % areas = zeros(n_members,n_load_cases);
    % for load_case_index = 1:n_load_cases
    %     dvs_load_case_i = dvs(load_case_index,:);
    %     for member_index = 1:n_members
    %         areas(member_index,load_case_index) = surfpack_eval(['area_',num2str(member_index)],dvs_load_case_i);
    %     end
    % end

    % Thickness

    thicknesses = zeros(n_members,n_load_cases);
    for load_case_index = 1:n_load_cases
        dvs_load_case_i = dvs(load_case_index,:);
        for member_index = 1:n_members
            if any([29,87,88,119,122] == member_index)
                thicknesses(member_index,load_case_index) = 1.600000e-03;
            else
                thicknesses(member_index,load_case_index) = surfpack_eval(['thickness_',num2str(member_index,'%03d')],dvs_load_case_i);
            end
        end
    end

    % Structure

    half_structure_mass = 0.0;
    for member_index = 1:n_members
        area = areas(member_index);
        thickness = max(thicknesses(member_index,:));
        half_structure_mass = half_structure_mass + area*thickness*material_density;
    end

    % Additional

    half_additional_mass = compute_half_additional_mass(dvs, dry_mass_current_step, auxdata);

    % Dry Mass

    % BE CONSISTENT: DRY MASS: WITHOUT PAYLOAD MASS !!!!!!!!!!!!!!!!!!!!!!!

    dry_mass = 2.0*half_structure_mass + 2.0*half_additional_mass;    

    dry_mass = dry_mass*0.5; % !!!!!!!!!!!!!!!! OTHERWISE TOO HEAVY AND CHANGE ALL MISSION

end

%-----------------------------------------------%
% End Function:  dry_mass                       %
%-----------------------------------------------%


function mass = compute_half_additional_mass(dvs, dry_mass_current_step, auxdata)

    % % GO THROUGH IT AGAIN TO MAKE SURE CLEAR DIFFERENCE BETWEEN CURRENT MASS AND MASS WHEN ALL FUEL IS IN

    thrust_index = 6;
    pdyn_index = 7;
    fuel_mass_index = 8;

    dv1_index = 9;
    dv2_index = 10;
    dv3_index = 11;
    dv4_index = 12;
    dv5_index = 13;
    dv6_index = 14;

    half_dry_mass_current_step = 0.5*dry_mass_current_step;
    half_mass_payload = 0.5 * auxdata.us_mass;

    max_fuel_mass = auxdata.spaceplane_fuel_mass_max;
    %fuel_percentage = dvs(uel_mass_index) / max_fuel_mass;
    half_mass_fuel_kero = 0.5 * max_fuel_mass * 0.4;
    half_mass_fuel_lox = 0.5 * max_fuel_mass * 0.6;

    half_max_thrust_newtons =  0.5 * max(dvs(:,thrust_index));
    traj_max_pdyn_inf = max(dvs(:,pdyn_index));

    pounds_to_kg = 0.453592;
    newtons_to_pounds = 0.224809;
    kg_to_pounds = 2.20462;
    meters_to_feet = 3.28084;

    wet_weight_pounds_current_step = 2.0*(half_dry_mass_current_step+half_mass_payload+half_mass_fuel_kero+half_mass_fuel_lox)*kg_to_pounds;

    max_thrust_pounds = 2.0*half_max_thrust_newtons*newtons_to_pounds;
    body_length_feet = 18.0*meters_to_feet;
    lox_density = 1141; % kg/m3
    kero_density = 810; % kg/m3
    N_engines = 1;
    rocket_expansion_ratio = 77.5;
    sts_tank_radius = 4.2; % m
    sts_tank_height = 46.9; % m
    sts_tank_area = 2.0*pi*sts_tank_radius*sts_tank_height + 2.0*pi*sts_tank_radius^2.0; % m2 # Assume: Right Cylinder
    sts_tank_empty_mass = 26500.0; % kg
    tank_mass_per_area = sts_tank_empty_mass/sts_tank_area; % kg/m2
    tank_mass_per_area = tank_mass_per_area*0.6; % COMPOSITE REDUCE MASS BY 40 %
    technology_improvement = 0.2; % KEEP ONLY 20 % OF MASS

    % Landing Gear Weight

    weight_gear = 0.00916*wet_weight_pounds_current_step^1.124;
    half_mass_gear = weight_gear*0.5*pounds_to_kg;

    % Avionics Weight

    weight_avionics = 66.37*wet_weight_pounds_current_step^0.361;
    half_mass_avionics = weight_avionics*0.5*pounds_to_kg * technology_improvement;

    % Electrical System Weight

    weight_elec = 1.167*wet_weight_pounds_current_step^0.5*body_length_feet^0.25;
    half_mass_elec = weight_elec*0.5*pounds_to_kg * technology_improvement;

    % Equipment Weight

    weight_equip = 1000.0 + 0.01*wet_weight_pounds_current_step;  % Maybe reduce fixed value !!!!!!!!
    half_mass_equip = weight_equip*0.5*pounds_to_kg * technology_improvement;

    % Tank LOX Weight

    volume_lox = 2.0*half_mass_fuel_lox/lox_density; % m3
    lox_tank_radius = 1.65; % m
    lox_tank_height = volume_lox/pi/lox_tank_radius/lox_tank_radius;
    lox_tank_area = 2.0*pi*lox_tank_radius*lox_tank_height + 2.0*pi*lox_tank_radius^2.0; % m2 % Assume: Right Cylinder
    half_mass_tank_lox = 0.5*tank_mass_per_area*lox_tank_area; % kg

    % Tank KERO Weight

    volume_kero = 2.0*half_mass_fuel_kero/kero_density; % m3
    kero_tank_radius = 1.65; % m
    kero_tank_height = volume_kero/pi/kero_tank_radius/kero_tank_radius;
    kero_tank_area = 2.0*pi*kero_tank_radius*kero_tank_height + 2.0*pi*kero_tank_radius^2.0; % m2 % Assume: Right Cylinder
    half_mass_tank_kero = 0.5*tank_mass_per_area*kero_tank_area; % kg

    % Engine Weight

    weight_engine = 0.00766*max_thrust_pounds + 0.00033*max_thrust_pounds*rocket_expansion_ratio^0.5 + 130.0*N_engines;
    half_mass_engine = weight_engine*0.5*pounds_to_kg;


    % Surface Dependant

    % TODO: RESPONSE SURFACE !!!!!!!!!!!!!!!!!!!!!  IS IT HALF ???????
    half_mass_tps = 100.0; % surfpack_eval(['tps',dvs);               
    half_mass_hydraulic = 100.0; %surfpack_eval(['hydraulic',dvs);







    mass = 0.0;
    mass = mass + half_mass_gear;
    mass = mass + half_mass_avionics;
    mass = mass + half_mass_elec;
    mass = mass + half_mass_equip;
    mass = mass + half_mass_tank_lox;
    mass = mass + half_mass_tank_kero;
    mass = mass + half_mass_engine;
    mass = mass + half_mass_hydraulic;
    mass = mass + half_mass_tps;

end








