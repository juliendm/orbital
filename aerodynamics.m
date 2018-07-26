
%-----------------------------------------------%
% Begin Function: aerodynamics                  %
%-----------------------------------------------%

function [cd,cl,rho,p,Tenv,mach,rey1m,el_def,bf_def,trim_fwd,trim_aft,ka_fwd,ka_aft] = aerodynamics(k,n,h,lon,glat,aoa_deg,v,tt,re,date0_doy,date0_sec,atm_model,ar_flag,dv1,dv2,dv3,dv4,dv5,dv6);

    max_mach_eval_sub = 0.95;
    min_mach_eval_sup = 1.1;

    bf_bound_lower = -0.25;
    bf_bound_upper = 0.5;

    [rho,p,ss] = get_rho_p_ss(int32(n),double(h),double(lon),double(glat),...
      double(tt),double(re),int32(date0_doy),double(date0_sec),int32(atm_model));

    rho = reshape(rho,n,1);
    p = reshape(p,n,1);
    ss = reshape(ss,n,1);

    mach = v./ss;

    Tenv = p./rho/287.15;
    Cs = 120.0; T0 = 291.15; mu0 = 18.27e-6;   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHANGE COEFFS
    visc = mu0*(T0+Cs)./(Tenv+Cs).*((Tenv/T0).^(3/2));
    rey1m = (rho.*v)./visc;

    cd = zeros(size(mach));
    cl = zeros(size(mach));
    el_def = zeros(size(mach));
    bf_def = zeros(size(mach));
    trim_fwd = zeros(size(mach));
    trim_aft = zeros(size(mach));
    ka_fwd = zeros(size(mach));
    ka_aft = zeros(size(mach));

    for i = 1:length(mach)

        % Change of Variables

        dv_mach = mach(i);
        dv_rey = 0.0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   log10(rey1m(i)) + 3.0/8.0*dv_mach - 7.0;  

        if (dv_mach >= 8.0)
            dv_mach = 8.0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATABASE NOT VALID BEYOND THAT: SPACE SHUTTLE BOOK FIX BEYOND MACH 10 ANYWAY ...
        end

        if (dv_mach >= 1.0)
            dv_aoa = (aoa_deg(i) - (3.92857*dv_mach+3.57143)) / (1.78571*dv_mach+5.71429); % Supersonic
        else
            dv_aoa = aoa_deg(i)/7.5-1.0; % Subsonic
        end

        dv_geo1 = dv1(i);
        dv_geo2 = dv2(i);
        dv_geo3 = dv3(i);
        dv_geo4 = dv4(i);
        dv_geo5 = dv5(i);
        dv_geo6 = dv6(i);


        %%%%%%%%%
        % Perfo %
        %%%%%%%%%

        if k >= 5 && dv_mach < 5.0

            dvs_perfo = [dv_mach,dv_rey,dv_aoa,dv_geo1,dv_geo2,dv_geo3,dv_geo4,dv_geo5,dv_geo6];

            val = surfpack_eval('trim_fwd',dvs_perfo);
            if val < bf_bound_lower
                el_def(i) = val - bf_bound_lower;
                bf_def(i) = bf_bound_lower;
            elseif val > bf_bound_upper
                el_def(i) = val - bf_bound_upper;
                bf_def(i) = bf_bound_upper;
            else
                el_def(i) = 0.0;
                bf_def(i) = val;
            end
            trim_fwd(i) = el_def(i)*180.0/pi;

            val = surfpack_eval('trim_aft',dvs_perfo);
            if val < bf_bound_lower
                el_def(i) = val - bf_bound_lower;
                bf_def(i) = bf_bound_lower;
            elseif val > bf_bound_upper
                el_def(i) = val - bf_bound_upper;
                bf_def(i) = bf_bound_upper;
            else
                el_def(i) = 0.0;
                bf_def(i) = val;
            end
            trim_aft(i) = el_def(i)*180.0/pi;

            ka_fwd(i) = surfpack_eval('k_alpha_fwd',dvs_perfo);
            ka_aft(i) = surfpack_eval('k_alpha_aft',dvs_perfo);

        else

            el_def(i) = 0.0;
            bf_def(i) = 0.0;

            trim_fwd(i) = 0.0;
            trim_aft(i) = 0.0;

            ka_fwd(i) = 0.0;
            ka_aft(i) = 0.0;

        end

        %%%%%%%%
        % Aero %
        %%%%%%%%

        dvs = [dv_mach,dv_rey,dv_aoa,el_def(i),bf_def(i),dv_geo1,dv_geo2,dv_geo3,dv_geo4,dv_geo5,dv_geo6];

        if (dv_mach <= max_mach_eval_sub)

            cd(i) = surfpack_eval('drag_sub',dvs);
            cl(i) = surfpack_eval('lift_sub',dvs);

        elseif (dv_mach >= min_mach_eval_sup)

            cd(i) = surfpack_eval('drag_sup',dvs);
            cl(i) = surfpack_eval('lift_sup',dvs);

        else

            dvs_sub = dvs;
            dvs_sub(1) = max_mach_eval_sub;
            dvs_sup = dvs;
            dvs_sup(1) = min_mach_eval_sup;

            cd_sub = surfpack_eval('drag_sub',dvs_sub);
            cd_sup = surfpack_eval('drag_sup',dvs_sup);
            cd(i) = cd_sub + (dv_mach-max_mach_eval_sub) * (cd_sup-cd_sub)/(min_mach_eval_sup-max_mach_eval_sub);

            cl_sub = surfpack_eval('lift_sub',dvs_sub);
            cl_sup = surfpack_eval('lift_sup',dvs_sup);
            cl(i) = cl_sub + (dv_mach-max_mach_eval_sub) * (cl_sup-cl_sub)/(min_mach_eval_sup-max_mach_eval_sub);

        end

        
    end

end

%-----------------------------------------------%
% End Function:  aerodynamics                   %
%-----------------------------------------------%
