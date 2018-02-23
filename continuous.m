
%--------------------------------------%
% BEGIN: function continuous           %
%--------------------------------------%

function phaseout = s3toContinuous(input)

  Nph = length(input.phase);

  for k=1:Nph

    [dr,dlon,dlat,dv,dgam,dal,dm,da,db,pdyn,hr,nx,ny,nz,thu,cd,cl,Fdrag,Flift,rho,p,Tenv,mach,rey1m,trim_fwd,trim_aft,ka_fwd,ka_aft] = dynamics(input.phase(k),input.phase(k).parameter,input.auxdata,k);

    phaseout(k).dynamics = [dr,dlon,dlat,dv,dgam,dal,dm,da,db];


    phaseout(k).path = [pdyn,hr,abs(nz)];
    % phaseout(k).path = [pdyn,hr,abs(nz),trim_fwd,trim_aft];
    % phaseout(k).path = [pdyn,hr,abs(nz),trim_fwd,trim_aft,ka_fwd,ka_aft];

  end

  for k=4:6

    daoa=input.phase(k).control(:,1);
    dbnk=input.phase(k).control(:,2);

    phaseout(k).integrand = abs(daoa)+abs(dbnk);

  end

end

%------------------------------------%
% END: function continuous           %
%------------------------------------%


