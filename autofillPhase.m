%-----------------------------------------------%
% Begin Function:  autofillPhase                %
%-----------------------------------------------%
function [] = autofillPhase ()

	global auxdata ph N bounds guess mesh ...
	       t_0 t_f t_min t_max r_0 r_f r_min r_max h_0 h_f h_min h_max ...
	       lon_0 lon_f lon_min lon_max glat_0 glat_f glat_min glat_max ...
	       lat_0 lat_f lat_min lat_max v_0 v_f v_min v_max gam_0 gam_f gam_min ...
	       gam_max al_0 al_f al_min al_max m_0 m_f m_min m_max aoa_0 aoa_f ...
	       aoa_min aoa_max bank_0 bank_f bank_min bank_max daoa_0 daoa_f ...
	       daoa_min daoa_max dbank_0 dbank_f dbank_min dbank_max ...
           tva_0 tva_f tva_min tva_max dtva_0 dtva_f dtva_min dtva_max ...
	       pdyn_min pdyn_max hr_min hr_max gacc_min gacc_max

	%
	% bounds
	%
	clearvars xxx;

	h0 = [h_0,h_f,h_min,h_max];
	lon0 = [lon_0,lon_f,lon_min,lon_max];
	glat0 = [glat_0,glat_f,glat_min,glat_max];

	[r,lon,lat] = geolat(int32(4),double(auxdata.re),double(auxdata.f_ell),double(h0),double(lon0),double(glat0));
	r_0=r(1);r_f=r(2);r_min=r(3);r_max=r(4);lat_0=lat(1);lat_f=lat(2);lat_min=lat(3);lat_max=lat(4);

	xxx.initialtime.lower = t_min;
	xxx.initialtime.upper = t_max;
	xxx.finaltime.lower = t_min;
	xxx.finaltime.upper = t_max;

	xxx.initialstate.lower = [r_min,lon_min,lat_min,v_min,gam_min,al_min,m_min,aoa_min,bank_min,tva_min]; 
	xxx.initialstate.upper = [r_max,lon_max,lat_max,v_max,gam_max,al_max,m_max,aoa_max,bank_max,tva_max];
	xxx.state.lower =        [r_min,lon_min,lat_min,v_min,gam_min,al_min,m_min,aoa_min,bank_min,tva_min];
	xxx.state.upper =        [r_max,lon_max,lat_max,v_max,gam_max,al_max,m_max,aoa_max,bank_max,tva_max];
	xxx.finalstate.lower =   [r_min,lon_min,lat_min,v_min,gam_min,al_min,m_min,aoa_min,bank_min,tva_min];
	xxx.finalstate.upper =   [r_max,lon_max,lat_max,v_max,gam_max,al_max,m_max,aoa_max,bank_max,tva_max];

	xxx.control.lower = [daoa_min,dbank_min,dtva_min];
	xxx.control.upper = [daoa_max,dbank_max,dtva_max];

	bounds.phase(ph)=xxx;

	%
	% guess
	%
	clearvars xxx;

	xxx.time = [t_0;t_f];

	xxx.state = [...
	[r_0;r_f],...
	[lon_0;lon_f],...
	[lat_0;lat_f],...
	[v_0;v_f],...
	[gam_0;gam_f],...
	[al_0;al_f],...
	[m_0;m_f],...
	[aoa_0;aoa_f],...
	[bank_0;bank_f],...
	[tva_0;tva_f]];

	xxx.control = [[daoa_0;daoa_f],[dbank_0;dbank_f],[dtva_0;dtva_f]];

	guess.phase(ph)=xxx;

	%
	% mesh
	%
	mesh(ph).colpoints= N;
	mesh(ph).fraction = 1;

end
%-----------------------------------------------%
% End Function:  autofillPhase                  %
%-----------------------------------------------%
