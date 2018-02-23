
%-----------------------------------------------%
% Begin Function: make_kml                      %
%-----------------------------------------------%

function [] = make_kml(sol)

	load('utils/20140902_geoid_height');

	r2d = 180/pi;
	d2r = pi/180;

	int = sol.interpsolution;
	aux = sol.setup.auxdata;

	re = aux.re;
	f_ell = aux.f_ell;


	s = geoshape();
	s.Geometry = 'line';

	Nph = length(int.phase);

	for k=1:Nph

	    traj = int.phase(k).state';

		r = traj(1,:);
		lon = traj(2,:);
		lat = traj(3,:);
		bnk = traj(9,:);
		n = length(r);
		[h,lon,glat] = latgeo(int32(n),double(re),double(f_ell),...
		     double(r),double(lon),double(lat));
		h = reshape(h,n,1);
		lon = reshape(lon,n,1);
		glat = reshape(glat,n,1);

		galt = geoid_height(lat'*r2d,wrapTo2Pi(lon)*r2d);

		s(k).Latitude = glat'*r2d;
		s(k).Longitude = lon'*r2d;
		s(k).Altitude = h'-galt';
		%s(1).Metadata.bnk = bnk*r2d;

	end

	velocity = traj(4,end);

	% colors = jet(1);
	% colors = 'ymcrgbwk';
	% rng('shuffle');
	kmlwrite('output/traj.kml',s,'Name',num2str(velocity),'Color',{'blue','red','blue','blue','blue','blue'},'LineWidth',2);

end

%-----------------------------------------------%
% End Function: make_kml                        %
%-----------------------------------------------%
