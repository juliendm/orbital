%-----------------------------------------------%
% Begin Function:  geo2lat                      %
%-----------------------------------------------%
function [r_0,r_f,r_min,r_max,lat_0,lat_f,lat_min,lat_max] =... 
    geo2lat (h_0,h_f,h_min,h_max,lon_0,lon_f,lon_min,lon_max,...
             glat_0,glat_f,glat_min,glat_max,re,f_ell)




h0 = [h_0,h_f,h_min,h_max];
lon0 = [lon_0,lon_f,lon_min,lon_max];
glat0 = [glat_0,glat_f,glat_min,glat_max];

[r,lon,lat] = geolat(int32(4),double(re),double(f_ell),...
    double(h0),double(lon0),double(glat0));

r_0=r(1)    ; r_f=r(2)    ; r_min=r(3)    ; r_max=r(4);
lat_0=lat(1); lat_f=lat(2); lat_min=lat(3); lat_max=lat(4);

end
%-----------------------------------------------%
% End Function:  geo2lat                        %
%-----------------------------------------------%
