function [x,y] = ll2psn(lat,lon,varargin)
% ll2psn (lat/lon to polarstereographic) transforms lat/lon coordinates to north
% polar stereographic coordinates. This function does NOT require Matlab's Mapping
% Toolbox. 
%   
%%% Syntax
% [x,y] = ll2psn(lat,lon) 
% [x,y] = ll2psn(lat,lon,'TrueLat',ReferenceLatitude) 
% [x,y] = ll2psn(lat,lon,'EarthRadius',RadiusInMeters) 
% [x,y] = ll2psn(lat,lon,'Eccentricity',EarthsMisshapenness) 
% [x,y] = ll2psn(lat,lon,'meridian',MeridianInDegrees) 
% 
%%% Description  
% [x,y] = ll2psn(lat,lon) transforms georeferenced coordinates to
% polarstereographic x,y coordinates referenced to 70 N. Inputs lat and lon
% can be scalar, vecotr, or matrices of equal size. 
% 
% [x,y] = ll2psn(lat,lon,'TrueLat',ReferenceLatitude) secifies a reference
% latitude of true scale in degrees; also known as the standard parallel.
% Default is 70 N. 
% 
% [x,y] = ll2psn(lat,lon,'EarthRadius',RadiusInMeters) specifies Earth's
% radius in meters. Default is 6378137.0 m, WGS84.
% 
% [x,y] = ll2psn(lat,lon,'Eccentricity',EarthsMisshapenness) specifies
% Earth's eccentricity or misshappenness.  Default values is 0.08181919. 
% 
% [x,y] = ll2psn(lat,lon,'meridian',MeridianInDegrees) specifies the meridian in 
% degrees along the positive Y axis of the map. Default value is -45.
% 

%%% Input checks: 

assert(nargin>1,'The ll2psn function requires at least two inputs: lat and lon.')
if any(lat(:))<0
   warning('You are transforming at least some datapoints that are in the southern hemisphere, but using a north polar stereographic projection.') 
end
assert(sum(any(lon>360))==0,'Input longitudes greater than 360 degrees? This seems wrong.')
assert(sum(any(lon<-360))==0,'Input longitudes less than -360 degrees? This seems wrong.')

%%% Set defaults: 

phi_c = 70;   % standard parallel - this is different from Andy Bliss' function, which uses -70! 
a = 6378137.0; % radius of ellipsoid, WGS84
e = 0.08181919;% eccentricity, WGS84
lambda_0 = -45;  % meridian along positive Y axis

%%% Parse user inputs: 

tmp = strncmpi(varargin,'true',4)|strncmpi(varargin,'lat',3)|strncmpi(varargin,'ref',3)|...
    strncmpi(varargin,'earthrad',8)|strncmpi(varargin,'rad',3)|...
    strncmpi(varargin,'ecc',3)|strncmpi(varargin,'PosYLon',4)|...
    strncmpi(varargin,'lon',3)|strncmpi(varargin,'merid',5); 
assert(sum(tmp)==(nargin-2)/2,'There seems to be at least one invalid input string. Are you trying to declare options that do not exist?')

tmp = strncmpi(varargin,'true',4)|strncmpi(varargin,'lat',3)|strncmpi(varargin,'ref',3); 
if any(tmp)
    phi_c = varargin{find(tmp)+1}; 
    assert(isscalar(phi_c)==1,'True lat must be a scalar.')
end

tmp = strncmpi(varargin,'earthrad',8)|strncmpi(varargin,'rad',3); 
if any(tmp)
    a = varargin{find(tmp)+1}; 
    assert(isscalar(a)==1,'Earth radius must be a scalar.')
    assert(a > 7e+3,'Earth radius should be something like 6378137 in meters.')
end

tmp = strncmpi(varargin,'ecc',3); 
if any(tmp)
    e = varargin{find(tmp)+1}; 
    assert(isscalar(e)==1,'Earth eccentricity must be a scalar.')
    assert(e>0 & e<1,'Earth eccentricity does not seem like a logical value.')
end

tmp = strncmpi(varargin,'PosYLon',4)|strncmpi(varargin,'lon',3)|...
    strncmpi(varargin,'merid',5); 
if any(tmp)
    lambda_0 = varargin{find(tmp)+1}; 
    assert(isscalar(lambda_0)==1,'PosYLon must be a scalar.')
    assert(lambda_0>=-180 & lambda_0<=360,'PsosYLon does not seem like a logical value.')
end

%%% Let the transformation begin
%%% Convert to radians:    
phi=lat*pi/180;
phi_c=phi_c*pi/180;
lambda=lon*pi/180;
lambda_0=lambda_0*pi/180;

%%%this is not commented very well. See Snyder for details.
t=tan(pi/4-phi/2)./((1-e*sin(phi))./(1+e*sin(phi))).^(e/2);
t_c=tan(pi/4 - phi_c/2)./((1-e*sin(phi_c))./(1+e*sin(phi_c))).^(e/2);
m_c=cos(phi_c)./sqrt(1-e^2*(sin(phi_c)).^2);
rho=a*m_c*t/t_c; %true scale at lat phi_c

x=rho.*sin(lambda-lambda_0);
y=-rho.*cos(lambda - lambda_0);

%%% Make two-column format if user requested fewer than two outputs: 
if nargout<2
    x = [x(:) y(:)]; 
end
