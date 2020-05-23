function [phasemap, dx, buckets] = ReadZygoBinary(filename)
% [phasemap, dx, buckets] = ReadZygoBinary(filename) read ZYGO binay maps
%
% calls `LoadZygoBinary' and then rescales phase data accordingly. The
% output are:
% `phasemap' is the rescaled phase map
% `dx' is the lateral resolution
% `buckets' is a (width x height x n) matrix containing the n buckets
% (intensity data) used to compute the phase map; they are saved in the
% Zygo .dat file only if diagnostic options have been turned on in MetroPro
%
% Author: Massimo Galimberti 2010 (all the merit goes to him)
% Modified by: Jerome Degallaix 2013

data = LoadZygoBinary(filename);

% scale PhaseData
switch data.phase_res
    case 0
        res = 4096;
    case 1
        res = 32768;
    case 2
        res = 131072;
    otherwise
        res = 1;
end
        
% phase data in waves
phasemap = data.PhaseData * data.intf_scale_factor * data.obliquity_factor / res;
% phase data in meters
phasemap = phasemap * data.wavelength_in;

phasemap = rot90(phasemap,1);
phasemap = flipud(phasemap);

dx = data.lateral_res;

% intensity data (if present)
if ~isempty(data.IntensityData)
    buckets = zeros( size(data.IntensityData,2) , size(data.IntensityData,1), size(data.IntensityData,3) );
    for k=1:size(data.IntensityData,3)
        buckets(:,:,k) = rot90( data.IntensityData(:,:,k), 1);
        buckets(:,:,k) = flipud( buckets(:,:,k) );
    end
end
