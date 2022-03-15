function Eout = Normalise_E(Ein,varargin)
% Normalise a E_field
% Use Normalise_E(E_field)    Normalises the E_field to 1W
%        Normalise_E(E_field,'power',P)    Set the E_field to have a
%        power of P Watt

p  = inputParser;

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check which power we would like to normalise to, by default 1W
p.addParameter('Power',1, @(x)(isreal(x) && ge(x,0)) );



p.parse(Ein,varargin{:})

Eout = Ein;
Amplitude_scaling = sqrt( p.Results.Power / Calculate_Power(Ein) );

if exist('isgpuarray','file') % to be compatible with version < 2020b
    Run_on_GPU = isgpuarray(Ein.Field);
else
    Run_on_GPU = false;
end

if Run_on_GPU
    Eout.Field = arrayfun(@times,Eout.Field,Amplitude_scaling);
    if Amplitude_scaling == 0 % put all the SBs to zero otherwise do not change them
        for ii = 1:Ein.Nb_Pair_SB
            Eout.SB(ii).Field_lower = zeros(Ein.Grid.Num_point,'gpuArray');
            Eout.SB(ii).Field_upper = zeros(Ein.Grid.Num_point,'gpuArray');
        end
    end
else
    Eout.Field = Eout.Field * Amplitude_scaling;
    if Amplitude_scaling == 0 % put all the SBs to zero otherwise do not change them
        
        for ii = 1:Ein.Nb_Pair_SB
            Eout.SB(ii).Field_lower = zeros(Ein.Grid.Num_point);
            Eout.SB(ii).Field_upper = zeros(Ein.Grid.Num_point);
        end
    end
end

end

