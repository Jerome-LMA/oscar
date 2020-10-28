function varargout = Fit_TEM00(Ein,varargin)
%  Fit_TEM00() find the beam radius and the wavefront curvature for a TEM00
%  beam non astigmatic
%
% Required the 'Optimisation Toolbox' to work

if (~exist('lsqcurvefit','file'))
    error('Fit_TEM00() required the Optimisation Toolbox to run')
end

p = inputParser;
p.FunctionName = 'Fit the TEM00 beam parameters';

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));
p.parse(Ein,varargin{:})

if Calculate_Power(Ein)==0
    error('Fit_TEM00(): No power in the field')
end

% Find if the mode is a TEM00 or not
[~,m,n] = Read_mode_name(Ein.Mode_name);

if (m~=0) || (n~=0)
    error('Fit_TEM00(): Can only be used to fit the fundamental mode. For higher order mode use Fit_E_Field()')
end

% For better accuracy, normalize the beam
E =  Normalise_E(Ein);

Radius_guess = E.Grid.Length/2;

% Take only the central part of the grid
tmp_index = find (E.Grid.D2_square < Radius_guess^2);
tmp_grid = E.Grid.D2_square(tmp_index);
Power_distri = abs(E.Field(tmp_index)).^2;

% Fit the power first
func_gauss = @(c,xdata) c(1)*(exp(-2*xdata/(c(2)^2)) );

options = optimset('Display','off','MaxFunEvals',1E6,'TolFun',1E-8,'DiffMinChange',1E-8);
c0 = [max(max(Power_distri)) Radius_guess];
[Map.fit_para,~,~,~,~] = lsqcurvefit(func_gauss,c0,tmp_grid,Power_distri,[],[],options);

Beam_rad = abs(Map.fit_para(2));

% Now find the RofC of the complex wavefront
clear tmp_index
tmp_index = find (E.Grid.D2_square < (Beam_rad*3)^2);
tmp_grid = E.Grid.D2_square(tmp_index);
tmp_amp = E.Field(tmp_index);

% first find a starting wavefront RofC
% So take a cross section of the phase in one direction
Cross_sec_phase_x = unwrap(angle(E.Field(E.Grid.Half_num_point,:)));

% Will do the fit on the beam diameter, find the index
Cross_sec_index = intersect(find (E.Grid.Axis < Beam_rad),find (E.Grid.Axis > -Beam_rad));

%
% figure(2)
% plot(E.Grid.Axis(Cross_sec_index),Cross_sec_phase_x(Cross_sec_index))

[polyscale,~, mu]  = polyfit(E.Grid.Axis(Cross_sec_index),Cross_sec_phase_x(Cross_sec_index),2);
beam_radius_fit = -E.k_prop/(2*polyscale(1)/mu(2)^2); % first approximation in 1D


% Then fit the phase
func_gauss = @(c,xdata) c(1)*exp(-xdata/Beam_rad^2).*exp(-1i*E.k_prop*xdata/(2*c(2)));

options = optimset('Display','off','MaxFunEvals',1E6,'TolFun',1E-25,'DiffMinChange',1E-25);
c0 = [max(max(E.Field)) beam_radius_fit];
[Map.fit_para,~,~,~,~] = lsqcurvefit(func_gauss,c0,tmp_grid,tmp_amp,[],[],options);
Beam_RofC = real(Map.fit_para(2));

switch nargout
    case 0
        fprintf('Beam radius [m]: %g  \t \t Wavefront curvature [m]: %g  \n',Beam_rad,Beam_RofC)
    case 1
        varargout{1} = Beam_rad;
    case 2
        varargout{1} = Beam_rad;
        varargout{2} = Beam_RofC;
    otherwise
        error('Fit_TEM00(): Too many output argument')
end

end
