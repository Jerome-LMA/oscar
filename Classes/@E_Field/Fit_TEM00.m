function varargout = Fit_TEM00(Ein,varargin)
%  Fit_TEM00() find the beam radius and the wavefront curvature for a TEM00
%  beam non astigmatic
%
% Required the 'Optimisation Toolbox' to work

if (~exist('lsqcurvefit','file'))
    error('Fit_TEM00() required the Optimisation Toolbox to run')
end

p  = inputParser;
p.FunctionName = 'Fit the TEM00 beam parameters';

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check what what we should fit
p.addOptional('for','carrier', @(x)strcmpi(x,'carrier') | strcmpi(x,'SB_upper') | strcmpi(x,'SB_lower') );

% Check the number of the sidebands we want to deal with
p.addParameter('SB_num',1, @(x) isnumeric(x)  && (x>0) && (mod(x,1) == 0)); % check if the number of the SB pair is positive and integer

p.parse(Ein,varargin{:})

if strcmp(p.Results.for,'SB_upper')
    Ein.Field = Ein.SB(p.Results.SB_num).Field_upper;
elseif strcmp(p.Results.for,'SB_lower')
    Ein.Field = Ein.SB(p.Results.SB_num).Field_lower;
end


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

if Cross_sec_index < 2
    Beam_rad
    figure(201)
    plot(E.Grid.Axis,Cross_sec_phase_x)
    error('Fit_TEM00(): beam radius estimation too small to fit the RoC')
end
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

q = 1/ (1/Beam_RofC -1i*1064e-9/pi/Beam_rad^2);

Beam_z2 = real(q);
Beam_zR = imag(q);
Beam_w0 = sqrt(Beam_zR*E.Wavelength/pi);

switch nargout
    case 0
        fprintf('Beam radius [m]: %g  \t \t Wavefront curvature [m]: %g  \n',Beam_rad,Beam_RofC)
    case 1
        varargout{1} = Beam_rad;
    case 2
        varargout{1} = Beam_rad;
        varargout{2} = Beam_RofC;
    case 3
        varargout{1} = Beam_rad;
        varargout{2} = Beam_RofC;
        varargout{3} = Beam_w0;
    case 3
        varargout{1} = Beam_rad;
        varargout{2} = Beam_RofC;
        varargout{3} = Beam_w0;
        varargout{4} = Beam_z2;
    otherwise
        error('Fit_TEM00(): Too many output argument')
end

end
