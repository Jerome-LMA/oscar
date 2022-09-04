function varargout = Get_info(Cin)

% Get_info(C1) calculate some informations about the cavity
% Get_info(C1) calculate fields in the cavity, power buildup, diffraction loss
% Do not calculate the reflected field use the function 'Calculate_fields' instead.

if isempty(Cin.Resonance_phase)
    error(['Get_info(' inputname(1) '): The resonance position must be calculated first'])
end

if ~Cin.Laser_start_on_input
    
    if isa(Cin.I_input, 'Interface')
        Cin.Laser_in =  Change_E_n(Cin.Laser_in,Cin.I_input.n2);
    else
        Cin.Laser_in =  Cin.Laser_in;
    end
    
    [Cin.Laser_in,~] = Transmit_Reflect_Optic(Cin.Laser_in,Cin.I_input,'AR');
    
else
    Cin.Laser_in =  Cin.Laser_in * Cin.I_input.t;
end

%Calculate_power(Cin.Laser_in);

% Calculate the number of iteration to reach the steady state for a given
% accurary
Accuracy = 0.0001;

RT_loss = Cin.I_input.r*Cin.I_end.r;
% Have to solve RT_loss^num_iter < 0.5*accuracy
num_iter = log(0.5*Accuracy)/(log(RT_loss));
num_iter = round(num_iter);
%num_iter = 2000;

Power_buildup = zeros(1,num_iter,'double');
%Cin.Laser_in = Normalise_E(Cin.Laser_in);

Field_transient = Cin.Laser_in;
Field_total = Normalise_E(Cin.Laser_in,0);



for q = 1:num_iter
    Field_total = Field_total + Field_transient;
    Power_buildup(q) = Calculate_power(Field_total);
    Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat);
    Field_transient = Reflect_mirror(Field_transient,Cin.I_end);
    Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat);
    Field_transient = Reflect_mirror(Field_transient,Cin.I_input);
    Field_transient = Field_transient * Cin.Resonance_phase;
end


% Calculate the diffraction loss
field_tmp = Field_total;
field_tmp = Normalise_E(field_tmp);
field_tmp = Propagate_E(field_tmp,Cin.Propagation_mat);
field_tmp = Reflect_mirror(field_tmp,Cin.I_end,'Ref',1);
field_tmp = Propagate_E(field_tmp,Cin.Propagation_mat);
field_tmp = Reflect_mirror(field_tmp,Cin.I_input,'Ref',1);

Cavity_loss = (1 - Calculate_power(field_tmp));


disp([' ---- Display results for cavity ' inputname(1)   ' -----'])
fprintf(' Round trip diffraction loss: %g [ppm] \n',Cavity_loss*1E6)

fprintf(' Circulating power: %g [W] \n',Calculate_power(Field_total))

% Find if the mode is a TEM00 or not
[~, m n] = Read_mode_name(Cin.Laser_in.Mode_name);

if (m==0) && (n==0)
    [Beam_radius, ~] = Fit_TEM00(Propagate_E(Field_total,Cin.Length));
    fprintf(' Size of the beam on the end mirror: %g [m] \n',Beam_radius)
    
    [Beam_radius, Beam_RofC] = Fit_TEM00(Field_total);
    fprintf(' Size of the beam on the input mirror: %g [m] \n',Beam_radius)
    
    
    % Calculate the q parameter of the beam on the input mirror
    Cavity_q = 1/(1/Beam_RofC - 1i*Cin.Laser_in.Wavelength/(pi*Beam_radius^2));
    Cavity_waist_size = sqrt(( (imag(Cavity_q))*Cin.Laser_in.Wavelength/pi ));
    Arm.waist_position = real(Cavity_q);
    
    fprintf(' Size of the cavity waist: %g [m] \n',Cavity_waist_size)
    fprintf(' Distance of the cavity waist from the input mirror: %g [m] \n',Arm.waist_position)
end


figure(104)
clf;
subplot(2,2,1)
E_plot(Cin.Laser_in)
title('Input field')
subplot(2,2,2)
E_plot(Field_total)
title('Circulating field')
subplot(2,2,[3 4])
plot(Power_buildup)
title('Power buildup')
xlabel('Number of iteration')
ylabel('Power [W]')


switch nargout
    case 0
    case 1
        varargout{1} = Cavity_loss;
    case 2
        varargout{1} = Cavity_loss;
        varargout{2} = Field_total;
    otherwise
        error('Get_info(): Too many output argument')
end

end