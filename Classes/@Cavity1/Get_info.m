function varargout = get_info(obj)

% get_info(C1) calculate some informations about the cavity
% get_info(C1) calculate fields in the cavity, power buildup, diffraction loss
% Do not calculate the reflected field use the function 'Calculate_fields' instead.

if isempty(obj.resonance_phase)
    error(['get_info(' inputname(1) '): The resonance position must be calculated first'])
end

if ~obj.laser_start_on_input
    
    if isa(obj.i_input, 'Interface')
        obj.laser_in =  Change_E_n(obj.laser_in,obj.i_input.n2);
    else
        obj.laser_in =  obj.laser_in;
    end
    
    [obj.laser_in,~] = Transmit_Reflect_Optic(obj.laser_in,obj.i_input,'AR');
    
else
    obj.laser_in =  obj.laser_in * obj.i_input.t;
end

%Calculate_power(Cin.laser_in);

% Calculate the number of iteration to reach the steady state for a given
% accurary
Accuracy = 0.0001;

RT_loss = obj.i_input.r*obj.i_end.r;
% Have to solve RT_loss^num_iter < 0.5*accuracy
num_iter = log(0.5*Accuracy)/(log(RT_loss));
num_iter = round(num_iter);
%num_iter = 2000;

Power_buildup = zeros(1,num_iter,'double');
%Cin.laser_in = Normalise_E(Cin.laser_in);

field_transient = obj.laser_in;
Field_total = Normalise_E(obj.laser_in,0);



for q = 1:num_iter
    Field_total = Field_total + field_transient;
    Power_buildup(q) = Calculate_power(Field_total);
    field_transient = Propagate_E(field_transient,obj.propagation_mat);
    field_transient = reflect_mirror(field_transient,obj.i_end);
    field_transient = Propagate_E(field_transient,obj.propagation_mat);
    field_transient = reflect_mirror(field_transient,obj.i_input);
    field_transient = field_transient * obj.resonance_phase;
end


% Calculate the diffraction loss
field_tmp = Field_total;
field_tmp = Normalise_E(field_tmp);
field_tmp = Propagate_E(field_tmp,obj.propagation_mat);
field_tmp = reflect_mirror(field_tmp,obj.i_end,'Ref',1);
field_tmp = Propagate_E(field_tmp,obj.propagation_mat);
field_tmp = reflect_mirror(field_tmp,obj.i_input,'Ref',1);

Cavity_loss = (1 - Calculate_power(field_tmp));


disp([' ---- Display results for cavity ' inputname(1)   ' -----'])
fprintf(' Round trip diffraction loss: %g [ppm] \n',Cavity_loss*1E6)

fprintf(' Circulating power: %g [W] \n',Calculate_power(Field_total))

% Find if the mode is a TEM00 or not
[~, m n] = Read_mode_name(obj.laser_in.Mode_name);

if (m==0) && (n==0)
    [Beam_radius, ~] = Fit_TEM00(Propagate_E(Field_total,obj.Length));
    fprintf(' Size of the beam on the end mirror: %g [m] \n',Beam_radius)
    
    [Beam_radius, Beam_RofC] = Fit_TEM00(Field_total);
    fprintf(' Size of the beam on the input mirror: %g [m] \n',Beam_radius)
    
    
    % Calculate the q parameter of the beam on the input mirror
    Cavity_q = 1/(1/Beam_RofC - 1i*obj.laser_in.Wavelength/(pi*Beam_radius^2));
    Cavity_waist_size = sqrt(( (imag(Cavity_q))*obj.laser_in.Wavelength/pi ));
    Arm.waist_position = real(Cavity_q);
    
    fprintf(' Size of the cavity waist: %g [m] \n',Cavity_waist_size)
    fprintf(' Distance of the cavity waist from the input mirror: %g [m] \n',Arm.waist_position)
end


figure(104)
clf;
subplot(2,2,1)
plot(obj.laser_in)
title('Input field')
subplot(2,2,2)
plot(Field_total)
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
        error('get_info(): Too many output argument')
end

end