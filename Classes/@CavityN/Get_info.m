function Get_Info(Cin)

% get_info(C1) calculate some informations about the cavity
% get_info(C1) calculate fields in the cavity, power buildup, diffraction loss
% Do not calculate the reflected field use the function 'Calculate_fields' instead.

if isempty(Cin.resonance_phase)
    error(['get_info(' inputname(1) '): The resonance position must be calculated first'])
end

if ~Cin.laser_start_on_input
    Field_in =  Change_E_n(Cin.laser_in,Cin.I_array(1).n2);
    Field_in = Transmit_Reflect_Interface(Field_in,Cin.I_array(1));
    Field_Circ = Field_in;
else
    Field_Circ = Cin.laser_in * Cin.I_array(1).t;
end
% Cin.I_array(1)
% 
% plot(Field_Circ); pause
% Fit_TEM00(Field_Circ)

% Calculate the number of iteration to reach the steady state for a given
% accurary
Accuracy = 0.0001;

RT_loss = 1;
for pp=1:Cin.Nb_mirror
    RT_loss = RT_loss * Cin.I_array(pp).r;
end

% Have to solve RT_loss^num_iter < 0.5*accuracy
num_iter = log(0.5*Accuracy)/(log(RT_loss));
num_iter = round(num_iter);
%num_iter = 2000;

Power_buildup = zeros(1,num_iter,'double');
%Cin.laser_in = Normalise_E(Cin.laser_in);

field_transient = Field_Circ;
Field_total = Normalise_E(Cin.laser_in,0);

for q = 1:num_iter
    
    Field_total = Field_total + field_transient;
    Power_buildup(q) = Calculate_power(Field_total);
    
    for pp=1:Cin.Nb_mirror
        if pp ~= Cin.Nb_mirror % check we are not at the last iteration
            field_transient = Propagate_E(field_transient,Cin.propagation_mat_array(pp));
            field_transient = reflect_mirror(field_transient,Cin.I_array(pp+1));
        else
            field_transient = Propagate_E(field_transient,Cin.propagation_mat_array(pp)) * Cin.resonance_phase;
            field_transient = reflect_mirror(field_transient,Cin.I_array(1));
        end
    end
    
end


% Calculate the diffraction loss
field_tmp = Field_total;
field_tmp = Normalise_E(field_tmp);

for pp=1:Cin.Nb_mirror
    if pp ~= Cin.Nb_mirror % check we are not at the last iteration
        field_tmp = Propagate_E(field_tmp,Cin.propagation_mat_array(pp));
        field_tmp = reflect_mirror(field_tmp,Cin.I_array(pp+1),'Ref',1);
    else
        field_tmp = Propagate_E(field_tmp,Cin.propagation_mat_array(pp)) * Cin.resonance_phase;
        field_tmp = reflect_mirror(field_tmp,Cin.I_array(1),'Ref',1);
    end
end

Cavity_loss = (1 - Calculate_power(field_tmp));


disp([' ---- Display results for cavity ' inputname(1)   ' -----'])
fprintf(' Round trip diffraction loss: %g [ppm] \n',Cavity_loss*1E6)

fprintf(' Circulating power: %g [W] \n',Calculate_power(Field_total))

% Find if the mode is a TEM00 or not
[~, m n] = Read_mode_name(Cin.laser_in.Mode_name);

if (m==0) && (n==0)
    
    [Beam_radius, Beam_RofC] = Fit_TEM00(Field_total);
    fprintf(' Size of the beam on the input mirror: %g [m] \n',Beam_radius)
end


figure(104)
clf;
subplot(2,2,1)
plot(Cin.laser_in)
title('Input field')
subplot(2,2,2)
plot(Field_total)
title('Circulating field')
subplot(2,2,[3 4])
plot(Power_buildup)
title('Power buildup')
xlabel('Number of iteration')
ylabel('Power [W]')



end