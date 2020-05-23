function Cout = Calculate_fields_AC2(Cin,varargin)
% Cout = Calculate_fields_AC(Cin) calculate the circulating, reflected and transmitted fields
% Function used to calculated the fields inside the cavity. The laser beam must be defined outside the cavity in order to calculate the reflected field.
% Use the accelerated convergence scheme
% Based on "Accelerated convergence method for fast Fourier transform
% simulation of coupled cavities" from R. Day

p  = inputParser;
p.FunctionName = 'Calculate fields inside the cavity';

% Check if the first argument is an Cavity1 object
p.addRequired('Cin', @(x)isa(x, 'Cavity1'));

% Check if the resolution of the grid if given
p.addParameter('accuracy',[],@(x)isnumeric(x) && x>0);

% Check if the resolution of the grid if given
p.addParameter('iter',[],@(x)isnumeric(x) && x>0);

p.parse(Cin,varargin{:})

if isempty(Cin.Resonance_phase)
    error(['Calculate_fields(' inputname(1) '): The resonance position must be calculated first'])
end

if Cin.Laser_start_on_input
    error(['Calculate_fields(' inputname(1) '): To calculate the reflected beam, the beam must be defined outside the cavity, set Laser_start_on_input = false'])
end

Cout = Cin;

if ~isempty(p.Results.accuracy)
    Accuracy = p.Results.accuracy;
else
    Accuracy = 1E-12;
end

display_iter = false;

Nb_Average = 6;
Nb_Field_to_take = 6;

% Define the overlap function:
Raw_overlap = @(x,y) (sum(conj(x).*y, 'all') );

% The laser starts outside the input mirror, change n from 1 to mirror
% substrate refractive index

if isa(Cin.I_input, 'Interface')
    Field_in =  Change_E_n(Cin.Laser_in,Cin.I_input.n2);
else
    Field_in =  Cin.Laser_in;
end

[Field_in,Field_reflec] = Transmit_Reflect_Optic(Field_in,Cin.I_input,'AR');

ii = 1;
error_P(ii) = 1;
E_geuss = Cin.Field_reso_guess * sqrt(Calculate_power(Cin.Laser_in)); % The 'guess' field was calculated for 1W of input power, so it has to be normalised according to the current incident power

% Preallocate the matrices with dummy fields, suppose average is 2



E_1 = E_geuss;

while error_P(ii) > Accuracy;
    ii = ii +1;
    
    E_n(1:Nb_Field_to_take) = E_Field(Cin.Grid);
    D_n(1:Nb_Field_to_take) = E_Field(Cin.Grid);
    
    for n = 1:Nb_Field_to_take
        
        for m = 1:Nb_Average
            
            E_n(n) = E_n(n) + E_1;
            
            % Can save one FFT here!
            %             if n == 1   %% So we want to calculate P_1 for the first time
            %                 if ii == 2 % if first iteration for the big loop
            %                     P_1 = Do_RT(Cin,E_1);
            %                 else % We can calculate it from previous iteration (ii-1) to save one RT calculation
            %                     P_1 = c(1)*(E_1_previous - D_n(1)) + c(2)*(E_n(2) - D_n(2));
            %                 end
            %             else
            %                 P_1 = Do_RT(Cin,E_1);
            %             end
            
            P_1 = Do_RT(Cin,E_1);
            
%             P_1 = Do_RT(Cin.RT_fast,E_1);
%             P_1 = P_1 * Cin.Resonance_phase;
            
            D_n(n) = D_n(n) + (E_1 - P_1);
            E_1 = P_1 + Field_in;
            
        end
        D_n(n) = (1/Nb_Average) * D_n(n);
        E_n(n) = (1/Nb_Average) * E_n(n);
    end
    
    for pp = 1:Nb_Field_to_take
        for oo = 1:Nb_Field_to_take
            M(pp,oo) = Raw_overlap(D_n(pp).Field,D_n(oo).Field);
        end
    end
    
    for pp = 1:Nb_Field_to_take
        A(pp,1) = Raw_overlap(D_n(pp).Field,Field_in.Field);
    end
    
    % Classic inversion
    c = M\A;
    
    %    E_1_previous = E_n(1);
    E_3 = E_Field(Cin.Grid);
    
    for pp = 1:Nb_Field_to_take
        E_3 = E_3 + c(pp) * E_n(pp);
    end
    
    error_P(ii) = Calculate_power(E_3 - E_1) / Calculate_power(E_1);
    
    if display_iter
        if rem(ii,50) == 0
            fprintf('Nb iteration: %i, error: %g \n',ii,error_P(ii))
        end
    end
    %error_P(ii)
    E_1 = E_3;
    
    if ii > 200 % max iteration, break the loop
        fprintf('Number of iteration max reached, exiting now, error: %g \n',error_P(ii))
        break;
    end
    
    
end

Cout.Field_circ = E_3; % circulating field found

if display_iter
    figure(10); semilogy(error_P)
end

%---------------------------------------------------------------------------------------------------------------------------------
% Ok we got now the 3 circulating fields (carrier + 2 sidebands)

Field_temp = Propagate_E(Cout.Field_circ,Cin.Propagation_mat);
Cout.Field_trans = Transmit_Reflect_Optic(Field_temp,Cin.I_end);

Field_temp = Propagate_E(Cout.Field_circ,Cin.Propagation_mat);
Field_temp = Reflect_mirror(Field_temp,Cin.I_end);
Field_temp = Propagate_E(Field_temp,Cin.Propagation_mat);
Field_temp = Field_temp * Cin.Resonance_phase;

Field_temp = Transmit_Reflect_Optic(Field_temp,Cin.I_input);
Cout.Field_ref = Field_reflec + Field_temp;

% Calculate the RTL
field_tmp = Cout.Field_circ;
field_tmp = Normalise_E(field_tmp);
field_tmp = Propagate_E(field_tmp,Cin.Propagation_mat);
field_tmp = Reflect_mirror(field_tmp,Cin.I_end,'Ref',1);
field_tmp = Propagate_E(field_tmp,Cin.Propagation_mat);
field_tmp = Reflect_mirror(field_tmp,Cin.I_input,'Ref',1);

Cout.Loss_RTL = (1 - Calculate_power(field_tmp));
%

if isa(Cin.I_input, 'Interface')
    Cout.Field_ref =  Change_E_n(Cout.Field_ref,Cin.I_input.n1);
end

if isa(Cin.I_end, 'Interface')
   Cout.Field_trans =  Change_E_n(Cout.Field_trans,Cin.I_end.n1);
end


end