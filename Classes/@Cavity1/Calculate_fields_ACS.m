function Cout = Calculate_fields_ACS(Cin,varargin)
% Cout = Calculate_fields_AC(Cin) calculate the circulating, reflected and transmitted fields
% Function used to calculated the fields inside the cavity. The laser beam must be defined outside the cavity in order to calculate the reflected field.
% Use the accelerated convergence scheme

p  = inputParser;
p.FunctionName = 'Calculate fields inside the cavity';

% Check if the first argument is an interface
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

% Accelerated convergence only works for the carrier. Remove the SB fields
Cin.Laser_in.Field_SBl=[];
Cin.Laser_in.Field_SBu=[];
Cin.Field_reso_guess.Field_SBl=[];
Cin.Field_reso_guess.Field_SBu=[];

Cout = Cin;

if ~isempty(p.Results.accuracy)
    Accuracy = p.Results.accuracy;
else
    Accuracy = 1E-12;
end

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

ii = 0;
error_P = 1;
E1 = Cin.Field_reso_guess * sqrt(Calculate_power(Cin.Laser_in)); % The 'guess' field was calculated for 1W of input power, so it has to be normalised according to the current incident power

Nb_RT_smooth = 2;

% Start with E1, Calculate D1 = E1 - A E1
En(1) = E1;

E1_circ = Propagate_E(E1,Cin.Propagation_mat);
E1_circ = Reflect_mirror(E1_circ,Cin.I_end);
E1_circ = Propagate_E(E1_circ,Cin.Propagation_mat);

E1_circ = E1_circ * Cin.Resonance_phase;
E1_circ = Reflect_mirror(E1_circ,Cin.I_input);

Dn(1) = E1 -  E1_circ;

ii = 1;

while error_P > Accuracy;
%for ii = 1:10
    % ii = ii + 1
    for pp = 2:Nb_RT_smooth
        
        En(pp) = En(pp-1) - Dn(pp-1) + Field_in;
        E1 =  En(pp) ;
        
        E1_circ = Propagate_E(E1,Cin.Propagation_mat);
        E1_circ = Reflect_mirror(E1_circ,Cin.I_end);
        E1_circ = Propagate_E(E1_circ,Cin.Propagation_mat);
        
        E1_circ = E1_circ * Cin.Resonance_phase;
        E1_circ = Reflect_mirror(E1_circ,Cin.I_input);
        
        Dn(pp) = E1 -  E1_circ;
        
    end
    
    % Fill the matrices M and A
    M = zeros(Nb_RT_smooth);
    A = zeros(Nb_RT_smooth,1);
    
    for pp = 1:Nb_RT_smooth
        for oo = 1:pp
            M(pp,oo) = Raw_overlap(Dn(pp).Field,Dn(oo).Field);
            if oo ~= pp
                M(oo,pp) = conj(M(pp,oo));
            end
        end
        A(pp) = Raw_overlap(Dn(pp).Field,Field_in.Field);
    end
    
     c = M\A;
    
     [~,S,~] = svd(M);
     diag_S  = diag(S);
%     
%     if diag_S < 10
%         diag_S = inf;
%     end
%     
%     c = V*(diag(1./diag_S))*(U') *A;
    
    E2 = Normalise_E(E1,0);
    
    for pp=1:Nb_RT_smooth
        E2 = E2 + c(pp) * En(pp);
    end
    
    D2 = E2;
    for pp=1:Nb_RT_smooth
        D2 = D2 - c(pp) * (En(pp) - Dn(pp));
    end
    
    error_P = Calculate_power(E2 - En(1)) / Calculate_power(En(1));
    
    Nb_RT_smooth = sum(diag_S(1)./diag_S <2E10 )+1;
    
    En(1)  = E2;
    Dn(1)  = D2;
    
end

Field_total = E2;
Cout.Field_circ = Field_total;

%------------------------------------------------------------------
% Calculate the transmitted and reflected field

Field_temp = Propagate_E(Field_total,Cin.Propagation_mat);
Cout.Field_trans = Transmit_Reflect_Optic(Field_temp,Cin.I_end);

Field_temp = Propagate_E(Field_total,Cin.Propagation_mat);
Field_temp = Reflect_mirror(Field_temp,Cin.I_end);
Field_temp = Propagate_E(Field_temp,Cin.Propagation_mat);
Field_temp = Field_temp * Cin.Resonance_phase;

Field_temp = Transmit_Reflect_Optic(Field_temp,Cin.I_input);

Cout.Field_ref = Field_reflec + Field_temp;


if isa(Cin.I_input, 'Interface')
    Cout.Field_ref =  Change_E_n(Cout.Field_ref,Cin.I_input.n1);
end

if isa(Cin.I_end, 'Interface')
   Cout.Field_trans =  Change_E_n(Cout.Field_trans,Cin.I_end.n1);
end

end