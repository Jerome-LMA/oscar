function [Eout] = Reflect_Mirror(E_in,R_or_I,varargin)
%Reflect_mirror(E_Field,RofC) Reflect an E_Field from a spherical mirror
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Reflect_mirror(E_in,Rmir) reflect the field on a spherical mirror of
%     radius Rmir, Rmir > 0 for convergent mirror
%     Reflect_mirror(E_in,I1) pass the field E_in through an interface
%     I1 with the wavefront distortion given by k_prop * n * surface
%     Reflect_mirror(E_in,I1,'Ref',1) same as above but overwrite the
%     reflectivity of the mirror
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


p  = inputParser;
p.FunctionName = 'Function to reflect E_Field on mirror or interface';

% Check if the first argument is a grid
p.addRequired('E_in', @(x)isa(x, 'E_Field'));                               % Check if the first argument is an E_field
p.addRequired('R_or_I', @(x)isscalar(x) || isa(x, 'Interface') || isa(x, 'Mirror') );      % Check if the second argument is a number or an interface

p.addParameter('Ref',[],@(x)isreal(x) && x>=0 && x <=1); % As an option bypass the reflectivity of the surface

p.parse(E_in,R_or_I,varargin{:});

%p.Results

if isreal(R_or_I)
    RofC = R_or_I;
    Eout = E_in;
    
    if  isempty(p.Results.Ref)          % Check if the user has entered a reflectivity
        Reflectivity = 1;
    else
        Reflectivity = sqrt(p.Results.Ref);
    end
    
    
    % Defined the wavefront induced by the mirror
    WF_Mirror = (RofC - sign(RofC)*sqrt(RofC^2 - E_in.Grid.D2_square))*2;
    PF_Mirror =  exp(1i * WF_Mirror *E_in.k_prop);                                       % phase change induced by the mirror
    
    Eout.Field = E_in.Field .* PF_Mirror * Reflectivity;
    Eout.Field = fliplr(Eout.Field);
    
    if E_in.Nb_Pair_SB
        for ii=1:E_in.Nb_Pair_SB
            Eout.SB(ii).Field_lower = Eout.SB(ii).Field_lower  .* PF_Mirror * Reflectivity;
            Eout.SB(ii).Field_lower = fliplr(Eout.SB(ii).Field_lower);
            
            Eout.SB(ii).Field_upper = Eout.SB(ii).Field_upper .* PF_Mirror * Reflectivity;
            Eout.SB(ii).Field_upper = fliplr(Eout.SB(ii).Field_upper);
        end
    end
    
elseif isa(R_or_I, 'Interface')
    I =R_or_I;
    % consider the surface of the interface as the mirror
    
    if  isempty(p.Results.Ref)          % Check if the user has entered a reflectivity
        Reflectivity =  I.r;
    else
        Reflectivity = sqrt(p.Results.Ref);
    end
    
    if I.Run_on_GPU
        if (E_in.Refractive_index == I.n1)
            
            PF_Mirror_ref = I.WP_n1_GPU;
            
            if  ~isempty(p.Results.Ref)
                PF_Mirror_ref = PF_Mirror_ref * sqrt(p.Results.Ref)/ max(abs(PF_Mirror_ref(:)));
            end
            
        else
            
            PF_Mirror_ref = I.WP_n2_GPU;
            
            if  ~isempty(p.Results.Ref)
                PF_Mirror_ref = PF_Mirror_ref * sqrt(p.Results.Ref)/ max(abs(PF_Mirror_ref(:)));
            end
        end
    else
        if (E_in.Refractive_index==I.n1)
            PF_Mirror_ref = exp(-1i * E_in.k_prop * I.surface *2)  .* I.mask .* Reflectivity;
        else
            PF_Mirror_ref = exp(1i * E_in.k_prop * I.surface *2)  .* I.mask .* Reflectivity;
        end
    end
    
    %     figure(102)
    %     imagesc(angle(PF_Mirror_ref))
    
    Eout = E_in;
    if I.Run_on_GPU
        Eout.Field = arrayfun(@times,E_in.Field,PF_Mirror_ref);
    else
        Eout.Field = E_in.Field .* PF_Mirror_ref;
    end
    
    Eout.Field = fliplr(Eout.Field);
    
    if E_in.Nb_Pair_SB
        for ii=1:E_in.Nb_Pair_SB
            Eout.SB(ii).Field_lower = Eout.SB(ii).Field_lower  .* PF_Mirror_ref;
            Eout.SB(ii).Field_lower = fliplr(Eout.SB(ii).Field_lower);
            
            Eout.SB(ii).Field_upper = Eout.SB(ii).Field_upper .* PF_Mirror_ref;
            Eout.SB(ii).Field_upper = fliplr(Eout.SB(ii).Field_upper);
        end
    end
    
elseif isa(R_or_I, 'Mirror')
    
    Mir = R_or_I;
    % consider the reflection over a an object of class mirror. The
    % reflection is made on the HR side
    
    if  isempty(p.Results.Ref)          % Check if the user has entered a reflectivity
        [~,Eout] = Transmit_Reflect_Mirror(E_in,Mir,'HR');
    else
        Eout = Reflect_mirror(E_in,Mir.I_HR,'Ref',p.Results.Ref);
    end
    
else
    disp('Reflect_mirror(): The second argument must be a radius of curvature or an interface object')
end

end
