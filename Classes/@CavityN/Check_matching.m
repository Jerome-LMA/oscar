function Check_matching(varargin)
% Check_matching(C1) Do 2 round trip of the input field in the cavity and
% check the size and wavefront curvature of the beam on the mirrors to see
% if the matching is correct.
% Check_matching(C1,n) Do n round trip in the cavity



switch nargin
    case 0
        disp('Check_matching(): not enough arguments, at least an object CavityN must be given')
        return
        
    case 1
        if isa(varargin{1}, 'CavityN')
            Cin = varargin{1};
            Num_iter = 2;
        else
            disp('Check_matching(): The first argument must be an instance of CavityN')
            return
        end
        
    case 2
        if ~isa(varargin{1}, 'CavityN')
            disp('Check_matching(): The first argument must be an instance of CavityN')
            return
        end
        
        if ~real(varargin{2})
            disp('Check_matching(): if 2 arguments, the second one must be a number of iteration')
            return
        end
        Cin = varargin{1};
        Num_iter = varargin{2};
        
    otherwise
        disp('Check_matching(): Invalid number of input arguments')
        return
end

if ~Cin.Laser_start_on_input
    Cin.Laser_in = Change_E_n(Cin.Laser_in,Cin.I_array(1).n2);
    Cin.Laser_in = Transmit_Reflect_Interface(Cin.Laser_in,Cin.I_array(1));
end

Field_Circ = Cin.Laser_in;

for ii=1:Num_iter
    fprintf(' \n Round trip number: %i  \n',ii)
    
    if Cin.type == 'ring'
        for pp=1:Cin.Nb_mirror
            [Beam_rad, Beam_RofC] = Fit_TEM00(Field_Circ);
            fprintf('After the mirror %i,  beam radius [mm]: %7.4f \t wavefront RofC [m]: %5.2e \n',pp,Beam_rad*1E3,Beam_RofC)
            
            Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat_array(:,:,pp));
            %       [Beam_rad Beam_RofC] = Fit_TEM00(Field_Circ);
            %      fprintf('Before the mirror %i,   beam radius [m]: %7.4f \t wavefront RofC [m]: %5.2e \n',pp+1,Beam_rad,Beam_RofC)
            
            if pp ~= Cin.Nb_mirror % check we are not at the last iteration
                [Beam_rad, Beam_RofC] = Fit_TEM00(Field_Circ);
                fprintf('Before the mirror %i,   beam radius [mm]: %7.4f \t wavefront RofC [m]: %5.2e \n',pp+1,Beam_rad*1E3,Beam_RofC)
                Field_Circ = Reflect_mirror(Field_Circ,Cin.I_array(pp+1));
            end
            
        end
        
        fprintf('Before the mirror %i,   beam radius [mm]: %7.4f \t wavefront RofC [m]: %5.2e \n',1,Beam_rad*1E3,Beam_RofC)
        Field_Circ = Reflect_mirror(Field_Circ,Cin.I_array(1));
        
    elseif Cin.type == 'folded'
        for pp = 1:Cin.Nb_mirror-1 % do one way
            [Beam_rad, Beam_RofC] = Fit_TEM00(Field_Circ);
            fprintf('After the mirror %i,  beam radius [mm]: %7.4f \t wavefront RofC [m]: %5.2e \n',pp,Beam_rad*1E3,Beam_RofC)
            
            Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat_array(pp));
            Field_Circ = Reflect_mirror(Field_Circ,Cin.I_array(pp+1));
        end
        
        [Beam_rad, Beam_RofC] = Fit_TEM00(Field_Circ);
        fprintf('After the mirror %i,  beam radius [mm]: %7.4f \t wavefront RofC [m]: %5.2e \n',pp+1,Beam_rad*1E3,Beam_RofC)
        
        for pp=Cin.Nb_mirror-1:-1:1 % and do the round trip    
            Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat_array(pp));
            Field_Circ = Reflect_mirror(Field_Circ,Cin.I_array(pp));
            
            [Beam_rad, Beam_RofC] = Fit_TEM00(Field_Circ);
            fprintf('After the mirror %i,  beam radius [mm]: %7.4f \t wavefront RofC [m]: %5.2e \n',pp,Beam_rad*1E3,Beam_RofC)
        end
    end
    
    
    
end

end










