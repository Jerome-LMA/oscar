function Check_matching(varargin)
% Check_matching(C1) Do 2 round trip of the input field in the cavity and
% check the size and wavefront curvature of the beam on the mirrors to see
% if the matching is correct.
% Check_matching(C1,n) Do n round trip in the cavity 



switch nargin
    case 0
        disp('Check_matching(): not enough arguments, at least an object Cavity1 must be given')
        return
        
    case 1
        if isa(varargin{1}, 'Cavity1')
            Cin = varargin{1};
            Num_iter = 2;
        else
            disp('Check_matching(): The first argument must be an instance of Cavity1')
            return
        end
        
    case 2
        if ~isa(varargin{1}, 'Cavity1')
            disp('Check_matching(): The first argument must be an instance of Cavity1')
            return
        end
        
        if ~real(varargin{2})
            disp('Check_matching(): if 2 arguments, the second one must be a number of iteration')
            return
        end
        Cin = varargin{1};
        Num_iter = varargin{2};
        
    otherwise
        disp('Calculate_power(): Invalid number of input arguments, no power calculation is made')
        return
end

if ~Cin.Laser_start_on_input
    if isa(Cin.I_input, 'Interface')
        Cin.Laser_in =  Change_E_n(Cin.Laser_in,Cin.I_input.n2);
    end
    Cin.Laser_in = Transmit_Reflect_Optic(Cin.Laser_in,Cin.I_input);
end

Field_Circ = Cin.Laser_in;

for ii=1:Num_iter
    fprintf(' \n Round trip number: %i  \n',ii) 
    [Beam_rad,Beam_RofC] = Fit_TEM00(Field_Circ);    
    fprintf('After the input mirror,  beam radius [m]: %7.4f \t wavefront RofC [m]: %5.2e \n',Beam_rad,Beam_RofC)
    
    Field_Circ = Propagate_E(Field_Circ,Cin.Length);
    [Beam_rad,Beam_RofC] = Fit_TEM00(Field_Circ);    
    fprintf('Before the end mirror,   beam radius [m]: %7.4f \t wavefront RofC [m]: %5.2e \n',Beam_rad,Beam_RofC)     
    
    Field_Circ = Reflect_mirror(Field_Circ,Cin.I_end);
    [Beam_rad,Beam_RofC] = Fit_TEM00(Field_Circ);    
    fprintf('After the end mirror,    beam radius [m]: %7.4f \t wavefront RofC [m]: %5.2e \n',Beam_rad,Beam_RofC) 
  
    Field_Circ = Propagate_E(Field_Circ,Cin.Length);
    [Beam_rad,Beam_RofC] = Fit_TEM00(Field_Circ);    
    fprintf('Before the input mirror, beam radius [m]: %7.4f \t wavefront RofC [m]: %5.2e \n',Beam_rad,Beam_RofC)
    
    Field_Circ = Reflect_mirror(Field_Circ,Cin.I_input);   
end

end










