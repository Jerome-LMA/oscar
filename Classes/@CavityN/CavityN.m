classdef CavityN
    %CavityN the class used to define a cavity with an arbitrary number of
    %mirrors
    % C1 = Cavity1(Iarray,darray,E1) with C1 the cavity defined, Iarray an
    % array of interfaces, darray an array of the distance between the
    % interfaces.
    % interfaces to represent the mirrors of the cavity and E1 the input field (an instance of
    % the class E_field)
    
    properties
        I_array
        d_array
        
        Nb_mirror
        type                                    % kind of cavity linear folded or ring cavity
        laser_in
        laser_start_on_input = false
        resonance_phase = [];
        cavity_scan_all_field = [];
        cavity_scan_param = [1000 500 2E-9]; % Number of points for the scan over one  FSR, Number of points for the zoom, span of the zoom
        Cavity_phase_param = 200;  % Number of iteration to find the resonance phase of the cavity, 100 is usually enough
        cavity_scan_r = [];
        cavity_scan_rz = [];
        Cavity_EM_mat = [];
        
        propagation_mat_array;     % Pre-compute the complex matrix used for the propagation
        Field_reso_guess = [];       % Smart guess for the shape and amplitude of the fields
        
        Field_circ = [];
        Field_ref = [];
        Field_trans = [];
    end
    
    methods
        function C = CavityN(varargin)
            switch nargin
                case{0,1,2}
                    disp('Cavity1(): at least 3 arguments must be given: one array of interface, one array of distances and the input laser beam')
                    return
                case 3
                    if  ~isa(varargin{1}, 'Interface')
                        disp('Cavity1(): the first argument must be an instance of the class Interface')
                        return
                    end
                    
                    if  (~isreal(varargin{2})) && (varargin{2} <= 0)
                        disp('Cavity1(): the second argument, length of the cavity must be a real positive number')
                        return
                    end
                    
                    if  ~isa(varargin{3}, 'E_Field')
                        disp('Cavity1(): the third argument, the input laser beam must be an instance of the class E_field')
                        return
                    end
                    
                    
                    C.I_array = varargin{1};
                    C.d_array = varargin{2};
                    C.laser_in = varargin{3};
                    
                    if length(C.I_array) == length(C.d_array)
                        C.type = "ring";
                    elseif length(C.I_array) == (length(C.d_array)+1)
                        C.type = "folded";
                    else
                        error('CavityN(): number of interfaces and length are different')
                    end
                    
                    C.Nb_mirror = length(C.I_array);
                    
                    % Pre-allocate the propagation operator
                    
                    C.propagation_mat_array = Prop_operator.empty(0,C.Nb_mirror);
                    for pp = 1:length(C.d_array)
                        C.propagation_mat_array(pp) = Prop_operator(C.laser_in,C.d_array(pp));
                    end
                    
                    fprintf('You have created a %s cavity with %i mirrors \n',C.type,C.Nb_mirror)
                    
                    % Check if the input laser is optimally mode matched
                    if  C.laser_in.Optimal_mode_matching
                        Beam_paramater = check_stability(C,'Display',false);
                        New_input_field =  E_Field(Grid(C),'w',Beam_paramater(1),'R',Beam_paramater(2),'mode',C.laser_in.Mode_name);
                        % add the sidebands as it used to be
                        for ii = 1:C.laser_in.Nb_Pair_SB
                            New_input_field = Add_Sidebands(New_input_field,C.laser_in.SB(ii).Frequency_Offset,C.laser_in.SB(ii).Input_Mod_index);
                        end
                        C.laser_in = New_input_field;
                    end
                    
                otherwise
                    disp('Cavity1(): invalid number of input arguments, cavity not created')
                    
            end
        end
        
        function G = Grid(obj)
            G = obj.laser_in.Grid;
        end
    end
    
end
