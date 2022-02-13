classdef Cavity1 < handle
    %Cavity1 the class used to define cavity with thin mirrors
    % To define a cavity use
    % C1 = Cavity1(I1,I2,L,E1) with C1 the cavity defined, I1 and I2 2
    % interfaces to represent the mirrors of the cavity, L to represents
    % the length of the cavity in m and E1 the input field (an instance of
    % the class E_field)
    
    properties
        i_input
        i_end
        
        length
        laser_in
        laser_start_on_input = false;
        run_on_gpu = false;
        resonance_phase = [];
        cavity_scan_all_field = [];
        cavity_scan_param = [1000 500 2E-9 1000]; % Number of points for the scan over one FSR, Number of points for the zoom, span of the zoom, max number of iteration (if the cavity is high finesse)
        cavity_phase_param = 200;  % Number of iteration to find the resonance phase of the cavity, 100 is usually enough
        cavity_scan_r = [];
        cavity_scan_rz = [];
        cavity_em_mat = [];
        
        propagation_mat;     % Pre-compute the complex matrix used for the propagation
        
        field_circ = [];
        field_ref = [];
        field_trans = [];
        field_reso_guess = [];
        loss_rtl = [];                
        
    end
    
    methods
        % Declare methods which are implemented in external files
        [] = declare_on_gpu(obj)
        [power_buildup] = calculate_fields(obj, varargin)
        [] = calculate_fields_ac(obj, varargin)
        [length_scan, power_scan] = scan(obj, varargin)        
    end
    
    methods                
        function obj = Cavity1(varargin)
            switch nargin
                case{0,1,2,3}
                    error('Cavity1(): at least 4 arguments must be given: 2 interfaces, a length and the input laser beam')
                case {4,5}
                    if  ~(isa(varargin{1}, 'Interface') || isa(varargin{1}, 'Mirror'))
                        error('Cavity1(): the first argument must be an instance of the class Interface or Mirror')
                    end
                    
                    if  ~(isa(varargin{2}, 'Interface') || isa(varargin{2}, 'Mirror'))
                        error('Cavity1(): the second argument must be an instance of the class Interface or Mirror')
                    end
                    
                    if  (~isreal(varargin{3})) && (varargin{3} <= 0)
                        error('Cavity1(): the third argument, length of the cavity must be a real positive number')
                    end
                    
                    if  ~isa(varargin{4}, 'E_Field')
                        error('Cavity1(): the fourth argument, the input laser beam must be an instance of the class E_field')
                    end                                        
                    
                    obj.i_input = varargin{1};
                    obj.i_end = varargin{2};
                    obj.length = varargin{3};
                    obj.laser_in = varargin{4};
                    
                    obj.propagation_mat = Prop_operator(obj.laser_in, obj.length);
                 
                % Check if the input laser is optimally mode matched    
                if  obj.laser_in.Optimal_mode_matching
                    Beam_paramater = check_stability(obj,'Display',false);
                    New_input_field =  E_Field(Grid(obj),'w',Beam_paramater(1),'R',Beam_paramater(2),'mode',obj.laser_in.Mode_name);
                    % add the sidebands as it used to be
                    for ii = 1:obj.laser_in.Nb_Pair_SB
                        New_input_field = Add_Sidebands(New_input_field,'Mod_freq',obj.laser_in.SB(ii).Frequency_Offset,'Mod_index',obj.laser_in.SB(ii).Input_Mod_index);
                    end
                    obj.laser_in = New_input_field;
                end    
                
                if (nargin == 5)
                    if  ~islogical(varargin{5})
                        error('Cavity1(): the fifth optional argument must be bool and specifies, weather or not to use the GPU (defaults to false)')
                    end   
                    obj.declare_on_gpu()                    
                end
                    
                otherwise
                    disp('Cavity1(): invalid number of input arguments, cavity not created')
                   
            end
        end
                        
        % TODO: Implement somewhere else
        function G = Grid(obj)
            G = obj.laser_in.Grid;
        end                

        
    end
    
end
