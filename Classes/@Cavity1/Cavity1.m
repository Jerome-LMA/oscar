classdef Cavity1
    %Cavity1 the class used to define cavity with thin mirrors
    % To define a cavity use
    % C1 = Cavity1(I1,I2,L,E1) with C1 the cavity defined, I1 and I2 2
    % interfaces to represent the mirrors of the cavity, L to represents
    % the length of the cavity in m and E1 the input field (an instance of
    % the class E_field)

    properties
        I_input
        I_end

        Length
        Laser_in
        Laser_start_on_input = false;
        Run_on_GPU = false;
        Resonance_phase = [];
        Cavity_scan_all_field = [];
        Cavity_scan_param = [1000 500 2E-9 1000]; % Number of points for the scan over one FSR, Number of points for the zoom, span of the zoom, max number of iteration (if the cavity is high finesse)
        Cavity_phase_param = 200;  % Number of iteration to find the resonance phase of the cavity, 100 is usually enough
        Cavity_scan_R = [];
        Cavity_scan_RZ = [];
        Cavity_EM_mat = [];

        Propagation_mat;     % Pre-compute the complex matrix used for the propagation
        ABCD_RT_mat;

        Field_circ = [];
        Field_ref = [];
        Field_trans = [];
        Field_reso_guess = [];
        Loss_RTL = [];

    end

    properties (Dependent)
        Wavelength
        Grid
    end

    methods
        function C = Cavity1(varargin)
            switch nargin
                case{0,1,2,3}
                    error('Cavity1(): at least 4 arguments must be given: 2 interfaces, a length and the input laser beam')
                case 4
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


                    C.I_input = varargin{1};
                    C.I_end = varargin{2};
                    C.Length = varargin{3};
                    C.Laser_in = varargin{4};

                    C.Propagation_mat = Prop_operator(C.Laser_in,C.Length);

                    % Check if the input laser is optimally mode matched
                    if  C.Laser_in.Optimal_mode_matching
                        Beam_paramater = Check_Stability(C,'Display',false);
                        if isempty(Beam_paramater)
                            error('Cavity1(): please check that the cavity is stable to find the optimal input beam parameters')
                        end
                        New_input_field =  E_Field(Grid(C),'w',Beam_paramater(1),'R',Beam_paramater(2),'mode',C.Laser_in.Mode_name,'Wavelength',C.Laser_in.Wavelength);
                        % add the sidebands as it used to be
                        for ii = 1:C.Laser_in.Nb_Pair_SB
                            New_input_field = Add_Sidebands(New_input_field,'Mod_freq',C.Laser_in.SB(ii).Frequency_Offset,'Mod_depth',C.Laser_in.SB(ii).Input_Mod_index);
                        end
                        C.Laser_in = New_input_field;
                    end

                    % Calculate the ABCD round trip matrix starting from
                    % the input mirror and going toward the end mirror
                    C.ABCD_RT_mat = C.I_input.ABCD_ref_from_n1 * C.Propagation_mat.ABCD_mat * C.I_end.ABCD_ref_from_n1 * C.Propagation_mat.ABCD_mat;

                otherwise
                    disp('Cavity1(): invalid number of input arguments, cavity not created')

            end
        end

        function value = get.Wavelength(obj)
            value = obj.Laser_in.Wavelength;
        end

        function value = get.Grid(obj)
            value = obj.Laser_in.Grid;
        end

    end

end
