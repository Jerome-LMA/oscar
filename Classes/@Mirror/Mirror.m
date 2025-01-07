classdef Mirror
    %MIRROR Define the class mirrors which have 2 surfaces and a
    %substrate

    properties
        I_HR
        I_AR
        length_substrate

        RT_inside = 1;             % Number of round trip when the thick mirror is dealt as a cavity
        n_substrate
        Propagation_mat_sub
    end

    properties (Dependent)
        r
        t
        R
        T
        mask
        Grid
        surface
        ABCD_ref_from_n1  % will assume it from the HR surface

    end

    methods
        function M = Mirror(varargin)

            switch nargin
                case {0,1,2}
                    error('Mirror(): at least 3 arguments must be given: 2 interfaces and the length of the substrate')

                case 3
                    if  ~isa(varargin{1}, 'Interface')
                        error('Mirror(): the first argument must be an instance of the class interface')
                    end

                    if  ~isa(varargin{2}, 'Interface')
                        error('Mirror(): the second argument must be an instance of the class interface')
                    end

                    if  ~isscalar(varargin{3})
                        error('Mirror(): the third argument must be a scalar (the length of the substrate)')
                    end

                    if ( max(varargin{1}.n1,varargin{1}.n2) ~=  max(varargin{1}.n1,varargin{1}.n2) )
                        error('Mirror(): Inconsistencie inthe refractive index of the 2 surfaces')
                    end

                    M.I_HR = varargin{1};
                    M.I_AR = varargin{2};

                    if (M.I_HR.T  >  M.I_AR.T )
                        warning('Mirror(): Transmission of the HR surface (%g) higher than the AR one (%g)',M.I_HR.T , M.I_AR.T)
                    end

                    M.length_substrate =  varargin{3};
                    M.n_substrate = max(varargin{1}.n1,varargin{1}.n2);

                otherwise
                    disp('Mirror(): invalid number of input arguments, the mirror is not created')
                    return

            end
        end

        function value = get.r(obj)
            value = obj.I_HR.r;
        end

        function value = get.t(obj)
            value = obj.I_HR.t;
        end

        function value = get.R(obj)
            value = (obj.I_HR.r)^2;
        end

        function value = get.T(obj)
            value = (obj.I_HR.t)^2;
        end

        function value = get.mask(obj)
            value = obj.I_HR.mask;
        end

        function value = get.Grid(obj)
            value = obj.I_HR.Grid;
        end

        function value = get.surface(obj)
            value = obj.I_HR.surface;
        end

        function value = get.ABCD_ref_from_n1(obj)
            value = obj.I_HR.ABCD_ref_from_n1;
        end

    end

end

