classdef Mirror
    %MIRROR Define the class mirrors which have 2 surfaces and a
    %substrate.
    % Input the HR surface, then Ar surface and then the thickness
    
    properties
        I_HR
        I_AR
        length_substrate
        
        Propagation_mat_sub
        
        RT_inside = 1;             % Number of round trip when the thick mirror is dealt as a cavity
        n_substrate
        
        r
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
                        error('Mirror(): Inconsistencie index refractive index of the 2 surfaces')
                    end
                    
                    M.I_HR = varargin{1};
                    M.I_AR = varargin{2};
                    
                    if (M.I_HR.T  >  M.I_AR.T )
                        warning('Mirror(): Transmission of the HR surface (%g) higher than the AR one (%g)',M.I_HR.T , M.I_AR.T)
                    end
                    
                    M.length_substrate =  varargin{3};
                    M.n_substrate = max(varargin{1}.n1,varargin{1}.n2);
                    
                    M.r = M.I_HR.r;
                    
                    M.Propagation_mat_sub = M.length_substrate; % Could be overwrite if the propagation matrix is preallocated (absolutely recommended) 
                    
                otherwise
                    disp('Mirror(): invalid number of input arguments, the mirror is not created')
                    return
                    
            end
        end
        
        
        
    end
    
    
    
    
end

