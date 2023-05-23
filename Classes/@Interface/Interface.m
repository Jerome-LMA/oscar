classdef Interface
    %     Interface() Create a interface between 2 media of different refractive
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Interface(Grid) create a flat interface
    %
    %     Interface(Grid,'RoC',2E3) create an interface with a radius of curvature
    %     of 2 km.
    %
    %     Interface(Grid,'RoC',2E3,'CA',0.1) create an interface with a radius of curvature
    %     of 2 km with a clear aperture of 10 cm (usually the diameter of the optics)
    %
    %     Interface(Grid,'RoC',2E3,'CA',0.1,'T',0.1,'L',50E-6) same as
    %     before but with transmission of 10% and loss 50 ppm in power.
    %
    %     Interface(Grid,'RoC',2E3,'CA',0.1,'T',0.1,'L',50E-6,'n1',1,'n2',1.45) same as
    %     before but with transmission of 10% and loss 50 ppm in power. To
    %     add the refractive index
    %
    %     RoC > 0 for a concave mirror (surface view from n1 toward n2)
    %     RoC < 0 for a convex mirror  (surface view from n1 toward n2)
    %
    %     If the interface is used as a lens with RofC < 0 and n2 > n1 it
    %     will be equivalent to a convergent lens
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        Grid
        
        surface
        mask
        T
        L
        n1
        n2

        WP_n1_GPU
        WP_n2_GPU
        
        t
        r
    end
    
    methods
        
        function I = Interface(Grid_in,varargin)
            p  = inputParser;
            p.FunctionName = 'Interface creator';
            
            % Check if the first argument is a grid
            p.addRequired('Grid_in', @(x)isa(x, 'Grid'));
            
            % enter the radius of curvature of the mirror
            p.addParamValue('RoC',1E99,@(x)isnumeric(x) );
            
            % enter the clear aperture, i.e. diameter of the optic
            p.addParamValue('CA',1E99,  @(x)isnumeric(x) && x>0);
            
            % enter the transmission in power
            p.addParamValue('T',0.1,  @(x)isnumeric(x) && x>=0 && x<=1);
            
            % enter the loss in power
            p.addParamValue('L', 0, @(x)isnumeric(x) && x>=0 && x<=1);
            
            % enter the refractive index of the first media
            p.addParamValue('n1',1,  @(x)isnumeric(x) && x>=1);
            
            % enter the refractive index of the first media
            p.addParamValue('n2',1.45,  @(x)isnumeric(x) && x>=1);
            
            % enter the angle of incidence (in degree)
            p.addParamValue('AoI',[],  @(x)isnumeric(x) && x>=0);
            
            
            p.parse(Grid_in,varargin{:});
            
            % Create the interface
            I.Grid =Grid_in;
            
            if p.Results.RoC == 0 || p.Results.RoC == Inf
                RoC = 1E99;
            else
                RoC = p.Results.RoC;
            end
            
            if isempty(p.Results.AoI) % if we arrive normal to the mirror
                I.surface =  -(RoC - sign(RoC)*sqrt(RoC^2 - I.Grid.D2_square));
            else
                I.surface =  -(I.Grid.D2_X.^2/(RoC*cos(p.Results.AoI*pi/180))+I.Grid.D2_Y.^2*cos(p.Results.AoI*pi/180)/RoC)/2;
            end
            
            I.n1 = p.Results.n1;
            I.n2 = p.Results.n2;
            
            I.T = p.Results.T;
            I.L = p.Results.L;
            
            if I.T + I.L >= 1
                I.t = 1i*sqrt(I.T-I.L);
                I.r = 0;
            else               
                I.t = 1i*sqrt(I.T);
                I.r = sqrt(1-(I.T + I.L));
            end
            
            %Mirror mask
            if isempty(p.Results.AoI)
                mask_index = I.Grid.D2_r < (p.Results.CA/2);
                I.mask = zeros(I.Grid.Num_point,I.Grid.Num_point,'double');
                I.mask(mask_index) = 1;
            else
                % define the radius in x and y of the aperture
                Rad_aperture_x = (p.Results.CA/2) * cos(p.Results.AoI * pi / 180);
                Rad_aperture_y = (p.Results.CA/2);
                
                mask_index = (I.Grid.D2_X.^2/Rad_aperture_x^2 + I.Grid.D2_Y.^2/Rad_aperture_y^2) < 1;
                I.mask = zeros(I.Grid.Num_point,I.Grid.Num_point,'double');
                I.mask(mask_index) = 1;
            end
            
        end
        
        
    end
    
    
end

