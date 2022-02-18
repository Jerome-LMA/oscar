function Iout = Add_Astigmatism(varargin)
% Iout = Add_astigmatism(Iin,value,diameter)
% Add astigmatism on a surface Iin. The astigmatism'value' is defined as the
% amplitude of the Zernike polynomial 2,2 on a unity circle of diameter
% 'diameter'


switch nargin
    case {0,1,2}
        error('Add_astigmatism(): not enough arguments, at least an object E_field (or surface) and an angle must be given')
        
    case 3
        
        Iin = varargin{1};
        zer.amp =  varargin{2};
        zer.diam =  varargin{3};
        
        Iout = Iin;
        
        if  ~isa(Iin, 'Interface')
            error('Add_astigmatism(): the first argument must be an object Interface')
        end
        
        if  ~isscalar(zer.amp)
            error('Add_astigmatism(): the second argument must be a number: the height of the Zernike polynomial')
        end
        
        if  ~isscalar(zer.diam)
            error('Add_astigmatism(): the third argument must be a number: the diameter of the Zernike circle')
        end
         
        zer.plane2D = Iin.Grid.D2_X + 1i*Iin.Grid.D2_Y;
        zer.radius2D = abs(zer.plane2D/(zer.diam/2));
        zer.angle2D = angle(zer.plane2D);
      % imagesc(Iin.Grid.Axis,Iin.Grid.Axis,zer.radius2D);axis square 
        zer.map2D = zer.amp*(zer.radius2D.^2) .* cos(2*zer.angle2D);
    %    imagesc(Iin.Grid.Axis,Iin.Grid.Axis,zer.map2D);axis square
        
        Iout.surface =  Iin.surface + zer.map2D;
        
        
    otherwise
        error('Add_astigmatism(): Invalid number of input arguments, no tilt has been added')
        
end

