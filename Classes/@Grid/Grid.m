classdef Grid < handle
    %Create all the 2D grids in the natural space and frequency space
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Always the first object to be defined in OSCAR
    % G1 = Grid();      Create a grid of one meter wide with 256 points
    % G1 = Grid( number of points , length );      Create a grid of length 'length' and with 'number of points' points  %
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        Num_point
        Length
        
        Step
        Step_sq
        Half_num_point
        Vector
        Axis
        Axis_FFT
        D2_X
        D2_Y
        D2_square
        D2_r
        D2_FFT_X
        D2_FFT_Y
    end
    
    methods
        
        function G1 = Grid(varargin)
            switch nargin
                case 0
                    G1.Num_point =256;
                    G1.Length = 1;
                case 1
                    G1.Num_point = 256;
                    G1.Length = varargin{2};
                case 2
                    G1.Num_point = varargin{1};
                    G1.Length = varargin{2};
                otherwise
                    disp('Grid(): invalid number of input arguments, grid not created')
                    return
            end
            
            G1.Step = G1.Length/G1.Num_point;
            G1.Step_sq = G1.Step^2;
            G1.Half_num_point = G1.Num_point/2;
            G1.Vector = 1:1:G1.Num_point;
            G1.Axis =  -G1.Length/2 + G1.Step/2 + (G1.Vector-1)*G1.Step;
            G1.Axis_FFT =  -1/(2*G1.Step) + (G1.Vector-1)/G1.Length;
            [G1.D2_X,G1.D2_Y] = meshgrid(G1.Axis);
            G1.D2_square = G1.D2_X.^2 + G1.D2_Y.^2;
            G1.D2_r   = sqrt(G1.D2_square);
            
            [G1.D2_FFT_X,G1.D2_FFT_Y] = meshgrid(G1.Axis_FFT);
            
        end
    end
    
end

