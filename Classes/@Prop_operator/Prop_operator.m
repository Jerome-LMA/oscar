classdef Prop_operator
    % Prop = Prop_operator(E1,1000,1.45)
    % Class to store the pre-computer matrix used for the propagation of
    % the beam

    properties
        n                           % the refractive index
        mat                       % the matrix of propagation
        dist                       % the distance of propagation
        Grid                       % the grid on which the propagation is defined

        Use_DI                     % Boolean to use the digital integration or not
        mat_DI                     % Propagation matrix for the digital integration

        Use_GPU                    % Boolean to use the calculation on the GPU, by default not
        ABCD_mat                   % propagation matrix from the ABCD method

    end

    methods

        function Prop = Prop_operator(E_in,dist,varargin)
            p  = inputParser;
            p.FunctionName = 'Prop_operator creator';

            % Check if the first argument is a E_Field object
            p.addRequired('E_in', @(x)isa(x, 'E_Field'));

            % Check if the first argument is a positice distance
            p.addRequired('dist', @(x)isnumeric(x) && x>0);

            % Then a third optionnal argument, to specify the refractive
            % index
            p.addParameter('n',[],@(x)isnumeric(x) && x>=1);

            % Another optional argument, use the digital integration or not for the propagation
            %
            p.addParameter('use_DI',false,@(x)isa(x,'logical'));


            p.parse(E_in,dist,varargin{:});

            % Now create the propagation operator
            % if n is given change the refractive index of the input laser
            % beam

            if  ~isempty(p.Results.n)
                E_in = Change_E_n(E_in,p.Results.n);
            end

            Prop.n = E_in.Refractive_index;
            Prop.dist = dist;
            Prop.Use_DI = p.Results.use_DI;
            Prop.Grid = E_in.Grid;
            Prop.Use_GPU = false;

            Prop.mat = exp(1i*(-E_in.k_prop*dist + ...
                pi*(E_in.Wavelength/Prop.n)*( E_in.Grid.D2_FFT_X.^2 + E_in.Grid.D2_FFT_Y.^2)*dist));

            Prop.ABCD_mat = [1 Prop.dist/Prop.n; 0 1]; % include the interface to the refractive index m


            % also add the propagation matrix for FFT-DI

            Prop.mat_DI = zeros(2*E_in.Grid.Num_point-1);

            X_tmp_1D = linspace(2*E_in.Grid.Axis(1),-2*E_in.Grid.Axis(1),2*E_in.Grid.Num_point-1);
            Y_tmp_1D = linspace(2*E_in.Grid.Axis(1),-2*E_in.Grid.Axis(1),2*E_in.Grid.Num_point-1);

            [X_tmp_2D, Y_tmp_2D] = meshgrid(X_tmp_1D,Y_tmp_1D);
            tmp_rad = sqrt(X_tmp_2D.^2 + Y_tmp_2D.^2 + Prop.dist^2);

            Prop.mat_DI = 1/(2*pi) * (exp(-1i * E_in.k_prop *tmp_rad) ./ tmp_rad) .* (Prop.dist ./  tmp_rad)  .* (1 ./  tmp_rad + 1i * E_in.k_prop );
            Prop.mat_DI = Prop.mat_DI * (E_in.Grid.Step)^2;

            Prop.mat_DI = fft2(Prop.mat_DI);

        end

    end

end

