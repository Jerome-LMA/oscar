function [Eout] = Propagate_E(varargin)
% Propagate: propagate the E_field over a distance
% E2 = Propagate(E1,dist): propagate the field E1 over the distance dist in
% meter
% E2 = Propagate(E1,Prop): propagate the field E1 using an instance of the Propagation
% Operator

switch nargin
    case {0,1}
        error('Propagate_E(): not enough arguments, at least an object E_field and a distance must be given')
    case 2
        
        E = varargin{1};
        dist = varargin{2};
        
        if ~isa(E, 'E_Field')
            error('Propagate_E(): the first argument must be an E_Field' )
        end
        
        Eout = E;
        
        if isscalar(dist) && isreal(dist) % the user input a distance
            
            % Create the propagation matrix
            Propa_mat =  exp(1i*(-E.k_prop*dist  + ...
                pi*(E.Wavelength/E.Refractive_index)*( E.Grid.D2_FFT_X.^2 +E.Grid.D2_FFT_Y.^2)*dist));
            
            Length_propa = dist;
            Use_DI = false;
            
            %imagesc(angle(Propa_mat)); axis square
       
            
        elseif isa(dist, 'Prop_operator')
            Propa_mat = dist.mat;
            Length_propa = dist.dist;
            Use_DI = dist.Use_DI;
            
            if E.Refractive_index ~= dist.n
                error('Propagate_E(): the refractive index from the E_Field does not match the one from the propagation operator')
            end
            
        else
            error('Propagate_E(): wrong second argument' )
        end
        
        if isempty(E.Field_SBl) % if no sidebands are present
            
            if ~Use_DI
                Wave_fft = fftshift(fft2(E.Field));
                Wave_prop = Propa_mat .* Wave_fft;
           
                Eout.Field = ifft2(ifftshift(Wave_prop));
                
            else                            
                Mat_U = zeros(2*E.Grid.Num_point - 1);
                Mat_U(1:E.Grid.Num_point,1:E.Grid.Num_point) = E.Field;
                
                S = ifft2(fft2(Mat_U) .* dist.mat_DI);
                Eout.Field = S(E.Grid.Num_point:end,E.Grid.Num_point:end);
                
            end
            
        else % if we have also sidebands
            
            D_phi = (2*pi*E.Frequency_Offset/2.99792E8) * Length_propa;
            
            if ~Use_DI
                
                % for the carrier
                Wave_fft = fftshift(fft2(E.Field));
                Wave_prop = Propa_mat .* Wave_fft;
                Eout.Field = ifft2(ifftshift(Wave_prop));
                
                % for the lower sideband
                Wave_fft = fftshift(fft2(E.Field_SBl));
                Wave_prop = Propa_mat .* Wave_fft;
                Eout.Field_SBl = ifft2(ifftshift(Wave_prop));
                Eout.Field_SBl = Eout.Field_SBl * exp(1i*D_phi);
                
                % for the upper sideband
                Wave_fft = fftshift(fft2(E.Field_SBu));
                Wave_prop = Propa_mat .* Wave_fft;
                Eout.Field_SBu = ifft2(ifftshift(Wave_prop));
                Eout.Field_SBu = Eout.Field_SBu * exp(-1i*D_phi);
                
            else
                
                % for the carrier
                Mat_U = zeros(2*E.Grid.Num_point - 1);
                Mat_U(1:E.Grid.Num_point,1:E.Grid.Num_point) = E.Field;
                S = ifft2(fft2(Mat_U) .* dist.mat_DI);
                Eout.Field = S(E.Grid.Num_point:end,E.Grid.Num_point:end);
                
                % for the lower sideband
                Mat_U = zeros(2*E.Grid.Num_point - 1);
                Mat_U(1:E.Grid.Num_point,1:E.Grid.Num_point) = E.Field_SBl;
                S = ifft2(fft2(Mat_U) .* dist.mat_DI);
                Eout.Field_SBl = S(E.Grid.Num_point:end,E.Grid.Num_point:end)* exp(1i*D_phi);
                
                % for the upper sideband
                Mat_U = zeros(2*E.Grid.Num_point - 1);
                Mat_U(1:E.Grid.Num_point,1:E.Grid.Num_point) = E.Field_SBu;
                S = ifft2(fft2(Mat_U) .* dist.mat_DI);
                Eout.Field_SBu = S(E.Grid.Num_point:end,E.Grid.Num_point:end)* exp(-1i*D_phi);
                                
            end
            
        end
        
        
    otherwise
        error('Propagate_E():Invalid number of input arguments')
        
end
end

