function Eout = Propagate_E(Ein,varargin)
% Propagate: propagate the E_field over a distance
% E2 = Propagate(E1,dist): propagate the field E1 over the distance dist in
% meter, the distance could be negative if you what you are doing.
% E2 = Propagate(E1,Prop): propagate the field E1 using an instance of the Propagation
% Operator

p  = inputParser;

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check if the second argument is a distance or a propagation operator
p.addRequired('dist',@(x) (isscalar(x) && isreal(x) || isa(x,'Prop_operator'))  );

p.parse(Ein,varargin{:})

dist = p.Results.dist;

E = Ein;
Eout = Ein;

if isscalar(dist) && isreal(dist) % the user input a distance
    
    % Create the propagation matrix
    Propa_mat =  exp(1i*(-E.k_prop*dist  + ...
        pi*(E.Wavelength/E.Refractive_index)*( E.Grid.D2_FFT_X.^2 +E.Grid.D2_FFT_Y.^2)*dist));
    
    Length_propa = dist;
    Use_DI = false;
    
    %imagesc(angle(Propa_mat)); axis square
     
else 
    Propa_mat = dist.mat;
    Length_propa = dist.dist;
    Use_DI = dist.Use_DI;
    
    if E.Refractive_index ~= dist.n
        error('Propagate_E(): the refractive index from the E_Field does not match the one from the propagation operator')
    end
    
end


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

if E.Nb_Pair_SB % if we have also sidebands
    for ii=1:E.Nb_Pair_SB
        D_phi = (2*pi*E.SB(ii).Frequency_Offset/2.99792E8) * Length_propa;
        
        if ~Use_DI
            % for the lower sideband
            Wave_fft = fftshift(fft2(E.SB(ii).Field_lower));
            Wave_prop = Propa_mat .* Wave_fft;
            Eout.SB(ii).Field_lower = ifft2(ifftshift(Wave_prop));
            Eout.SB(ii).Field_lower = Eout.SB(ii).Field_lower * exp(1i*D_phi);
            
            % for the upper sideband
            Wave_fft = fftshift(fft2(E.SB(ii).Field_upper));
            Wave_prop = Propa_mat .* Wave_fft;
            Eout.SB(ii).Field_upper = ifft2(ifftshift(Wave_prop));
            Eout.SB(ii).Field_upper = Eout.SB(ii).Field_upper * exp(-1i*D_phi);
        else
            % for the lower sideband
            Mat_U = zeros(2*E.Grid.Num_point - 1);
            Mat_U(1:E.Grid.Num_point,1:E.Grid.Num_point) = E.SB(ii).Field_lower;
            S = ifft2(fft2(Mat_U) .* dist.mat_DI);
            Eout.SB(ii).Field_lower = S(E.Grid.Num_point:end,E.Grid.Num_point:end) * exp(1i*D_phi);
            
            
            % for the upper sideband
            Mat_U = zeros(2*E.Grid.Num_point - 1);
            Mat_U(1:E.Grid.Num_point,1:E.Grid.Num_point) = E.SB(ii).Field_upper;
            S = ifft2(fft2(Mat_U) .* dist.mat_DI);
            Eout.SB(ii).Field_upper = S(E.Grid.Num_point:end,E.Grid.Num_point:end)* exp(-1i*D_phi);
            
        end
    end
end


end

