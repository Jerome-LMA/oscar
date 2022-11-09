function [Eout,Gout] = Focus_Beam_With_Telescope(Ein,array_L_D,varargin)
% Focus_beam_with_lens() use this function for lens with short focal length
%  This function re-adapt the grid size since the focusing may lead to a
%  small beam size
%  array_L_D is a vector representing the focal length of the lens and the
%  distance of propagation, several lens and distance can be simulated. For
%  example if array_L_D = [100 5 -10 4] means that the beam first cross a
%  lens of focal length 100m then propagate 5m, then meet another lens of
%  focal length -10m and finally propagate 4m.

% Possibility to know add the aperture of the telescope length (aperture = diameter of the optics)

p = inputParser;
p.FunctionName = 'Simulate a telescope';

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check if the secondt argument is a vector
p.addRequired('array_L_D', @(x)isvector(x));

% Check the optional third argument (the magnification)
p.addParameter('magnification',[], @(x)isvector(x));

% Use the Digital Integration for the propagation
p.addParameter('Use_DI',false, @(x)islogical(x));

% Possibility to set the diameter of the optics
p.addParameter('aperture',[], @(x)isvector(x));

% Possibility to add a surface map for the optic (map of wavefront transmission, so twice the surface, - tilt/curvature removed)
p.addParameter('map',[], @(x)isvector(x));


p.parse(Ein,array_L_D,varargin{:})

% Nomber of iteration for the focusing and propagation:
num_iter = length(array_L_D) / 2;
Focal_length = array_L_D(1+(0:num_iter-1)*2);
Distance = array_L_D(2+(0:num_iter-1)*2);
Aperture = p.Results.aperture;
Map_list = p.Results.map;

if ~isempty(Aperture)
    if length(Aperture) ~= num_iter
        error('Focus_beam_with_telescope(): incorrect number of aperture')
    end
end

if ~isempty(Map_list)
    if length(Map_list) ~= num_iter
        error('Focus_beam_with_telescope(): incorrect number of maps')
    end
end

if isempty(p.Results.magnification)
    %------------------------------------------------------------------------------------
    % If no magnification is set, use the optimal one
    
    % first derive the parameters of the input beam and remove the wavefront
    [~,R_in] = Fit_TEM00(Ein);
    
    if isinf(R_in)
        R_in = 1E99;
    end
    
    WF_change = exp(1i * Ein.k_prop *  Ein.Grid.D2_square * (1/(2*R_in)) );
    Ein = Ein .* WF_change; % will do for SB at the same time
    
    for pp=1:num_iter
        % Add the aperture here
        if ~isempty(Aperture)
            Ein = Transmit_Aperture(Ein,Aperture(pp));
        end
        % Add the map there
        if ~isempty(Map_list)
            % First rescale the map on the calculation grid
            I_temp = Resize_Interface(Map_list(pp),Ein.Grid);
            Ein = Ein .* (exp(1i * Ein.k_prop * I_temp.surface) .* I_temp.mask);
        end
        
        
        %figure(1); E_plot(Ein); figure(2); I_plot(I_temp); pause()
        
        % Change the sign of the lens to be consitent with the following
        f_lens = - Focal_length(pp);
        
        if (f_lens + R_in) ~= 0
            f_new = -(f_lens * R_in)/(f_lens + R_in);
        else
            f_new = 1E99; % So the lens will cancel the incident wavefront of curvature
        end
        
        % Corrected distance to propagate:
        dist_prop = - Distance(pp) * f_new / (Distance(pp) - f_new);
        
        if isinf(dist_prop)
            dist_prop = 1E99;
        end
        
        %dist_prop = - Distance(pp) * 1 / (-Distance(pp)*(f_lens + R_in)/(f_lens * R_in) -1);
        
        if p.Results.Use_DI
            Propagation_OP = Prop_operator(Ein,dist_prop,'use_DI',true);
            E_prop = Propagate_E(Ein,Propagation_OP);
        else
            E_prop = Propagate_E(Ein,dist_prop);
        end
        
        Scaling_factor = (f_new - Distance(pp) ) / f_new;
        
        Gout = Grid(Ein.Grid.Num_point,Ein.Grid.Length*abs(Scaling_factor));
        
        Eout = Ein;
        Eout.Grid = Gout;
        Eout.Field = E_prop.Field * (1/Scaling_factor);
        
        %         if ~isempty(Eout.Field_SBl)
        %             Eout.Field_SBl = Eout.Field_SBl / Scaling_factor;
        %             Eout.Field_SBu = Eout.Field_SBu / Scaling_factor;
        %         end
        
        %E_Plot(Eout)
        
        % Check if no power fall outside the grid
        Check_Grid_size(E_prop,0.10)
        New_wavefront = -(Distance(pp)  - f_new);
        
        % Remove the wavefront curvature from the beam
        %E_Plot(Eout)
        [~,R_in2] = Fit_TEM00(Eout);
        WF_change = exp(1i * Ein.k_prop *  Eout.Grid.D2_square * (1/(2*R_in2)) );
        Eout = Eout .* WF_change;
        
        %         if ~isempty(Eout.Field_SBl)
        %             Eout.Field_SBl = Eout.Field_SBl .* WF_change;
        %             Eout.Field_SBu = Eout.Field_SBu .* WF_change;
        %         end
        
        % Calculate the new wavefront
        R_in = 1/(1/R_in2 - 1/New_wavefront);
        Ein = Eout;
        
    end
    
    % At the end, add the final wavefront
    WF_final = exp(1i * Eout.k_prop *  Eout.Grid.D2_square * (-1/(2*R_in)) );
    Eout = Eout.* WF_final;
    
    %     if ~isempty(Eout.Field_SBl)
    %         Eout.Field_SBl = Eout.Field_SBl .* WF_final;
    %         Eout.Field_SBu = Eout.Field_SBu .* WF_final;
    %     end
    
else
    % If the user set the magnification for the simulation
    Vec_mag = p.Results.magnification;
    
    % If magnification is one, no change of the Grid, for the procedure it
    % should be -1. So here the line to interprate correctly:
    Vec_mag(Vec_mag == 1) = -1;
    % Check we have the right number of magnification (must be equal to the number of lenses)
    if num_iter ~= length(Vec_mag)
        error('Focus_beam_with_telescope(): number of lenses and magnification do not match')
    end
    
    Focal_length = array_L_D(1+(0:num_iter-1)*2);
    Distance = array_L_D(2+(0:num_iter-1)*2);
    
    %  Roc_from_previous_iter = 1E99;
    
    for pp=1:num_iter
        
        % Add the aperture here
        if ~isempty(Aperture)
            Ein = Transmit_Aperture(Ein,Aperture(pp));
        end
        % Add the map there
        if ~isempty(Map_list)
            % First rescale the map on the calculation grid
            I_temp = Resize_Interface(Map_list(pp),Ein.Grid);
            Ein = Ein .* (exp(1i * Ein.k_prop * I_temp.surface) .* I_temp.mask);
        end
        %figure(1); E_plot(Ein); figure(2); I_plot(I_temp); pause()
        %figure(1); imagesc(I_temp.mask);axis square; figure(2); I_plot(I_temp); pause()
        
        New_RoC_lens = Focal_length(pp);
        New_RoC_mag = Distance(pp) / (Vec_mag(pp) - 1);
        
        New_RoC_total = 1/ ( 1/New_RoC_lens + 1/New_RoC_mag);
        WF_change = exp(1i * Ein.k_prop *  Ein.Grid.D2_square * (1/(2*New_RoC_total)) );
        Ein = Ein .* WF_change;
        
        dist_prop =   Distance(pp) / Vec_mag(pp);
        
        if p.Results.Use_DI
            Propagation_OP = Prop_operator(Ein,dist_prop,'use_DI',true);
            E_prop = Propagate_E(Ein,Propagation_OP);
        else
            E_prop = Propagate_E(Ein,dist_prop);
        end
        
        Eout = Ein;
        
        if Vec_mag(pp) == -1 % same magnification, so do not change the grid
            Eout.Grid = Ein.Grid;
            Gout = Eout.Grid;
            disp('test1')
        else
            Gout = Grid(Ein.Grid.Num_point,Ein.Grid.Length*abs(Vec_mag(pp)));
            Eout.Grid = Gout;
        end
        
        Eout.Field = E_prop.Field / Vec_mag(pp);
        
        %New_RoC = -New_RoC_mag^2 / (New_RoC_mag + dist);
        New_RoC = (Vec_mag(pp) / (1 - Vec_mag(pp)) * Distance(pp));
        %Roc_from_previous_iter = New_RoC;
        Eout.Field = Eout.Field .* exp(1i * Eout.k_prop *  Eout.Grid.D2_square * (1/(2*New_RoC) ));
        Ein = Eout;
        
    end
    
    %    Eout.Field = Eout.Field .* exp(1i * Eout.k_prop *  Eout.Grid.D2_square * (1/(2*New_RoC) ));
    
    
end

end

