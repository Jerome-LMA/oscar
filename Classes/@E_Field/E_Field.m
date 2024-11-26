classdef E_Field
    % E_Field, this class represents an electric field defined on a grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % E = E_Field( Grid_name , 'w', 0.1)     Create a fundamental
    % Gaussian beam at the waist with a waist size of 10 cm
    % E = E_Field(Grid_name , 'w', 0.1,'R',-2000)    Create a fundamental mode with a defined beam radius of 10 cm and
    % wavefront radius of curvature of -2 km, so the beam propagates toward the waist
    % E = E_Field(Grid_name , 'w0', 0.1,'Z',-60)    Create a fundamental mode with a defined beam waist of 10 cm and
    % we are are 60m from the waist, so the beam propagates toward the waist
    % E = E_Field(Grid_name , 'q',-100 + 425*1i)    Create a fundamental
    % mode defined with the complex q parameter
    %  E = E_Field(Grid_name , 'w0', 0.1,'mode','LG m n') Create a mode Hermitte Gauss of order m,n
    %  E = E_Field(Grid_name , 'w0', 0.1,'mode','LG m n')  Create a mode Laguerre Gauss helicoidal d'ordre
    % p,l
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        Grid
        Field
        SB
        
        Refractive_index = 1;
        Wavelength = 1064E-9;
        Mode_name
        Optimal_mode_matching = false;
        
        ABCD_q
        
        Nb_Pair_SB = 0;
        k_prop
    end
    
    %     events
    %         Nan_present
    %     end
    
    methods
        function E = E_Field(Grid_in,varargin)
            p  = inputParser;
            p.FunctionName = 'E_field creator';
            
            % Check if the first argument is a grid
            p.addRequired('Grid_in', @(x)isa(x, 'Grid'));
            
            % Check if the wavelength of the laser beam is given
            p.addParameter('Wavelength',1064E-9, @(x)isnumeric(x) && x>0);
            
            % enter either the size of the beam radius and the complex
            % radius of curvature
            p.addParameter('w',[],@(x)isnumeric(x) && x>0);
            p.addParameter('R',1E99,  @(x)isnumeric(x));
            
            % can also enter the size of the waist and the distance from
            % the waist
            p.addParameter('w0',[],  @(x)isnumeric(x) && x>0);
            p.addParameter('z', 0,  @(x)isnumeric(x));
            
            % or enter the complex radius of curvature
            p.addParameter('q',[],  @(x)isnumeric(x));
            
            % Check if the power of the beam is given
            p.addParameter('Power',1, @(x)isnumeric(x) && x>0);
            
            % Check if we want the best mode matching, once the cavity is
            % defined
            p.addParameter('Optimal_MM',false,@(x)islogical(x));
            
            % Mode
            p.addParameter('Mode','HG 0 0', @(x)ischar(x));
            
            if nargin > 0
                p.parse(Grid_in,varargin{:});
                
                E.Wavelength = p.Results.Wavelength;
                E.k_prop = (2*pi/E.Wavelength);
                
                % check of what is entered
                
                if ( isempty(p.Results.q) && isempty(p.Results.w) && isempty(p.Results.w0) && ~p.Results.Optimal_MM )
                    error('E_Field(): at least the parameter w,w0 or q must be given to define the laser beam')
                end
                
                if  ~isempty(p.Results.w)
                    q_start = 1/(1/p.Results.R - 1i*(E.Wavelength)/(pi*p.Results.w^2));
                    beam_radius = p.Results.w;
                elseif ~isempty(p.Results.w0)
                    q_start = p.Results.z + 1i*pi*p.Results.w0^2/E.Wavelength;
                    beam_radius = sqrt( 1/(-imag(1/q_start)*pi/(E.Wavelength)) );
                elseif ~isempty(p.Results.q)
                    q_start =  p.Results.q;
                    beam_radius = sqrt( 1/(-imag(1/q_start)*pi/(E.Wavelength)) );
                elseif p.Results.Optimal_MM
                    q_start = 1/(- 1i*(E.Wavelength)/((pi*Grid_in.Length/10)^2));% put a dummy value there
                    beam_radius = sqrt( 1/(-imag(1/q_start)*pi/(E.Wavelength)) );
                    E.Optimal_mode_matching = true;
                end
                
                [family,m,n] = Read_mode_name(p.Results.Mode);
                
                E.Grid = Grid_in;
                
                if strcmp(family,'LG') % for the backward compatibility with previous version
                    family = 'LG_HELI';
                end
                
                if strcmp(family,'HG')
                    
                    E.Field = exp(-1i*E.k_prop*E.Grid.D2_square/(2*q_start));
                    
                    E.Field =  E.Field  .* HermitePolynomial(m, sqrt(2)/beam_radius * E.Grid.D2_X) .*...
                        HermitePolynomial(n, sqrt(2)/beam_radius * E.Grid.D2_Y);
                    
                elseif strcmp(family,'LG_HELI') % the one before V3.17
                    
                    E.Field = exp(-1i*E.k_prop*E.Grid.D2_square/(2*q_start));
                    
                    E.Field =  E.Field  .* (2* E.Grid.D2_square / beam_radius^2) .^ (abs(n)/2);
                    E.Field =  E.Field  .* LaguerrePolynomial(m, abs(n), 2* E.Grid.D2_square / beam_radius^2);
                    E.Field =  E.Field  .* exp(1i * n*atan2(E.Grid.D2_Y,E.Grid.D2_X));
                    
                elseif strcmp(family,'LG_SIN')
                    
                    E.Field = exp(-1i*E.k_prop*E.Grid.D2_square/(2*q_start));
                    
                    E.Field =  E.Field  .* (2* E.Grid.D2_square / beam_radius^2) .^ (abs(n)/2);
                    E.Field =  E.Field  .* LaguerrePolynomial(m, abs(n), 2* E.Grid.D2_square / beam_radius^2);
                    E.Field =  E.Field  .* exp(1i * n*atan2(E.Grid.D2_Y,E.Grid.D2_X)) + E.Field  .* exp(-1i * n*atan2(E.Grid.D2_Y,E.Grid.D2_X));
                    
                else
                    error('E_Field():the mode name must be HG or LG_HELI or LG_SIN')
                end
                
                E.Mode_name = p.Results.Mode;
                E = Normalise_E(E,'Power',p.Results.Power);
                
                E.ABCD_q = q_start;
                
            else
            end
        end
        
        function value = get.Field(obj)
            value = obj.Field;
        end
    end
    
end
