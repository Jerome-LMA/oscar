function varargout = Fit_E_Field(varargin)
%  Fit_E_Field() find the beam radius and the wavefront curvature for a
%  E_field object.
%  Fit the mode according to the variable Mode_name
%  !! Fitting is not an exact science check carefully the results
% Required the 'Optimisation Toolbox' to work

if (~exist('lsqcurvefit','file'))
    error('Fit_TEM00() required the Optimisation Toolbox to run')
end


switch nargin
    case {0}
        error('Fit_E_Field(): Not enough arguments, one object E_field must be given')
    case 1
        E = varargin{1};
        
        if isa(E,'E_Field')
            
            if Calculate_power(E)==0
                error('Fit_E_Field():: No power in the field')
            end
            
            % Create the vectors to fit
            V_power =abs(E.Field(:)).^2;
            V_grid(:,1) = E.Grid.D2_X(:);
            V_grid(:,2) = E.Grid.D2_Y(:);
            
            % Find which mode to fit
            [family m n] = Read_mode_name(E.Mode_name);
            
            if strcmp(family,'HG')
                func_mode = @(c,xdata) c(1)*(exp(-2*(xdata(:,1).^2+xdata(:,2).^2)/(c(2)^2)) )  .* ...
                    ( HermitePolynomial(m, sqrt(2)/c(2) * xdata(:,1)) .*...
                    HermitePolynomial(n, sqrt(2)/c(2) * xdata(:,2)) ).^2;
            else
               func_mode = @(c,xdata) c(1)*(exp(-2*(xdata(:,1).^2+xdata(:,2).^2)/(c(2)^2)) )  .* ...
                   ( LaguerrePolynomial(m, abs(n), 2*(xdata(:,1).^2+xdata(:,2).^2)/c(2)^2) ).^2 .*...
                   ( ((xdata(:,1).^2+xdata(:,2).^2)/c(2)^2) .^ (abs(n)/2) ).^2;               
            end
            
            % Option for the fit
            options = optimset('Display','off','MaxFunEvals',1E6,'TolFun',1E-8,'DiffMinChange',1E-8);
            
            % Try to find an initial guess, scan some possible radius
            
            first_try.radius = (1:1:10)*E.Grid.Length/50;
            
            for jj=1:length(first_try.radius)
                c0 = [max(max(V_power)) first_try.radius(jj)];
                [Map.fit_para,~,fit_residual,~,~] = lsqcurvefit(func_mode,c0,V_grid,V_power,[],[],options);
                first_try.residual(jj) = sum(abs(fit_residual), 'all');
                %Map.fit_para(2)
            end
            
             [~,first_try.ind] = min(first_try.residual);         
             
            c0 = [max(max(V_power)) first_try.radius(first_try.ind)];
            [Map.fit_para,~,~,~,~] = lsqcurvefit(func_mode,c0,V_grid,V_power,[],[],options);
            
            
            %figure(101)
            %imagesc(reshape(V_power,size(E.Field))); axis square
            
            %figure(102)
            % imagesc(reshape(func_mode(Map.fit_para,V_grid),size(E.Field))); axis square
            
            Beam_rad = abs(Map.fit_para(2));
       
            
            % Now find the RofC of the complex wavefront
            % So propagate the field along some distance and watch for the
            % change in beam radius
                      
            dist_propa = E.Grid.Length*100;
            beam_rad_vec = Beam_rad;
            dist_vec = 0;
            
            for jj=1:20;
                dist_vec = [dist_vec jj*dist_propa];
                
                E2 = Propagate_E(E,dist_vec(jj+1));
                V_power =abs(E2.Field(:)).^2;
                
                c0 = [max(max(V_power)) beam_rad_vec(jj)];
                [Map.fit_para,~,~,~,~] = lsqcurvefit(func_mode,c0,V_grid,V_power,[],[],options);
                
                beam_rad_vec = [beam_rad_vec abs(Map.fit_para(2))];              
            end
            
             func_beam_evo = @(c,xdata) c(1) * sqrt( 1 + ( (xdata - c(2)) / ((pi * c(1)^2 / E.Wavelength )) ).^2 );
            
             d0 = [Beam_rad 100];
             
             [beam_fit_para,~,~,~,~] = lsqcurvefit(func_beam_evo,d0,dist_vec,beam_rad_vec,[],[],options);
           
             waist_size = beam_fit_para(1);
             dist_wait = - beam_fit_para(2);  % minus sign to be compatible with OSCAR definition
   
             
             %            figure(103)
           % plot(dist_vec,beam_rad_vec,dist_vec,)
             
             
            Beam_RofC = dist_wait * (1 +   ( (pi * waist_size^2 / E.Wavelength ) / dist_wait    )^2  );

            switch nargout
                case 0
                    fprintf('Beam radius [m]: %g  \t \t Wavefront curvature [m]: %g  \n',Beam_rad,Beam_RofC)
                case 1
                    varargout{1} = Beam_rad;
                case 2
                    varargout{1} = Beam_rad;
                    varargout{2} = Beam_RofC;
                otherwise
                    error('Fit_E_Field(): Too many output argument')
            end
            
            
        else
            error('Fit_E_Field(): The first argument must be a E_Field')
        end
        
    otherwise
        error('Fit_TEM00():: invalid number of input arguments, no calculation is made')
end
