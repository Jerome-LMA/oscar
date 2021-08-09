function varargout = Check_Stability(Cin,varargin)
%% Check_stability(Cin) calculate RofC of the mirrors and the cavity stability

% First calculate the RofC of the mirrors, do a fit for the curvature of
% the surface.

p  = inputParser;

% Check if the first argument is a cavity
p.addRequired('Cin', @(x)isa(x, 'CavityN'));

% Check if the results have to be displayed
p.addParameter('Display',true,@(x)islogical(x));

% Check if the results have to be displayed
p.addParameter('diam',0,@(x)isnumeric(x) && x>0);

p.parse(Cin,varargin{:})

Display = p.Results.Display;
diam_mir = p.Results.diam;



Nb_mirrors = length(Cin.I_array);
Roc_mirrrors = zeros(1,Nb_mirrors);
Ref_mir = zeros(1,Nb_mirrors);

for ii=1:Nb_mirrors
    IorM = Cin.I_array(ii);
    
    if isa(IorM, 'Interface')
        I1 = IorM;
    elseif isa(Cin.I_input, 'Mirror')
        I1 = Cin.I_array(ii).I_HR;
    else
        error('Check_stability(): serious problem, contact the administrator!')
    end
    
    % Check whether the mirror mirror is defined
    if diam_mir
        I1_mask_index = find(I1.Grid.D2_r < diam_mir/2);
    else
        I1_mask_index = find(I1.mask == 1);
    end
    
    Ref_mir(ii) = I1.r;
    
    % Define the variables for the fit
    
    I1_grid2D(:,1) = I1.Grid.D2_X(I1_mask_index);
    I1_grid2D(:,2) = I1.Grid.D2_Y(I1_mask_index);
    I1_surf =  I1.surface(I1_mask_index);
    
    % Define the quadratic fit, take into mirror offcenter, tilt and piston
    
    func_curv = @(c,xdata) c(1)*((xdata(:,1)-c(2)).^2 +(xdata(:,2)-c(3)).^2) + c(4) +...
        c(5)*xdata(:,1) + c(6)*xdata(:,2);
    
    % Find a rough initial guess, take one point at the center and then one
    % slightly offset and guess the radius from there
    
    sagitta_change = I1.surface(I1.Grid.Half_num_point,I1.Grid.Half_num_point) - I1.surface(I1.Grid.Half_num_point,I1.Grid.Half_num_point+10);
    RofC_guess = (10*I1.Grid.Step)^2 / (2*sagitta_change);
    
    if  isinf(RofC_guess)
        RofC_guess = 1E9;
    end
    
    % Initial guess
    c0 = [1/RofC_guess 0 0 min(I1_surf) 0 0];
    
    % Option for the fit
    options = optimset('Display','off','MaxFunEvals',1E6,'TolFun',1E-12,'DiffMinChange',1E-12);
    [Map.fit_para,~,residual,~,~] = lsqcurvefit(func_curv,c0,I1_grid2D,I1_surf,[],[],options);
    
    I1_RofC = -1/(2*Map.fit_para(1));
    if Display
        fprintf('----------------- For the surface %i  ----------------- \n',ii)
        
        fprintf('RofC fitted (m): %g \n',I1_RofC)
        fprintf('Center of the map, horizontal (mm): %g \n',Map.fit_para(2)*1E3)
        fprintf('Center of the map, vertical (mm): %g \n',Map.fit_para(3)*1E3)
        fprintf('Tilt horizontal (nrad): %g \n',Map.fit_para(5)*1E9)
        fprintf('Tilt vertical (nrad): %g \n',Map.fit_para(6)*1E9)
        fprintf('Flatness RMS (nm): %g,\n', std(residual)*1E9)
        disp('  ')
    end
    Roc_mirrrors(ii) = I1_RofC;
    
    clear I1_mask_index  I1_grid2D  I1_surf
    
end

% Now check the stability, create the ABCD matrix for one round trip
Mat_RT = [1 0;0 1];

if Cin.type == 'ring'
    for ii = 1:(Nb_mirrors - 1)
        Mat_RT = [1 0;-2/Roc_mirrrors(ii+1) 1]*[1 Cin.d_array(ii);0 1]*Mat_RT;
    end
    
    Mat_RT = [1 0;-2/Roc_mirrrors(1) 1]*[1 Cin.d_array(ii+1);0 1]*Mat_RT;
    
elseif Cin.type == 'folded'
    for ii = 1:(Nb_mirrors-1)
        Mat_RT = [1 0;-2/Roc_mirrrors(ii+1) 1]*[1 Cin.d_array(ii);0 1]*Mat_RT;
    end
    
    for pp= Cin.Nb_mirror-1:-1:1 % and do the round trip
        Mat_RT = [1 0;-2/Roc_mirrrors(pp) 1]*[1 Cin.d_array(pp);0 1]*Mat_RT;
        
    end
end


if abs( 0.5*(Mat_RT(1,1) + Mat_RT(2,2) ) ) > 1
    fprintf('!!   Cavity %s not stable   !!\n',inputname(1))
    if nargout == 1
        varargout{1} = {};
    end
else
    Lambda = Cin.Laser_in.Wavelength;
    
    % Calculate the beam parameter on the first mirror
    Roc_input = (2 * Mat_RT(1,2)) / (Mat_RT(2,2) - Mat_RT(1,1));
    beam_size = sqrt( (Lambda * abs(Mat_RT(1,2)) / pi) / ( (1 - ( 0.5*(Mat_RT(1,1)+Mat_RT(2,2)) )^2 ))^0.5  );
    if Display
        fprintf('  Beam radius on the first mirror: %g \n',beam_size)
    end
    
    q1 = 1/ ( 1/Roc_input - 1i*Lambda / (pi * beam_size^2) );
    if Cin.type == 'ring'
        
        for ii=2:Cin.Nb_mirror
            Mat_propa = [1 0;-2/Roc_mirrrors(ii) 1]*[1 Cin.d_array(ii-1);0 1];
            q2 = (Mat_propa(1,1)*q1 + Mat_propa(1,2))/(Mat_propa(2,1)*q1 + Mat_propa(2,2));
            q2_inv = 1/q2;
            Beam_rad = sqrt( 1/(-imag(q2_inv)*pi/(Lambda)));
            if Display
                fprintf('  Beam radius on mirror %i [m]: %g \n',ii,Beam_rad)
            end
            q1 = q2;
            
        end
        
    elseif Cin.type == 'folded'
        for ii = 1:(Nb_mirrors-1)
            Mat_propa = [1 0;-2/Roc_mirrrors(ii+1) 1]*[1 Cin.d_array(ii);0 1];
            q2 = (Mat_propa(1,1)*q1 + Mat_propa(1,2))/(Mat_propa(2,1)*q1 + Mat_propa(2,2));
            q2_inv = 1/q2;
            Beam_rad = sqrt( 1/(-imag(q2_inv)*pi/(Lambda)));
            if Display
                fprintf('  Beam radius on mirror %i [m]: %g \n',ii+1,Beam_rad)
            end
            q1 = q2;
        end
    end
    
    % Calculate the finesse and gain
    if Display
        if Cin.type == 'ring'
            
            cav_finesse = pi * ( (prod(Ref_mir))^.5 ) / (1 - prod(Ref_mir));
            fprintf('  Cavity finesse: %g \n', cav_finesse);
            fprintf('  Cavity gain: %g \n', (sqrt(1-Cin.I_array(1).r^2) / (1 - prod(Ref_mir)))^2 )
            
        elseif Cin.type == 'folded'
            RT_r = prod(Ref_mir) * prod(Ref_mir(2:end-1));
            
            cav_finesse = pi * ( RT_r^.5 ) / (1 - RT_r);
            fprintf('  Cavity finesse: %g \n', cav_finesse);
            fprintf('  Cavity gain: %g \n', (sqrt(1-Cin.I_array(1).r^2) / (1 - RT_r))^2 )
        end
        
    end
    
    % Now make it pass the input optics
    Field_in = Transmit_Reflect_Optic(E_Field(Cin.Laser_in.Grid,'w',beam_size,'R',-Roc_input),Cin.I_array(1),'HR');
    
    if isa(Cin.I_array(1), 'Interface')
        Field_in =  Change_E_n(Field_in,Cin.I_array(1).n1);
    end
    
    
    [wb,rb] = Fit_TEM00(Field_in);
    %E_plot(Field_in)
    
    if Display
        disp('  Mode matched input beam parameters:')
        fprintf('  Beam radius [m]: %g  \t \t Wavefront curvature [m]: %g  \n',wb,-rb)
    end
    
    if nargout == 1
        varargout{1} = [wb,-rb];
    end
end



end