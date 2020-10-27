function  varargout = Check_Stability(Cin,varargin)
%% Check_stability(C1) calculate RofC of the mirrors and the cavity stability
% Can do the RoC over a smaller area if a diameter is specified ( for example with Check_stability(C1,0.150) )
% First calculate the RofC of the mirrors, do a fit for the curvature of
% the surface.

p  = inputParser;

% Check if the first argument is a cavity
p.addRequired('Cin', @(x)isa(x, 'Cavity1'));

% Check if the results have to be displayed
p.addParameter('Display',true,@(x)islogical(x));

% Check if the results have to be displayed
p.addParameter('diam',0,@(x)isnumeric(x) && x>0);

p.parse(Cin,varargin{:})

Display = p.Results.Display;
diam_mir = p.Results.diam;

I1 = Cin.I_input;

% Check where the mirror is defined
if diam_mir
    I1_mask_index = find(I1.Grid.D2_r < diam_mir/2);
else
    I1_mask_index = find(I1.mask == 1);
end

% Define the variables for the fit

I1_grid2D(:,1) = I1.Grid.D2_X(I1_mask_index);
I1_grid2D(:,2) = I1.Grid.D2_Y(I1_mask_index);
I1_surf =  I1.surface(I1_mask_index);

% Define the quadratic fit, take into mirror offcenter, tilt and piston

func_curv = @(c,xdata) c(1)*((xdata(:,1)).^2 +(xdata(:,2)).^2) + c(2) +...
    c(3)*xdata(:,1) + c(4)*xdata(:,2);

% Find a rough initial guess, take one point at the center and then one
% slightly offset and guess the radius from there

sagitta_change = I1.surface(I1.Grid.Half_num_point,I1.Grid.Half_num_point) - I1.surface(I1.Grid.Half_num_point,I1.Grid.Half_num_point+10);
RofC_guess = (10*I1.Grid.Step)^2 / (2*sagitta_change);

if  isinf(RofC_guess)
    RofC_guess = 1E9;
end

% Initial guess
c0 = [-1/RofC_guess min(I1_surf) 0 0];

% tmp_g(:,1) = I1.Grid.D2_X(:);
% tmp_g(:,2) = I1.Grid.D2_Y(:);
%
% figure(2);
% imagesc(reshape(func_curv(c0,tmp_g),size( I1.Grid.D2_X)))


% Option for the fit
options = optimset('Display','off','MaxFunEvals',1E6,'TolFun',1E-12,'DiffMinChange',1E-12);
[Map.fit_para,~,residual,~,~] = lsqcurvefit(func_curv,c0,I1_grid2D,I1_surf,[],[],options);

I1_RofC = -1/(2*Map.fit_para(1));

if Display
    disp('----------------- For the input mirror -----------------')
    fprintf('RofC fitted (m): %g \n',I1_RofC)
    fprintf('Tilt horizontal (nrad): %g \n',Map.fit_para(3)*1E9)
    fprintf('Tilt vertical (nrad): %g \n',Map.fit_para(4)*1E9)
    fprintf('Flatness RMS (nm): %g,\n', std(residual)*1E9)
    disp('  ')
end

% Do the same thing for the second surface

I2 = Cin.I_end;

if diam_mir
    I2_mask_index = find(I2.Grid.D2_r < diam_mir/2);
else
    I2_mask_index = find(I2.mask == 1);
end

I2_grid2D(:,1) = I2.Grid.D2_X(I2_mask_index);
I2_grid2D(:,2) = I2.Grid.D2_Y(I2_mask_index);
I2_surf =  I2.surface(I2_mask_index);

func_curv = @(c,xdata) c(1)*((xdata(:,1)).^2 +(xdata(:,2)).^2) + c(2) +...
    c(3)*xdata(:,1) + c(4)*xdata(:,2);

sagitta_change = I2.surface(I2.Grid.Half_num_point,I2.Grid.Half_num_point) - I2.surface(I2.Grid.Half_num_point,I2.Grid.Half_num_point+10);
RofC_guess = (10*I2.Grid.Step)^2 / (2*sagitta_change);

if  isinf(RofC_guess)
    RofC_guess = 1E9;
end

c0 = [-1/RofC_guess min(I2_surf) 0 0];

options = optimset('Display','off','MaxFunEvals',1E6,'TolFun',1E-12,'DiffMinChange',1E-12);
[Map.fit_para,~,residual,~,~] = lsqcurvefit(func_curv,c0,I2_grid2D,I2_surf,[],[],options);

I2_RofC = -1/(2*Map.fit_para(1));

if Display
    disp('----------------- For the end mirror -----------------')
    fprintf('RofC fitted (m): %g \n',I2_RofC)
    fprintf('Tilt horizontal (nrad): %g \n',Map.fit_para(3)*1E9)
    fprintf('Tilt vertical (nrad): %g \n',Map.fit_para(4)*1E9)
    fprintf('Flatness RMS (nm): %g,\n', std(residual)*1E9)
    disp('  ')
end
% Now check the stability

g1 = 1 -  Cin.Length/I1_RofC;
g2 = 1 -  Cin.Length/I2_RofC;
g_factor_cavity = g1 * g2;

fprintf('The g-factor of the cavity is: %g \n',g_factor_cavity)


% Calculate various parameters if the cavity is stable
if (g_factor_cavity > 0) && (g_factor_cavity < 1)
    
    Lambda = Cin.Laser_in.Wavelength;
    Length_Cavity = Cin.Length;
    
    w0 = sqrt((Lambda * Length_Cavity / pi) * sqrt ((g1*g2*(1-g1*g2))/(g1+g2-2*g1*g2)^2));
    distITM_waist = (g2*(1 - g1)*Length_Cavity )/(g1+g2 - 2*g1*g2);
    
    W_onITM = sqrt( (Lambda * Length_Cavity / pi) * sqrt (g2 /(g1*(1-g1*g2))));
    W_onETM = sqrt( (Lambda * Length_Cavity / pi) * sqrt (g1 /(g2*(1-g1*g2))));
    
    % Mode separation in unit of the FSR
    Mode_sep = (1/pi)*acos(sqrt(g1*g2));
    
    cav_finesse = pi * ( (Cin.I_input.r * Cin.I_end.r)^.5 ) / (1 - Cin.I_input.r * Cin.I_end.r);
    if Display
        fprintf('The beam waist size in the cavity: %g \n',w0)
        fprintf('distance with the ITM: %g \n',distITM_waist)
        fprintf('Consecutive optical mode separation [1/FSR]: %g \n',Mode_sep)
        fprintf('\nBeam radius on ITM: %g \n',W_onITM)
        fprintf('Beam radius on ETM: %g \n',W_onETM)
        fprintf('\n Cavity finesse: %g \n', cav_finesse)
        fprintf(' Cavity gain: %g \n', abs(Cin.I_input.t)^2  /((1 - Cin.I_input.r * Cin.I_end.r))^2 )
    end
    %     % Calculate the parameters of the beam entering the cavity
    %     % Matrice ABCD of the propagation
    %     Cavity_q = 1/(-1/I1_RofC - 1i*Cin.Laser_in.Wavelength/(pi*W_onITM^2));
    %
    %     Mat_propa = [1 Cin.Length; 0 1]*[1 0;-2/I2_RofC 1]*[1 Cin.Length; 0 1];
    %
    %     q_propa = (Mat_propa(1,1)*Cavity_q + Mat_propa(1,2))/(Mat_propa(2,1)*Cavity_q + Mat_propa(2,2));
    %
    %     n =1;
    %     q_circ_inv = 1/(q_propa);
    %
    %     % Which gives us the beam parameter on the input mirror (toward the input mirror):
    %     RofC = 1/real(q_circ_inv);
    %     Beam_rad = sqrt( 1/(-imag(q_circ_inv)*pi/(1064E-9/n)));
    
    % Now make it pass the input optics
    Field_in = Transmit_Reflect_Optic(E_Field(Cin.Laser_in.Grid,'w',W_onITM,'R',I1_RofC),Cin.I_input,'HR');
    
    if isa(Cin.I_input, 'Interface')
        Field_in =  Change_E_n(Field_in,Cin.I_input.n1);
    end
    
    disp('Mode matched input beam parameters:')
    [wb, rb] = Fit_TEM00(Field_in);
    
    if Display
        fprintf('Beam radius [m]: %g  \t \t Wavefront curvature [m]: %g  \n',wb,-rb)
    end
    
    if nargout == 1
        varargout{1} = [wb,-rb];
    end
else
    if nargout == 1
        varargout{1} = {};
    end
    
end



end