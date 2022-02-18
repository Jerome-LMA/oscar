function  varargout = Display_Cavity_Modes(Cin,varargin)
% Display interactively the cavity modes

p = inputParser;
p.FunctionName = 'Display eigen modes';

% Check if the first argument is a cavity
p.addRequired('Cin', @(x)isa(x, 'CavityN'));

% Check if there is the number of modes to display
p.addParameter('N',20,@(x)isnumeric(x) );

% Check if the Airy peaks have to be displayed
p.addParameter('Airy',false,@(x)isnumeric(x) );

p.parse(Cin,varargin{:});

%p.Results

if isempty(Cin.Cavity_EM_mat)
    error('Display_cavity_modes(): cavity modes must first be calculated using the function Calculate_RT_mat()   ')
end

Cin = p.Results.Cin;
Nb_eigenvalue = p.Results.N;
Draw_Airy = p.Results.Airy;

% Calculate the Eigen vectors and eigen values
[V,D] = eigs(Cin.Cavity_EM_mat,Nb_eigenvalue);
Eigen_value = max((D));

% Tidy the eigen modes by diffraction loss
[Eigen_value,IX] = sort(Eigen_value,'descend');

% Check the number of point for the grid
[m,~] = size(Cin.Cavity_EM_mat);
Num_point = sqrt(m);

G_new = Grid(Num_point,Cin.I_array(1).Grid.Length);

Eigen_mode = complex(zeros(Num_point,Num_point,Nb_eigenvalue));
Reso_angle = zeros(Nb_eigenvalue,1);

for pp = 1:Nb_eigenvalue
    tmp_mode = reshape(V(:,IX(pp)),[Num_point Num_point]);
    Eigen_mode(:,:,pp) = tmp_mode;
    
    mode_brt = Normalise_E(E_Field(G_new,'w0',0.02),'Power',0);
    mode_brt.Field = Eigen_mode(:,:,pp);
    
    mode_brt = Resample_E(mode_brt,Cin.Laser_in.Grid);
    Field_Circ = mode_brt;
    
    for ii=1:Cin.Nb_mirror
        if ii ~= Cin.Nb_mirror % check we are not at the last iteration
            Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat_array(ii));
            Field_Circ = Reflect_Mirror(Field_Circ,Cin.I_array(ii+1),'Ref',1);
        else
            Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat_array(ii));
            Field_Circ = Reflect_Mirror(Field_Circ,Cin.I_array(1),'Ref',1);
        end
    end
    
    mode_art = Field_Circ;    
    Reso_angle(pp) = angle(Calculate_Overlap(mode_art,mode_brt));
end

% Bring the angle from 0 to pi, with the TEM00 at 0
Ind_00 = 1;

if (Reso_angle(Ind_00+1) - Reso_angle(Ind_00)) > pi
    Flip_sign = -1;
else
    Flip_sign = 1;
end

Vec_res_freq2 = mod(Reso_angle - Reso_angle(Ind_00),Flip_sign*2*pi);
Vec_res_loss =  abs(Eigen_value).^2;

if sum(Vec_res_loss > 1) >= 1
    error('Display_cavity_modes(): eigen value(s) superior to 1! likely cause: grid resolution too large to use digital integration.')
end


%  Initialize and hide the GUI as it is being constructed.
f = figure('Visible','off','Position',[100,100,1150,400]);

hap1 = axes('Units','pixels','Position',[100,50,600,320]);datacursormode on;
hai1 = axes('Units','pixels','Position',[750,50,320,320]);

scatter(Vec_res_freq2,Vec_res_loss,'filled','Parent',hap1)
ylim([min(Vec_res_loss)*0.95 1])
xlabel(hap1,'Cavity detuning [rad]')
ylabel(hap1,'1 - Round trip loss []')
title(hap1,'Eigen modes losses [m]')

dcm_obj = datacursormode(f);
set(dcm_obj,'UpdateFcn', @myupdatefcn )

%set(f,'Name',['Cavity scan of ' inputname(1)])
% Move the GUI to the center of the screen.
movegui(f,'center')

set(f,'Visible','on')

% Draw also the Airy peak

if Draw_Airy
    
    Phase_scan = linspace(-pi,pi,2000);
    
    r_loss = 1;
    for pp=1:Cin.Nb_mirror
        r_loss = r_loss * Cin.I_array(pp).r;
    end
        
    T = Cin.I_array(1).T;
    r = r_loss * abs(Eigen_value(1));
    
    figure;
    semilogy(Phase_scan,T./(abs(1 - r*exp(1i*(Phase_scan -Vec_res_freq2(1)) ) )).^2)
    %plot(Phase_scan+Vec_res_freq(1))
    
    hold all
    for pp=2:Nb_eigenvalue
        T = Cin.I_array(1).T;
        r =  r_loss * abs(Eigen_value(pp));       
        semilogy(Phase_scan,T./(abs(1 - r*exp(1i*(Phase_scan -Vec_res_freq2(pp)) ) )).^2)
    end
    hold off
    xlim([-pi pi])
    xlabel('Cavity detuning [rad]')
    ylabel('Cavity gain []')
    
end

if nargout == 1
    
    %  E_input = E_Field(G1,'w',0.02,'R',-2500);
    
    for pp = 1:Nb_eigenvalue
        tmp_E_field = E_Field(Cin.I_array(1).Grid,'w',0.02,'R',-2500); % define a
        tmp_E_field.Field = Eigen_mode(:,:,pp);
        E_out(pp) =  tmp_E_field;
    end
    varargout = {E_out};
end


    function txt = myupdatefcn(~, event_obj)
        pos = event_obj.Position;
        
        % find the index of the mode
        [~,ind] = min(abs(Vec_res_freq2 - pos(1)));
        
        imagesc(G_new.Axis,G_new.Axis,abs(Eigen_mode(:,:,ind)),'Parent',hai1)
        title(hai1,'Cavity circulating beam')
        %disp(['You clicked X:',num2str(pos(1)),', Y:',num2str(pos(2))]);
        txt = {['Detuning: ' num2str(pos(1))],['1 - Loss: ' num2str(pos(2))]};
    end

end

