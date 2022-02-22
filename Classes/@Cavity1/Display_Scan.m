function Display_Scan(Cin,varargin)
% Display_scan() Display the scan of the cavity. Can show the different
% mode excited in the cavity

p  = inputParser;

% Check if the first argument is an E_Field
p.addRequired('Cin', @(x)isa(x, 'Cavity1'));

% Horizontal axis of the scan
p.addParameter('scan','length', @(x)strcmpi(x,'length') | ...
    strcmpi(x,'frequency') | strcmpi(x,'RT phase') );

% Check if the circulating is displayed with the SB or not
p.addParameter('Display_with_SB',true,@(x)islogical(x));

p.parse(Cin,varargin{:})

%  Initialize and hide the GUI as it is being constructed.
f = figure('Visible','off','Position',[100,300,800,700]);

hap1 = axes('Units','pixels','Position',[50,380,700,280]);datacursormode on;
hai1 = axes('Units','pixels','Position',[70,50,250,250]);
hai2 = axes('Units','pixels','Position',[470,50,250,250]);

% Check if we display the SB or not
Display_with_SB = p.Results.Display_with_SB;

% Build the horizontal axis
if strcmp(p.Results.scan,'length') % display the horizontal axis as cavity length change
    Coef_scale = 2;
    Hor_axis = Cin.Cavity_scan_R(:,1) / Coef_scale;
    Hor_title = 'Cavity length detuning [m]';
elseif strcmp(p.Results.scan,'frequency') % display the horizontal axis as laser frequency change
    Coef_scale = Cin.Laser_in.Wavelength / (299792458/(2*Cin.Length));
    %Hor_axis = (Cin.Cavity_scan_R(:,1)/(Cin.Laser_in.Wavelength))*(299792458/(2*Cin.Length));
    Hor_axis = Cin.Cavity_scan_R(:,1) / Coef_scale;
    Hor_title = 'Cavity frequency shift [Hz]';
else  % display the horizontal axis as cavity round trip phase change
    %Hor_axis = (Cin.Cavity_scan_R(:,1)/(Cin.Laser_in.Wavelength))*2;
    Coef_scale = Cin.Laser_in.Wavelength /2;
    Hor_axis = Cin.Cavity_scan_R(:,1) / Coef_scale;
    Hor_title = 'Cavity round trip phase shift [pi]';
end

% Top plot: the cavity scan
semilogy(Hor_axis,Cin.Cavity_scan_R(:,2),'Parent',hap1,'LineWidth',2)
axis(hap1,'tight')
title(hap1,'Circulating power vs resonance length')
dcm_obj = datacursormode(f);
set(dcm_obj,'UpdateFcn', @myupdatefcn )
xlabel(hap1,Hor_title)


% Bottom left: the input beam
imagesc(Cin.Laser_in.Grid.Axis,Cin.Laser_in.Grid.Axis,abs(Cin.Laser_in.Field).^2,'Parent',hai1)
title(hai1,'Cavity input beam')

% Bottom left: the circulating beam

[~,Length_scan_max_in] = max(Cin.Cavity_scan_R(:,2));
imagesc(Cin.Laser_in.Grid.Axis,Cin.Laser_in.Grid.Axis,abs(SumField(Cin,Cin.Cavity_scan_R(Length_scan_max_in,1))),'Parent',hai2)
title(hai2,'Cavity circulating beam')


% Change units to normalized so components resize
% automatically.
set([f,hap1,hai1,hai2],...
    'Units','normalized');

set(f,'Name',['Cavity scan of ' inputname(1)])
% Move the GUI to the center of the screen.
movegui(f,'center')
% Make the GUI visible.
set(f,'Visible','on')


%imagesc(peaks(40),'Parent',yourGUIaxehandle)

    function txt = myupdatefcn(~, event_obj)
        pos = event_obj.Position;
        
        imagesc(Cin.Laser_in.Grid.Axis,Cin.Laser_in.Grid.Axis,SumField(Cin,pos(1)),'Parent',hai2)
        title(hai2,'Cavity circulating beam')
        %disp(['You clicked X:',num2str(pos(1)),', Y:',num2str(pos(2))]);
        txt = {['Position: ' num2str(pos(1))],['Power: ' num2str(pos(2))]};
    end

    function Field_reconstructed = SumField(Cin,length_reso) % Return the intensity of the field
        Grid_num_point = Cin.Laser_in.Grid.Num_point;
        tmp = size(Cin.Cavity_scan_all_field);
        num_iter = tmp(3);
        
        if Display_with_SB   % Calculate the the total field (carrier + SB) for display
            Field_reconstructed_car = complex(zeros(Grid_num_point,Grid_num_point,'double'));
            Field_reconstructed_SBu = Field_reconstructed_car; % do it even if no SB
            Field_reconstructed_SBl = Field_reconstructed_car;
            
            
            if Cin.Laser_in.Nb_Pair_SB
                D_phi = (2*pi*Cin.Laser_in.SB(1).Frequency_Offset/2.99792E8) * Cin.Length;
            end
            
            for ii=1:num_iter
                Field_reconstructed_car = Field_reconstructed_car + Cin.Cavity_scan_all_field(:,:,ii) * exp(-1i*Cin.Laser_in.k_prop* length_reso*ii*Coef_scale);
                if Cin.Laser_in.Nb_Pair_SB
                    Field_reconstructed_SBu = Field_reconstructed_SBu + Cin.Cavity_scan_all_field(:,:,ii) * exp(-1i*Cin.Laser_in.k_prop* length_reso*ii*Coef_scale) * exp(1i*D_phi*ii);
                    Field_reconstructed_SBl = Field_reconstructed_SBl + Cin.Cavity_scan_all_field(:,:,ii) * exp(-1i*Cin.Laser_in.k_prop* length_reso*ii*Coef_scale) * exp(-1i*D_phi*ii);
                end
            end
            power_norm_SB =  Calculate_Power(Cin.Laser_in,'SB')/2;
            Field_reconstructed = abs(Field_reconstructed_car).^2 + (abs(Field_reconstructed_SBu).^2)*power_norm_SB +...
                (abs(Field_reconstructed_SBl).^2)*power_norm_SB;
            
        else
            Field_reconstructed = complex(zeros(Grid_num_point,Grid_num_point,'double'));
            for ii=1:num_iter
                Field_reconstructed = Field_reconstructed + Cin.Cavity_scan_all_field(:,:,ii) * exp(-1i*Cin.Laser_in.k_prop* length_reso*ii*Coef_scale);
            end
            Field_reconstructed = abs(Field_reconstructed).^2;
        end
    end


end