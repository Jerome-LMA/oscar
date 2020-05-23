function Display_scan(Cin)
% Display_scan() Display the scan of the cavity. Can show the different
% mode excited in the cavity

if ~isa(Cin, 'Cavity1')
    error('Display_scan(C): The first and only argument must be an instance of the class Cavity1  ')
end

%  Initialize and hide the GUI as it is being constructed.
f = figure('Visible','off','Position',[100,300,800,700]);

hap1 = axes('Units','pixels','Position',[50,380,700,280]);datacursormode on;
hai1 = axes('Units','pixels','Position',[70,50,250,250]);
hai2 = axes('Units','pixels','Position',[470,50,250,250]);

% Check if we display the SB or not
Display_with_SB = true;

% Top plot: the cavity scan
semilogy(Cin.Cavity_scan_R(:,1),Cin.Cavity_scan_R(:,2),'Parent',hap1,'LineWidth',2)
axis(hap1,'tight')
title(hap1,'Circulating power vs resonance length')
dcm_obj = datacursormode(f);
set(dcm_obj,'UpdateFcn', @myupdatefcn )
xlabel(hap1,'Cavity detuning [m]')

% Bottom left: the input beam
imagesc(Cin.Laser_in.Grid.Axis,Cin.Laser_in.Grid.Axis,abs(Cin.Laser_in.Field),'Parent',hai1)
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
            
            D_phi = (2*pi*Cin.Laser_in.Frequency_Offset/2.99792E8) * Cin.Length;
            
            for ii=1:num_iter
                Field_reconstructed_car = Field_reconstructed_car + Cin.Cavity_scan_all_field(:,:,ii) * exp(1i*Cin.Laser_in.k_prop* length_reso*ii);
                Field_reconstructed_SBu = Field_reconstructed_SBu + Cin.Cavity_scan_all_field(:,:,ii) * exp(1i*Cin.Laser_in.k_prop* length_reso*ii) * exp(1i*D_phi*ii);
                Field_reconstructed_SBl = Field_reconstructed_SBl + Cin.Cavity_scan_all_field(:,:,ii) * exp(1i*Cin.Laser_in.k_prop* length_reso*ii) * exp(-1i*D_phi*ii);
            end
            power_norm_SB =  Calculate_power(Cin.Laser_in,'SB')/2;
            Field_reconstructed = abs(Field_reconstructed_car).^2 + (abs(Field_reconstructed_SBu).^2)*power_norm_SB +...
                (abs(Field_reconstructed_SBl).^2)*power_norm_SB;   
            
        else
            Field_reconstructed = complex(zeros(Grid_num_point,Grid_num_point,'double'));
            for ii=1:num_iter
                Field_reconstructed = Field_reconstructed + Cin.Cavity_scan_all_field(:,:,ii) * exp(1i*Cin.Laser_in.k_prop* length_reso*ii);
            end
            Field_reconstructed = abs(Field_reconstructed).^2;
        end
    end


end