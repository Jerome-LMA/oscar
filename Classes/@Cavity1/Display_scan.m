function display_scan(obj)
% display_scan() Display the scan of the cavity. Can show the different
% mode excited in the cavity

if ~isa(obj, 'Cavity1')
    error('display_scan(C): The first and only argument must be an instance of the class Cavity1  ')
end

%  Initialize and hide the GUI as it is being constructed.
f = figure('Visible','off','Position',[100,300,800,700]);

hap1 = axes('Units','pixels','Position',[50,380,700,280]);datacursormode on;
hai1 = axes('Units','pixels','Position',[70,50,250,250]);
hai2 = axes('Units','pixels','Position',[470,50,250,250]);

% Check if we display the SB or not
Display_with_SB = true;

% Top plot: the cavity scan
semilogy(obj.cavity_scan_r(:,1),obj.cavity_scan_r(:,2),'Parent',hap1,'LineWidth',2)
axis(hap1,'tight')
title(hap1,'Circulating power vs resonance length')
dcm_obj = datacursormode(f);
set(dcm_obj,'UpdateFcn', @myupdatefcn )
xlabel(hap1,'Cavity detuning [m]')

% Bottom left: the input beam
imagesc(obj.laser_in.Grid.Axis,obj.laser_in.Grid.Axis,abs(obj.laser_in.Field).^2,'Parent',hai1)
title(hai1,'Cavity input beam')

% Bottom left: the circulating beam

[~,Length_scan_max_in] = max(obj.cavity_scan_r(:,2));
imagesc(obj.laser_in.Grid.Axis,obj.laser_in.Grid.Axis,abs(SumField(obj,obj.cavity_scan_r(Length_scan_max_in,1))),'Parent',hai2)
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
        
        imagesc(obj.laser_in.Grid.Axis,obj.laser_in.Grid.Axis,SumField(obj,pos(1)),'Parent',hai2)
        title(hai2,'Cavity circulating beam')
        %disp(['You clicked X:',num2str(pos(1)),', Y:',num2str(pos(2))]);
        txt = {['Position: ' num2str(pos(1))],['Power: ' num2str(pos(2))]};
    end

    function Field_reconstructed = SumField(Cin,length_reso) % Return the intensity of the field
        Grid_num_point = Cin.laser_in.Grid.Num_point;
        tmp = size(Cin.cavity_scan_all_field);
        num_iter = tmp(3);
        
        if Display_with_SB   % Calculate the the total field (carrier + SB) for display
            Field_reconstructed_car = complex(zeros(Grid_num_point,Grid_num_point,'double'));
            Field_reconstructed_SBu = Field_reconstructed_car; % do it even if no SB
            Field_reconstructed_SBl = Field_reconstructed_car;
            
            
            if Cin.laser_in.Nb_Pair_SB
                D_phi = (2*pi*Cin.laser_in.SB(1).Frequency_Offset/2.99792E8) * Cin.Length;
            end
            
            for ii=1:num_iter
                Field_reconstructed_car = Field_reconstructed_car + Cin.cavity_scan_all_field(:,:,ii) * exp(1i*Cin.laser_in.k_prop* length_reso*ii);
                if Cin.laser_in.Nb_Pair_SB
                Field_reconstructed_SBu = Field_reconstructed_SBu + Cin.cavity_scan_all_field(:,:,ii) * exp(1i*Cin.laser_in.k_prop* length_reso*ii) * exp(1i*D_phi*ii);
                Field_reconstructed_SBl = Field_reconstructed_SBl + Cin.cavity_scan_all_field(:,:,ii) * exp(1i*Cin.laser_in.k_prop* length_reso*ii) * exp(-1i*D_phi*ii);
                end
            end
            power_norm_SB =  calculate_power(Cin.laser_in,'SB')/2;
            Field_reconstructed = abs(Field_reconstructed_car).^2 + (abs(Field_reconstructed_SBu).^2)*power_norm_SB +...
                (abs(Field_reconstructed_SBl).^2)*power_norm_SB;
            
        else
            Field_reconstructed = complex(zeros(Grid_num_point,Grid_num_point,'double'));
            for ii=1:num_iter
                Field_reconstructed = Field_reconstructed + Cin.cavity_scan_all_field(:,:,ii) * exp(1i*Cin.laser_in.k_prop* length_reso*ii);
            end
            Field_reconstructed = abs(Field_reconstructed).^2;
        end
    end


end