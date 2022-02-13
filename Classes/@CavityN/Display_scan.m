function Display_Scan(Cin)
% display_scan() Display the scan of the cavity. Can show the different
% mode excited in the cavity

if ~isa(Cin, 'CavityN')
    error('display_scan(): The first and only argument must be an instance of the class Cavity1  ')
end

%  Initialize and hide the GUI as it is being constructed.
f = figure('Visible','off','Position',[100,300,800,700]);

hap1 = axes('Units','pixels','Position',[50,380,700,280]);datacursormode on;
hai1 = axes('Units','pixels','Position',[70,50,250,250]);
hai2 = axes('Units','pixels','Position',[470,50,250,250]);


% Top plot: the cavity scan
semilogy(Cin.cavity_scan_r(:,1),Cin.cavity_scan_r(:,2),'Parent',hap1,'LineWidth',2)
axis(hap1,'tight')
title(hap1,'Circulating power vs resonance length')
dcm_obj = datacursormode(f);
set(dcm_obj,'UpdateFcn', @myupdatefcn )
xlabel(hap1,'Cavity detuning [m]')

% Bottom left: the input beam
imagesc(Cin.laser_in.Grid.Axis,Cin.laser_in.Grid.Axis,abs(Cin.laser_in.Field).^2,'Parent',hai1)
title(hai1,'Cavity input beam')

% Bottom left: the circulating beam

[~,Length_scan_max_in] = max(Cin.cavity_scan_r(:,2));
imagesc(Cin.laser_in.Grid.Axis,Cin.laser_in.Grid.Axis,abs(SumField(Cin,Cin.cavity_scan_r(Length_scan_max_in,1))).^2,'Parent',hai2)
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
  
  imagesc(Cin.laser_in.Grid.Axis,Cin.laser_in.Grid.Axis,abs(SumField(Cin,pos(1))).^2,'Parent',hai2)
  title(hai2,'Cavity circulating beam')
  %disp(['You clicked X:',num2str(pos(1)),', Y:',num2str(pos(2))]);
  txt = {['Position: ' num2str(pos(1))],['Power: ' num2str(pos(2))]};
end

    function Field_reconstructed = SumField(Cin,length_reso)
        Grid_num_point = Cin.laser_in.Grid.Num_point;
        tmp = size(Cin.cavity_scan_all_field);
        num_iter = tmp(3);
        
        Field_reconstructed = complex(zeros(Grid_num_point,Grid_num_point,'double'));
        for ii=1:num_iter
            Field_reconstructed = Field_reconstructed + Cin.cavity_scan_all_field(:,:,ii) * exp(1i*Cin.laser_in.k_prop* length_reso*ii);
        end       
    end


end