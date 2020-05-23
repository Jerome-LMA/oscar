function Display_results(Cin)
% Display_results(Cin) display the results of the cavity calculations.
% The function Calculate_fields() must have been run first

if isempty(Cin.Field_ref)
    error('Display_results(): Before displaying the results, the function Calculate_fields() must be run     ')
end

if ~isempty(Cin.Laser_in.Field_SBl)
    disp('---------- For the carrier ---------------')
end

% Calculate the total transmitted power
tmp_power = 0;

if Cin.type == 'ring'
    for pp=2:Cin.Nb_mirror
        tmp_power = tmp_power + Calculate_power(Cin.Field_trans(pp));
    end
elseif Cin.type == 'folded'
    tmp_power = Calculate_power(Cin.Field_trans(end));
end

fprintf(' Power in the input beam %g [W] \n',Calculate_power(Cin.Laser_in))
fprintf(' Circulating power %g [W] \n',Calculate_power(Cin.Field_circ))
if Cin.type == 'ring'
    fprintf(' Total transmitted power %g [W] \n',tmp_power)
elseif Cin.type == 'folded'
    fprintf(' End transmitted power %g [W] \n',tmp_power)    
end
fprintf(' Reflected power %g [W] \n\n',Calculate_power(Cin.Field_ref))

if ~isempty(Cin.Laser_in.Field_SBl)
    
    [Pin1,Pin2] = Calculate_power_SB(Cin.Laser_in);
    [Pcirc1,Pcirc2] = Calculate_power_SB(Cin.Field_circ);
    
    Ptrans1 = 0; Ptrans2 = 0;
    for pp=2:Cin.Nb_mirror
        [tmp_power_lsb,tmp_power_usb] = Calculate_power_SB(Cin.Field_trans(pp));
        Ptrans1 = tmp_power_lsb + Ptrans1;
        Ptrans2 = tmp_power_usb + Ptrans2;
    end
    
    [Pref1 Pref2] = Calculate_power_SB(Cin.Field_ref);
    
    disp('---------- For the lower sideband ---------------')
    
    fprintf(' Power in the input beam %g [W] \n',Pin1)
    fprintf(' Circulating power %g [W] \n',Pcirc1)
    fprintf(' Transmitted power %g [W] \n',Ptrans1)
    fprintf(' Reflected power %g [W] \n\n',Pref1)
    
    disp('---------- For the upper sideband ---------------')
    
    fprintf(' Power in the input beam %g [W] \n',Pin2)
    fprintf(' Circulating power %g [W] \n',Pcirc2)
    fprintf(' Transmitted power %g [W] \n',Ptrans2)
    fprintf(' Reflected power %g [W] \n\n',Pref2)
    
end


figure(105)
clf;
subplot(2,2,1)
E_plot(Cin.Laser_in)
title('Input field')
subplot(2,2,2)
E_plot(Cin.Field_circ)
title('Circulating field')
subplot(2,2,3)
E_plot(Cin.Field_ref)
title('Reflected field')
subplot(2,2,4)
E_plot(Cin.Field_trans(end))
title('Transmitted field')






end