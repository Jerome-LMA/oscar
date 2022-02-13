function display_results(Cin)
% Display_results(Cin) display the results of the cavity calculations.
% The function Calculate_fields() must have been run first

if isempty(Cin.Field_ref)
    error('Display_results(): Before displaying the results, the function Calculate_fields() must be run     ')
end

if Cin.laser_in.Nb_Pair_SB % if there is no SB no need to mention for the carrier
    disp('---------- For the carrier ---------------')
end


% Calculate the total transmitted power
tmp_power = 0;

if Cin.type == 'ring'
    for pp=2:Cin.Nb_mirror
        tmp_power = tmp_power + calculate_power(Cin.Field_trans(pp));
    end
elseif Cin.type == 'folded'
    tmp_power = calculate_power(Cin.Field_trans(end));
end

fprintf(' Power in the input beam %g [W] \n',calculate_power(Cin.laser_in))
fprintf(' Circulating power %g [W] \n',calculate_power(Cin.Field_circ))
if Cin.type == 'ring'
    fprintf(' Total transmitted power %g [W] \n',tmp_power)
elseif Cin.type == 'folded'
    fprintf(' End transmitted power %g [W] \n',tmp_power)
end
fprintf(' Reflected power %g [W] \n\n',calculate_power(Cin.Field_ref))

for ii=1:Cin.laser_in.Nb_Pair_SB
    
    [Pin1, Pin2] = calculate_power(Cin.laser_in,'include','SB','SB_num',ii);
    [Pcirc1, Pcirc2] = calculate_power(Cin.Field_circ,'include','SB','SB_num',ii);
    
    Ptrans1 = 0; Ptrans2 = 0;
    for pp=2:Cin.Nb_mirror
        [tmp_power_lsb,tmp_power_usb] = calculate_power(CCin.Field_trans(pp),'include','SB','SB_num',ii);
        Ptrans1 = tmp_power_lsb + Ptrans1;
        Ptrans2 = tmp_power_usb + Ptrans2;
    end
    
    [Pref1, Pref2] = calculate_power(Cin.Field_ref,'include','SB','SB_num',ii);
    
    fprintf('---------- For the sidebands %i ---------------\n',ii)
    fprintf(' for the lower and upper sidebands respectively \n')
    fprintf(' Power in the input beam: \t %6g \t %6g \t [W] \n',Pin1,Pin2)
    fprintf(' Circulating power: \t\t %6g \t %6g \t [W] \n',Pcirc1,Pcirc2)
    fprintf(' Transmitted power: \t\t %6g \t %6g \t [W] \n',Ptrans1,Ptrans2)
    fprintf(' Reflected power: \t\t\t %6g \t %6g \t [W] \n\n',Pref1,Pref2)
    
end


figure(105)
clf;
subplot(2,2,1)
E_Plot(Cin.laser_in)
title('Input field')
subplot(2,2,2)
E_Plot(Cin.Field_circ)
title('Circulating field')
subplot(2,2,3)
E_Plot(Cin.Field_ref)
title('Reflected field')
subplot(2,2,4)
E_Plot(Cin.Field_trans(end))
title('Transmitted field')






end