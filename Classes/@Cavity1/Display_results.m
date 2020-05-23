function Display_results(Cin)
% Display_results(Cin) display the results of the cavity calculations.
% The function Calculate_fields() must have been run first

if isempty(Cin.Field_ref)
    error('Display_results(): Before displaying the results, the function Calculate_fields() must be run     ')
end

if ~isempty(Cin.Laser_in.Field_SBl)
    disp('---------- For the carrier ---------------')
end

fprintf(' Power in the input beam %g [W] \n',Calculate_power(Cin.Laser_in))
fprintf(' Circulating power %g [W] \n',Calculate_power(Cin.Field_circ))
fprintf(' Transmitted power %g [W] \n',Calculate_power(Cin.Field_trans))
fprintf(' Reflected power %g [W] \n\n',Calculate_power(Cin.Field_ref))

% Round trip losses calculation, using the energy conservation
% See virgo note VIR-0706A-10
%Loss = ( Calculate_power(Cin.Laser_in) - Calculate_power(Cin.Field_trans) - Calculate_power(Cin.Field_ref) )/ (Calculate_power(Cin.Field_circ) * (Cin.I_end.r)^2);


fprintf('Round trip losses [ppm]: %g \n',Cin.Loss_RTL *1E6)


if ~isempty(Cin.Laser_in.Field_SBl)
    
    [Pin1, Pin2] = Calculate_power_SB(Cin.Laser_in);
    [Pcirc1, Pcirc2] = Calculate_power_SB(Cin.Field_circ);
    [Ptrans1, Ptrans2] = Calculate_power_SB(Cin.Field_trans);
    [Pref1, Pref2] = Calculate_power_SB(Cin.Field_ref);
    
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
E_plot(Cin.Field_trans)
title('Transmitted field')






end