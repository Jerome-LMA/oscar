function Display_results(Cin)
% Display_results(Cin) display the results of the cavity calculations.
% The function Calculate_fields() must have been run first

if isempty(Cin.Field_ref)
    error('Display_results(): Before displaying the results, the function Calculate_fields() must be run     ')
end

disp('---------- For the carrier ---------------')

fprintf('  Power in the input beam: \t %6.6g \t [W] \n',Calculate_Power(Cin.Laser_in))
fprintf('  Circulating power: \t\t %6.6g \t [W] \n',Calculate_Power(Cin.Field_circ))
fprintf('  Transmitted power: \t\t %6.6g \t [W] \n',Calculate_Power(Cin.Field_trans))
fprintf('  Reflected power: \t\t\t %6.6g \t [W] \n',Calculate_Power(Cin.Field_ref))
fprintf('  Round trip losses [ppm]: \t %6.6g \n\n',Cin.Loss_RTL *1E6)

% Round trip losses calculation, using the energy conservation
% See virgo note VIR-0706A-10
%Loss = ( Calculate_power(Cin.Laser_in) - Calculate_power(Cin.Field_trans) - Calculate_power(Cin.Field_ref) )/ (Calculate_power(Cin.Field_circ) * (Cin.I_end.r)^2);

for ii=1:Cin.Laser_in.Nb_Pair_SB
    
    [Pin1, Pin2] = Calculate_Power(Cin.Laser_in,'include','SB','SB_num',ii);
    [Pcirc1, Pcirc2] = Calculate_Power(Cin.Field_circ,'include','SB','SB_num',ii);
    [Ptrans1, Ptrans2] = Calculate_Power(Cin.Field_trans,'include','SB','SB_num',ii);
    [Pref1, Pref2] = Calculate_Power(Cin.Field_ref,'include','SB','SB_num',ii);
    
    fprintf('---------- For the sidebands %i ---------------\n',ii)
    fprintf(' for the lower and upper sidebands respectively \n')
    fprintf(' Power in the input beam: \t %6g \t %6g \t [W] \n',Pin1,Pin2)
    fprintf(' Circulating power: \t\t %6g \t %6g \t [W] \n',Pcirc1,Pcirc2)
    fprintf(' Transmitted power: \t\t %6g \t %6g \t [W] \n',Ptrans1,Ptrans2)
    fprintf(' Reflected power: \t\t\t %6g \t %6g \t [W] \n\n',Pref1,Pref2)
       
end


%
%
% figure(105)
% clf;
% subplot(2,2,1)
% E_plot(Cin.Laser_in)
% title('Input field')
% subplot(2,2,2)
% E_plot(Cin.Field_circ)
% title('Circulating field')
% subplot(2,2,3)
% E_plot(Cin.Field_ref)
% title('Reflected field')
% subplot(2,2,4)
% E_plot(Cin.Field_trans)
% title('Transmitted field')
%
%
%



end