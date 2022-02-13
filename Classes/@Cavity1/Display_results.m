function display_results(Cin,varargin)
% Display_results(Cin) display the results of the cavity calculations.
% The function Calculate_fields() must have been run first

p  = inputParser;

% Check if the first argument is an E_Field
p.addRequired('Cin', @(x)isa(x, 'Cavity1'));

% Display or not the E_field
p.addParameter('display',true,@(x)isa(x,'logical'));

p.parse(Cin,varargin{:})

if isempty(Cin.Field_ref)
    error('Display_results(): Before displaying the results, the function Calculate_fields() must be run     ')
end

if Cin.laser_in.Nb_Pair_SB % if there is no SB no need to mention for the carrier
    disp('---------- For the carrier ---------------')
end

fprintf('  Power in the input beam: \t %6.6g \t [W] \n',calculate_power(Cin.laser_in))
fprintf('  Circulating power: \t\t %10.6g \t [W] \n',calculate_power(Cin.Field_circ))
fprintf('  Transmitted power: \t\t %10.6g \t [W] \n',calculate_power(Cin.Field_trans))
fprintf('  Reflected power: \t\t\t %10.6g \t [W] \n',calculate_power(Cin.Field_ref))
fprintf('  Round trip losses: \t\t %10.6g \t [ppm] \n\n',Cin.Loss_RTL *1E6)

% Round trip losses calculation, using the energy conservation
% See virgo note VIR-0706A-10
%Loss = ( Calculate_power(Cin.laser_in) - Calculate_power(Cin.Field_trans) - Calculate_power(Cin.Field_ref) )/ (Calculate_power(Cin.Field_circ) * (Cin.i_end.r)^2);

for ii=1:Cin.laser_in.Nb_Pair_SB
    
    [Pin1, Pin2] = calculate_power(Cin.laser_in,'include','SB','SB_num',ii);
    [Pcirc1, Pcirc2] = calculate_power(Cin.Field_circ,'include','SB','SB_num',ii);
    [Ptrans1, Ptrans2] = calculate_power(Cin.Field_trans,'include','SB','SB_num',ii);
    [Pref1, Pref2] = calculate_power(Cin.Field_ref,'include','SB','SB_num',ii);
    
    fprintf('---------- For the sidebands %i, frequency: %5.4g [MHz] ---------------\n',ii,Cin.laser_in.SB(ii).Frequency_Offset/1E6)
    fprintf(' for the lower and upper sidebands respectively \n')
    fprintf(' Power in the input beam: \t %6g \t %6g \t [W] \n',Pin1,Pin2)
    fprintf(' Circulating power: \t\t %6g \t %6g \t [W] \n',Pcirc1,Pcirc2)
    fprintf(' Transmitted power: \t\t %6g \t %6g \t [W] \n',Ptrans1,Ptrans2)
    fprintf(' Reflected power: \t\t\t %6g \t %6g \t [W] \n\n',Pref1,Pref2)
    
end


if p.Results.display
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
    E_Plot(Cin.Field_trans)
    title('Transmitted field')
end




end