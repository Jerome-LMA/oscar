clearvars; close all;
addpath(genpath([pwd filesep '..' filesep 'Classes']));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.31                                  ')
disp('  ')

% To be noted not all the functions are supported by the ABCD matrix but it
% can already help for some quick development

G1 = Grid(256,0.5);

E_input = E_Field(G1,'w',0.04,'R',-600);

disp('  ')
disp('Comparison of the input beam:')

Fit_TEM00(E_input)
Fit_TEM00(E_input,'With_ABCD',true)

disp('  ')
disp('Comparison after a propagation:')

E2 = Propagate_E(E_input,800);

Fit_TEM00(E2)
Fit_TEM00(E2,'With_ABCD',true)

disp('  ')
disp('Comparison after the transmission from a lens:')

E2 = Transmit_Lens(E_input,2E3/0.45);
E2 = Propagate_E(E2,100);

Fit_TEM00(E2)
Fit_TEM00(E2,'With_ABCD',true)
% 

disp('  ')
disp('Comparison after the transmission from a curved interface and the change of refractive index:')

I1 = Interface(G1,'RoC',-2E3);
E2 = Change_E_n(E_input,1.45);
[E3,E4] = Transmit_Reflect_Interface(E2,I1);
E3 = Propagate_E(E3,100);

E4 = Propagate_E(E4,100);
E4 = Change_E_n(E4,1);

Fit_TEM00(E3)
Fit_TEM00(E3,'With_ABCD',true)

disp('  ')
disp('Comparison after the reflection from a curved interface and the change of refractive index:')

Fit_TEM00(E4)
Fit_TEM00(E4,'With_ABCD',true)











