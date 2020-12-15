clearvars; close all;
addpath(genpath(strcat(pwd, '\..\Classes')));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.21                                  ')
disp('  ')

%% Test Add_Tilt()
disp('----------   Test Add_Tilt()  -----------')

G1 = Grid(256,0.4);
E1 = E_Field(G1,'w',0.043,'R',-1034);
E2 = Add_Tilt(E1,1E-6,'dir','y');
Check_Position_Tilt(E2)

disp('')

%% Test Add_Sidebands()
disp('----------   Test Add_Sidebands()  -----------')

G1 = Grid(256,0.4);
E1 = E_Field(G1,'w',0.043,'R',-1034);
E2 = Add_Sidebands(E1,3.4E6,0.2);

disp('')

%% Test Calculate_Overlap()
disp('----------   Test Calculate_Overlap()  -----------')

G1 = Grid(256,0.4);
E1 = E_Field(G1,'w',0.043,'R',-1034);
E2 = Add_Sidebands(E1,3.4E6,0.2);

Calculate_Overlap(E2)
Calculate_Overlap(E1,E2)

disp('')

%% Test  Calculate_Power()
disp('----------   Test Calculate_Power  -----------')


G1 = Grid(256,0.4);
E1 = E_Field(G1,'w',0.043,'R',-1034);
E2 = Add_Sidebands(E1,3.4E6,0.2);

Calculate_Power(E1,'all')
[a,b] = Calculate_Power(E2,'SB','SB_num',1);

disp('')

%% Test E_Plot()
disp('----------   Test E_Plot()  -----------')


G1 = Grid(256,0.4);
E1 = E_Field(G1,'w',0.043,'R',-1034);
E2 = Add_Sidebands(E1,3.4E6,0.2);

E_Plot(E1)
E_Plot(E2,'all')
E_Plot(E2,'SB','SB_num',1,'display','intensity')

disp('')

%% Test Fit_TEM00() and other
disp('----------   Test Fit_TEM00(), Fit_TEM00xy, Fit_Efield  -----------')

G1 = Grid(256,0.4);
E1 = E_Field(G1,'w',0.043,'R',-1034);

Fit_TEM00(E1)
Fit_TEM00_xy(E1);


G1 = Grid(256,0.4);
E2 = E_Field(G1,'w',0.05,'R',-2000,'mode','HG 3 0');
Fit_E_Field(E2);

%% Test Focus_Beam_With_Telescope() Focus_Mirror();
disp('----------   Test Focus_Beam_With_Telescope()  -----------')

G1 = Grid(256,0.4);
E1 = E_Field(G1,'w',0.043,'R',inf);
E2 = Add_Sidebands(E1,3.4E6,0.2);

[E3,G3] = Focus_Beam_With_Telescope(E2,[100 102]);
Fit_TEM00(E3)
E_Plot(E3,'SB','SB_num',1,'display','intensity')

[E4,G4] = Focus_Mirror(E1,200,102);
Fit_TEM00(E4)

%%  Transmit_Aperture()
disp('----------   Test Transmit_Aperture()  -----------')

G1 = Grid(256,0.4);
E1 = E_Field(G1,'w',0.043,'R',inf);

E2 = Transmit_Aperture(E1,0.05);
figure(1); E_Plot(E2)

E3 = Transmit_Aperture(E1,0.05,'Shape','square');
figure(2); E_Plot(E3)

E4 = Transmit_Aperture(E1,0.2,'Shape','batman');
figure(3); E_Plot(E4)

%% Transmit_Lens()
disp('----------   Test Transmit_Lens()  -----------')

G1 = Grid(256,0.4);
E1 = E_Field(G1,'w',0.043,'R',inf);

E2 = Transmit_Lens(E1,200);
E2 = Propagate_E(E2,100);
Fit_TEM00(E2)

I1 = Interface(G1,'RoC',-0.45*200);
E2 = Transmit_Lens(E1,I1);
E2 = Propagate_E(E2,100);
Fit_TEM00(E2)

%% Transmit_Reflect_Interface()
disp('----------   Test Transmit_Reflect_Interface()  -----------')

G1 = Grid(256,0.4);
E1 = E_Field(G1,'w',0.043,'R',inf);
E1 = Add_Sidebands(E1,3.4E6,0.2);
I1 = Interface(G1,'RoC',-200);

E2 = Transmit_Reflect_Interface(E1,I1);
E2 = Propagate_E(E2,100);
Fit_TEM00(E2)

%% Transmit_Reflect_Mirror()
disp('----------   Test  Transmit_Reflect_Mirror()  -----------')

G1 = Grid(256,0.4);
E1 = E_Field(G1,'w',0.043,'R',inf);
E1 = Add_Sidebands(E1,3.4E6,0.2);

IM_AR = Interface(G1,'RoC',-1420,'CA',0.33,'T',1-100E-6,'L',0,'n1',1,'n2',1.45);
IM_HR = Interface(G1,'RoC',1420,'CA',0.33,'T',0.014,'L',0,'n1',1,'n2',1.45);
IM = Mirror(IM_HR,IM_AR,0.2);

[E2,E3] = Transmit_Reflect_Optic(E1,IM,'AR');

Fit_TEM00(E2)
Fit_TEM00(E3)

%% Transmit_Reflect_Optic()
disp('----------   Test  Transmit_Reflect_Optic()  -----------')

G1 = Grid(256,0.4);
E1 = E_Field(G1,'w',0.043,'R',inf);
E1 = Add_Sidebands(E1,3.4E6,0.2);

IM_AR = Interface(G1,'RoC',-1420,'CA',0.33,'T',1-100E-6,'L',0,'n1',1,'n2',1.45);
IM_HR = Interface(G1,'RoC',1420,'CA',0.33,'T',0.014,'L',0,'n1',1,'n2',1.45);
IM = Mirror(IM_HR,IM_AR,0.2);

[E2,E3] = Transmit_Reflect_Optic(E1,IM,'AR');
[E2,E3] = Transmit_Reflect_Optic(E1,IM_HR);


%[E2,E3] = Transmit_Reflect_Optic(E1,IM);