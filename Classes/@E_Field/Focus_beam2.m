function [Eout Gout] = Focus_beam2(Ein,RofC_mir,dist)
%  [Eout Gout] = Focus_beam(Ein,RofC_mir,dist)
% Strongly focus a beam, for that need to resize the grid at the output
% Ein = Input E_field, define just before the mirror
% RofC_mir = Radius of curvature of the mirror
% dist = distance from the mirror where Eout will be defined (usually up to
% the next optic).

%!! The distance of propagation after the mirror must be less than half the
%mirror RofC

% For explanation see Hiro's talk: G1000001-v1

P_alpha =1 - 2*dist*1/RofC_mir;
%P_alpha = 1 - 2*dist*(1/RofC_mir + 1/(2*-2500));

% 2D wavefront change phi_in
WF_change = exp(1i * Ein.k_prop * ((P_alpha - 1)/(2*dist) + 1/RofC_mir ) *  Ein.Grid.D2_square );     
% figure(101)
% imagesc(abs(WF_change))

% Field E_0 tilt in the new coordinate system is:
Eout = Ein;
Eout.Field = Ein.Field .* WF_change;

% Propagate the new field over (z1 - Z0) / alpha
Eout = Propagate_E(Eout,dist/P_alpha);

Gout = Grid(Ein.Grid.Num_point,Ein.Grid.Length*abs(P_alpha));

Eout.Grid = Gout; 

% Multiply the results by C tilt
Eout.Field = Eout.Field .* (1/P_alpha) .* exp(1i * Ein.k_prop * ((1/P_alpha - 1)/(2*dist))  *  Eout.Grid.D2_square );     

end

