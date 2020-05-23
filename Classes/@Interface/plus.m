function [Iout] = 	plus(I1,I2)
%Plus: Overload the + function for the object of class interface
%simply add the interface together
% !!! keep the properties (mask, reflectivity,...) of  I1 and only add the
% surface of I2

if (I1.Grid ~= I2. Grid)
    error('plus(): error the 2 inputs interfaces are defined on a different grid')
end


% Check if the two masks have the same size
if (I1.mask ~= I2.mask)
    disp('Warning: addition of 2 maps with different clear aperture')
end


Iout = I1;
Iout.mask = I1.mask.*I2.mask;

Iout.surface = I1.surface + I2.surface;

end

