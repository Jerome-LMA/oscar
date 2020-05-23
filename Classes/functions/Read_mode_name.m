function [family,m,n] = Read_mode_name(name)
%Read_mode_name() = takes the name and returns the family name and the mode numbers
% family = 'HG' or 'LG'   m,n order of the mode

ind_space = strfind(name,' ');

family = name(1:ind_space(1)-1);

m = str2double( name(ind_space(1)+1:ind_space(2)-1) );
n = str2double( name(ind_space(2)+1:end) );


end