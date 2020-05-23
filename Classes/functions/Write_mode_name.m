function str_out = Write_mode_name(family,m,n)
%Write_mode_name() = write a string with the name of the mode
% family = 'HG' or 'LG'   m,n order of the mode

if strcmp(family, 'HG') || strcmp(family, 'LG') 
str_out = [family ' ' num2str(m) ' ' num2str(n)];

else
    error('The first argument must be either LG or HG')
end

end

