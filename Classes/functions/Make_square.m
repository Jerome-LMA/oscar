function map_out = Make_square(map_in)
%Make_square() Take a map and make it square by cutting the excess lines or
%row. Keep the map centered

[m,n] = size(map_in);
% m number of line, n number of column

map_out = map_in;

if m > n
    To_cut = m - n;
    if rem(To_cut,2) % odd number of line to cut
        cut1  = (To_cut-1)/2;
        cut2 =  cut1;
        
        map_out(1:cut1,:) = [];
        map_out(end:-1:end-cut2,:) = [];
        
    else % even number of line to cut
        cut1  = To_cut/2;
        cut2 =  cut1-1;
        
        map_out(1:cut1,:) = [];
        map_out(end:-1:end-cut2,:) = [];
        
    end
    
elseif n > m
    To_cut = n - m;
    if rem(To_cut,2) % odd number of line to cut
        cut1  = (To_cut-1)/2;
        cut2 =  cut1;
        
        map_out(:,1:cut1) = [];
        map_out(:,end:-1:end-cut2) = [];
        
    else % even number of line to cut
        cut1  = To_cut/2;
        cut2 =  cut1-1;
        
        map_out(:,1:cut1) = [];
        map_out(:,end:-1:end-cut2) = [];
        
    end
end

end

