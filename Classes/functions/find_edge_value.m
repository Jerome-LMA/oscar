function edge_value = find_edge_value(map)
% find_edge_value(): return the average values on the edge of a map. That
% avoids a jump on the edge of the map

[m n] = size(map);

%figure(2); imagesc(map); axis square
edgevalue_tmp = [];

for ii=1:m
    ind1 = find(~isnan(map(ii,:)),10,'first');
    ind2 = find(~isnan(map(ii,:)),10,'last');
    
    if ~isempty(ind1)
        edgevalue_tmp(end+1:end+length(ind1)) = map(ii,ind1);
    end
    
    if ~isempty(ind2)
        edgevalue_tmp(end+1:end+length(ind2)) = map(ii,ind2);
    end
    
end

%plot(edgevalue_tmp)
edge_value = median(edgevalue_tmp);

if isempty(edge_value)
    edge_value = 0;
    disp('Warning: Add_map() did not manage to find the edge values for the map')
    
end
end

