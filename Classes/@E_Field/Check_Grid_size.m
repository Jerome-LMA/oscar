function Check_Grid_size(Ein,tol_out)
% Check if the power of the beam is concentrated in the middle of the grid
% (expected case) or not. In the latter case, that means the user has to
% increase the grid size
% 'tol_out' is the relative maximum tolerated amount of power outside the central
% diameter which 90% of the size of the grid.

if isa(Ein, 'E_Field') && isa(tol_out, 'numeric')
    total_power = Calculate_power(Ein);
        
    mask_index = Ein.Grid.D2_r < (Ein.Grid.Length * 0.45);
    mask = ones(Ein.Grid.Num_point,Ein.Grid.Num_point,'double');
    mask(mask_index) = 0;
    
    Ein.Field = Ein.Field .* mask;
    power_outside = Calculate_power(Ein);
    
    if (power_outside/total_power) > tol_out
        warning('war:Grid_too_small','Grid size likely to be too small, check the results carefully')
        %E_plot(Ein)
    end
    
else
    error('Check_Grid_size(): wrong input argument')
end

end


