function  [] = declare_on_gpu(obj)


% Check if a GPU ressource is availabe
if ~parallel.gpu.GPUDevice.isAvailable
    error('Declare_on_GPU(): no GPU device available')
end

obj.run_on_gpu = true;

% Transfer the grid on the GPU, all the grids in the various instances of
% objects are linked to the same grid.
% Only do that for the relevant array; none found
obj.laser_in.Grid.Step_sq = gpuArray(obj.laser_in.Grid.Step_sq);


% Put the input beam in the GPU
obj.laser_in.Field = gpuArray(obj.laser_in.Field);

% Put the propagation matrix too
obj.propagation_mat.Use_GPU = true;
obj.propagation_mat.mat = gpuArray(obj.propagation_mat.mat);
obj.propagation_mat.mat_DI = gpuArray(obj.propagation_mat.mat_DI);

% Precompute also the wavefront distortions also in the GPU

obj.i_input.WP_n1_GPU = exp(-1i * obj.laser_in.k_prop * obj.i_input.surface *2) .* obj.i_input.mask .* obj.i_input.r;
obj.i_input.WP_n1_GPU = gpuArray(obj.i_input.WP_n1_GPU);

obj.i_input.WP_n2_GPU = exp(1i * obj.laser_in.k_prop * obj.i_input.surface *2) .* obj.i_input.mask .* obj.i_input.r;
obj.i_input.WP_n2_GPU = gpuArray(obj.i_input.WP_n2_GPU);

obj.i_end.WP_n1_GPU = exp(-1i * obj.laser_in.k_prop * obj.i_end.surface *2) .* obj.i_end.mask .* obj.i_end.r;
obj.i_end.WP_n1_GPU = gpuArray(obj.i_end.WP_n1_GPU);

obj.i_end.WP_n2_GPU = exp(1i * obj.laser_in.k_prop * obj.i_end.surface *2) .* obj.i_end.mask .* obj.i_end.r;
obj.i_end.WP_n2_GPU = gpuArray(obj.i_end.WP_n2_GPU);

end

