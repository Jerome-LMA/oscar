function  Cout = Declare_on_GPU(Cin)

Cout = Cin;
% Check if a GPU ressource is availabe
if ~parallel.gpu.GPUDevice.isAvailable
    error('Declare_on_GPU(): no GPU device available')
end

Cout.Run_on_GPU = true;

% Transfer the grid on the GPU, all the grids in the various instances of
% objects are linked to the same grid.
% Only do that for the relevant array; none found
Cout.Laser_in.Grid.Step_sq = gpuArray(Cout.Laser_in.Grid.Step_sq);


% Put the input beam in the GPU
Cout.Laser_in.Field = gpuArray(Cout.Laser_in.Field);

% Put the propagation matrix too
Cout.Propagation_mat.Use_GPU = true;
Cout.Propagation_mat.mat = gpuArray(Cout.Propagation_mat.mat);
Cout.Propagation_mat.mat_DI = gpuArray(Cout.Propagation_mat.mat_DI);

% Precompute also the wavefront distortions also in the GPU

Cout.I_input.WP_n1_GPU = exp(-1i * Cout.Laser_in.k_prop * Cout.I_input.surface *2) .* Cout.I_input.mask .* Cout.I_input.r;
Cout.I_input.WP_n1_GPU = gpuArray(Cout.I_input.WP_n1_GPU);

Cout.I_input.WP_n2_GPU = exp(1i * Cout.Laser_in.k_prop * Cout.I_input.surface *2) .* Cout.I_input.mask .* Cout.I_input.r;
Cout.I_input.WP_n2_GPU = gpuArray(Cout.I_input.WP_n2_GPU);

Cout.I_end.WP_n1_GPU = exp(-1i * Cout.Laser_in.k_prop * Cout.I_end.surface *2) .* Cout.I_end.mask .* Cout.I_end.r;
Cout.I_end.WP_n1_GPU = gpuArray(Cout.I_end.WP_n1_GPU);

Cout.I_end.WP_n2_GPU = exp(1i * Cout.Laser_in.k_prop * Cout.I_end.surface *2) .* Cout.I_end.mask .* Cout.I_end.r;
Cout.I_end.WP_n2_GPU = gpuArray(Cout.I_end.WP_n2_GPU);

end

