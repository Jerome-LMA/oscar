function [length_scan, power_scan] = scan(obj, varargin)
% Cout = Cavity_scan(obj) scan the cavity to find the position of the
% maximum power
% Scan the cavity over one FSR and then find the suitable resonance length
% to maximise the circulating power

p  = inputParser;

% Check if we use the parallel computing toolbox (or not)
p.addParameter('use_parallel',true,@(x)islogical(x));

% Check if we save the scan in a file
p.addParameter('save_scan',false,@(x) islogical(x) || ischar(x));

% Check if we have to zoom to calculate and save the resonance length
p.addParameter('Define_L_length',false,@(x)islogical(x));

% Check if we all include the SB in the calculation
p.addParameter('With_SB',false,@(x)islogical(x));

% Scan over how many FSR
p.addParameter('Nb_FSR',1,@(x)isscalar(x));


p.parse(varargin{:})

% Define some variables
Grid_num_point = obj.laser_in.Grid.Num_point;

if Grid_num_point > 256
    answer = input('Scanning a cavity with a large grid can take excessive time and memory \n Do you want to continue ? Y/N [N]: ', 's');
    if isempty(answer)
        answer = 'N';
    end
    if ~strcmp(answer,'Y')
        error('Scan aborted by the user')
    end
end

if isempty(obj.cavity_scan_all_field)
    obj.propagate_field();
end

% Save the scan
obj.cavity_scan_all_field = obj.cavity_scan_all_field;

tmp = size(obj.cavity_scan_all_field);
num_iter = tmp(3);
num_point_scan = obj.cavity_scan_param(1);

% Following needed if SB, only work for the first SB:
if obj.laser_in.Nb_Pair_SB
    D_phi = (2*pi*obj.laser_in.SB(1).Frequency_Offset/2.99792E8) * obj.Length;
    if obj.laser_in.Nb_Pair_SB > 1
        disp('Only the first pair of SB is taking into account for the cavity scan')
    end
else
    D_phi = 0;
end

% Define where we store the results of the scan
power_scan = zeros(1,num_point_scan,'double');

% Create the length vector to scan the cavity
length_scan = (1:num_point_scan) * obj.laser_in.Wavelength/num_point_scan * round(abs(p.Results.Nb_FSR));

fprintf(' Scanning the cavity ...       ')

if license('test','distrib_computing_toolbox') && p.Results.use_parallel          % check if the Parallel Computing Toolbox exists
    
    if gpuDeviceCount > 0 && parallel.gpu.GPUDevice.isAvailable(1) && isgpuarray(obj.laser_in.Field)
        
        disp('Found suitable GPU. Starting GPU-based scan.')
        gq = 1:num_point_scan;
        ii = 1:num_iter;
        gcavity_scan_all_field_arr = gpuArray(obj.cavity_scan_all_field(:,:,ii));
        gcavity_scan_all_field_arr_perm = permute(gcavity_scan_all_field_arr, [3,1,2]);
        gPhase_shifts = gpuArray(exp(1i*obj.laser_in.k_prop* length_scan(gq)'*ii));
        gFields_reconstructed = pagefun(@mtimes,gPhase_shifts, gcavity_scan_all_field_arr_perm);
        Fields_reconstructed = gather(gFields_reconstructed);
        
        for qqq = 1:num_point_scan
            Dummy_E = obj.laser_in;
            Dummy_E.Field = squeeze(Fields_reconstructed(qqq,:,:));
            power_scan(qqq) = calculate_power(Dummy_E);
        end
        
    else
        
        pool_obj = gcp('nocreate');
        if (isempty(pool_obj))
            disp('Parallel pool not initialized. Starting now...')
            is_par_pool_init = false;
            pool_obj = gcp();
        else
            is_par_pool_init = true;
        end
        
        parfor qq = 1:num_point_scan
            Field_reconstructed = complex(zeros(Grid_num_point,Grid_num_point,'double'));
            Field_reconstructed_SBu = Field_reconstructed; % do it even if no SB
            Field_reconstructed_SBl = Field_reconstructed;
            
            for ii=1:num_iter
                Field_reconstructed = Field_reconstructed + obj.cavity_scan_all_field(:,:,ii) * exp(1i*obj.laser_in.k_prop* length_scan(qq)*ii);
                if obj.laser_in.Nb_Pair_SB
                    Field_reconstructed_SBu = Field_reconstructed_SBu + obj.cavity_scan_all_field(:,:,ii) * exp(1i*obj.laser_in.k_prop* length_scan(qq)*ii) * exp(1i*2*D_phi*ii);
                    Field_reconstructed_SBl = Field_reconstructed_SBl + obj.cavity_scan_all_field(:,:,ii) * exp(1i*obj.laser_in.k_prop* length_scan(qq)*ii) * exp(-1i*2*D_phi*ii);
                end
            end
            
            % Create a dummy field with the right field inside and normalise
            % the SB by the input power
            Dummy_E = obj.laser_in;
            Dummy_E.Field = Field_reconstructed;
            if obj.laser_in.Nb_Pair_SB
                Dummy_E.SB(1).Field_upper = Field_reconstructed_SBu * sqrt(calculate_power(obj.laser_in,'SB')/2);
                Dummy_E.SB(1).Field_lower = Field_reconstructed_SBl * sqrt(calculate_power(obj.laser_in,'SB')/2);
            end
            
            if p.Results.With_SB
                power_scan(qq) = calculate_power(Dummy_E,'include','all');
            else
                power_scan(qq) = calculate_power(Dummy_E);
            end
        end
        
        if (~is_par_pool_init)
            disp('Shutting down pool because it was not initialized at startup.')
            delete(pool_obj);
        end
    end
    
else % if the PCT is not installed or we do not want to use the toolbox
    
    for qq = 1:num_point_scan
        
        Field_reconstructed = complex(zeros(Grid_num_point,Grid_num_point,'double'));
        Field_reconstructed_SBu = Field_reconstructed; % do it even if no SB
        Field_reconstructed_SBl = Field_reconstructed;
        
        for ii=1:num_iter
            Field_reconstructed = Field_reconstructed + obj.cavity_scan_all_field(:,:,ii) * exp(1i*obj.laser_in.k_prop* length_scan(qq)*ii);
            if obj.laser_in.Nb_Pair_SB
                Field_reconstructed_SBu = Field_reconstructed_SBu + obj.cavity_scan_all_field(:,:,ii) * exp(1i*obj.laser_in.k_prop* length_scan(qq)*ii) * exp(1i*D_phi*ii);
                Field_reconstructed_SBl = Field_reconstructed_SBl + obj.cavity_scan_all_field(:,:,ii) * exp(1i*obj.laser_in.k_prop* length_scan(qq)*ii) * exp(-1i*D_phi*ii);
            end
        end
        
        % Create a dummy field with the right field inside and normalise
        % the SB by the input power
        Dummy_E = obj.laser_in;
        Dummy_E.Field = Field_reconstructed;
        if obj.laser_in.Nb_Pair_SB
            Dummy_E.SB(1).Field_upper = Field_reconstructed_SBu * sqrt(calculate_power(obj.laser_in,'SB')/2);
            Dummy_E.SB(1).Field_lower = Field_reconstructed_SBl * sqrt(calculate_power(obj.laser_in,'SB')/2);
        end
        
        if p.Results.With_SB
            power_scan(qq) = calculate_power(Dummy_E,'include','all');
        else
            power_scan(qq) = calculate_power(Dummy_E);
        end
        
        if (rem(qq,num_point_scan/100) == 0)
            fprintf('\b\b\b\b\b\b\b\b\b   %-3.0i %% ',100*qq/num_point_scan)
        end
    end
    
end

fprintf('Finished ... \n')

if p.Results.save_scan
    tmp_save(:,1) = length_scan;
    tmp_save(:,2) = power_scan;
    save(['Cavity_scan_' p.Results.save_scan '.txt'],'tmp_save','-ASCII');
end
%
% figure(101)
% semilogy(length_scan,power_scan)
% title('Cavity scan over one FSR')


obj.cavity_scan_r(:,1) = length_scan;
obj.cavity_scan_r(:,2) = power_scan;

if p.Results.Define_L_length
    %------------ Zoom on the maximum ---------------------
    
    [~,ind_max] = max(power_scan);
    max_pos = length_scan(ind_max);
    
    num_point_scan = obj.cavity_scan_param(2);
    span_scan = obj.cavity_scan_param(3);
    
    
    %Define the zoom length vector
    power_scan = zeros(1,num_point_scan,'double');
    
    
    % Create the length vector to scan the cavity
    length_scan = max_pos - span_scan/2 + ...
        (1:num_point_scan)*span_scan/num_point_scan;
    
    fprintf(' Zooming on the resonance peak ...       ')
    
    
    if license('test','distrib_computing_toolbox') && p.Results.use_parallel       % check if the Parallel Computing Toolbox exists
        
        if gpuDeviceCount > 0 && parallel.gpu.GPUDevice.isAvailable(1) && isgpuarray(obj.laser_in.Field)
            
            disp('Found suitable GPU. Starting GPU-based scan.')
            gq = 1:num_point_scan;
            ii = 1:num_iter;
            gcavity_scan_all_field_arr = gpuArray(obj.cavity_scan_all_field(:,:,ii));
            gcavity_scan_all_field_arr_perm = permute(gcavity_scan_all_field_arr, [3,1,2]);
            gPhase_shifts = gpuArray(exp(1i*obj.laser_in.k_prop* length_scan(gq)'*ii));
            gFields_reconstructed = pagefun(@mtimes,gPhase_shifts, gcavity_scan_all_field_arr_perm);
            Fields_reconstructed = gather(gFields_reconstructed);
            
            for qqq = 1:num_point_scan
                Dummy_E = obj.laser_in;
                Dummy_E.Field = squeeze(Fields_reconstructed(qqq,:,:));
                power_scan(qqq) = calculate_power(Dummy_E);
            end
            
        else
            
            pool_obj = parpool;
            
            parfor qq = 1:num_point_scan
                Field_reconstructed = complex(zeros(Grid_num_point,Grid_num_point,'double'));
                for ii=1:num_iter
                    Field_reconstructed = Field_reconstructed + obj.cavity_scan_all_field(:,:,ii) * exp(1i*obj.laser_in.k_prop* length_scan(qq)*ii);
                end
                Dummy_E = obj.laser_in;
                Dummy_E.Field = Field_reconstructed;
                power_scan(qq) = calculate_power(Dummy_E);
            end
            
            delete(pool_obj);
        end
        
    else % The PCT is not installed
        
        for qq = 1:num_point_scan
            Field_reconstructed = complex(zeros(Grid_num_point,Grid_num_point,'double'));
            for ii=1:num_iter
                Field_reconstructed = Field_reconstructed + obj.cavity_scan_all_field(:,:,ii) * exp(1i*obj.laser_in.k_prop* length_scan(qq)*ii);
            end
            Dummy_E = obj.laser_in;
            Dummy_E.Field = Field_reconstructed;
            
            power_scan(qq) = calculate_power(Dummy_E);
            
            if (rem(qq,num_point_scan/100) == 0)
                fprintf('\b\b\b\b\b\b\b\b\b   %-3.0i %% ',100*qq/num_point_scan)
            end
        end
    end
    
    
    fprintf('Finished ... \n')
    
    
    % figure(102)
    % semilogy(length_scan,power_scan)
    % title('Zoom on the resonance peak')
    
    [~,ind_max] = max(power_scan);
    max_pos = length_scan(ind_max);
    
    %
    % fprintf(' Microscopique resonance length: %-9.7d ... \n',max_pos);
    %
    
    obj.cavity_scan_rz(:,1) = length_scan;
    obj.cavity_scan_rz(:,2) = power_scan;
    
    % Find the additional round trip phase to put the cavity on resonance
    obj.resonance_phase = exp(1i*obj.laser_in.k_prop* max_pos);
    
end

end
