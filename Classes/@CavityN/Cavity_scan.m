function Cout = Cavity_Scan(Cin)
% Cout = Cavity_scan(Cin) scan the cavity to find the position of the
% maximum power
% Scan the cavity over one FSR and then find the suitable resonance length
% to maximise the circulating power


Cout = Cin;

% Define some variables
Grid_num_point = Cin.Laser_in.Grid.Num_point;

if Grid_num_point > 256
    answer = input('Scanning a cavity with a large grid can take excessive time and memory \n Do you want to continue ? Y/N [N]: ', 's');
    if isempty(answer)
        answer = 'N';
    end
    if ~strcmp(answer,'Y')
        error('Scan aborted by the user')
    end
end


if isempty(Cin.Cavity_scan_all_field)
    Cin = Cavity_propagate_field(Cin);
end

tmp = size(Cin.Cavity_scan_all_field);
num_iter = tmp(3);
num_point_scan = Cin.Cavity_scan_param(1);

% Define where we store the results of the scan
Length_scan = zeros(1,num_point_scan,'double');
Power_scan = zeros(1,num_point_scan,'double');

% Create the length vector to scan the cavity
Length_scan = (1:num_point_scan) * Cin.Laser_in.Wavelength/num_point_scan;

fprintf(' Scanning the cavity ...       ')

if license('test','distrib_computing_toolbox')         % check if the Parallel Computing Toolbox exists
    
    pool_obj = parpool;
    
    parfor qq = 1:num_point_scan
        Field_reconstructed = complex(zeros(Grid_num_point,Grid_num_point,'double'));
        for ii=1:num_iter
            Field_reconstructed = Field_reconstructed + Cin.Cavity_scan_all_field(:,:,ii) * exp(1i*Cin.Laser_in.k_prop* Length_scan(qq)*ii);
        end
        
        % Create a dummy field with the right field inside
        Dummy_E = Cin.Laser_in;
        Dummy_E.Field = Field_reconstructed;
        Power_scan(qq) = Calculate_power(Dummy_E);
    end
    
    delete(pool_obj);
    
else % if the PCT is not installed
    
    
    for qq = 1:num_point_scan
        Field_reconstructed = complex(zeros(Grid_num_point,Grid_num_point,'double'));
        for ii=1:num_iter
            Field_reconstructed = Field_reconstructed + Cin.Cavity_scan_all_field(:,:,ii) * exp(1i*Cin.Laser_in.k_prop* Length_scan(qq)*ii);
        end
        
        % Create a dummy field with the right field inside
        Dummy_E = Cin.Laser_in;
        Dummy_E.Field = Field_reconstructed;
        
        Power_scan(qq) = Calculate_power(Dummy_E);
        
        if (rem(qq,num_point_scan/100) == 0)
            fprintf('\b\b\b\b\b\b\b\b\b   %-3.0i %% ',100*qq/num_point_scan)
        end
    end
    
end



fprintf('Finished ... \n')

% tmp_save(:,1) = Length_scan;
% tmp_save(:,2) = Power_scan;
% save(['Cavity_scan_' inputname(1) '.txt'],'tmp_save','-ASCII');
%
% figure(101)
% semilogy(Length_scan,Power_scan)
% title('Cavity scan over one FSR')

Cout.Cavity_scan_R(:,1) = Length_scan;
Cout.Cavity_scan_R(:,2) = Power_scan;

%------------ Zoom on the maximum ---------------------

[~,ind_max] = max(Power_scan);
max_pos = Length_scan(ind_max);

num_point_scan = Cin.Cavity_scan_param(2);
span_scan = Cin.Cavity_scan_param(3);


%Define the zoom length vector
Length_scan = zeros(1,num_point_scan,'double');
Power_scan = zeros(1,num_point_scan,'double');


% Create the length vector to scan the cavity
Length_scan = max_pos - span_scan/2 + ...
    (1:num_point_scan)*span_scan/num_point_scan;

fprintf(' Zooming on the resonance peak ...       ')


if exist('matlabpool','builtin')          % check if the Parallel Computing Toolbox exists
    
    pool_obj = parpool;
    
    parfor qq = 1:num_point_scan
        Field_reconstructed = complex(zeros(Grid_num_point,Grid_num_point,'double'));
        for ii=1:num_iter
            Field_reconstructed = Field_reconstructed + Cin.Cavity_scan_all_field(:,:,ii) * exp(1i*Cin.Laser_in.k_prop* Length_scan(qq)*ii);
        end
        Dummy_E = Cin.Laser_in;
        Dummy_E.Field = Field_reconstructed;
        Power_scan(qq) = Calculate_power(Dummy_E);
    end
    
    delete(pool_obj);
    
    
else % The PCT is not installed
    
    for qq = 1:num_point_scan
        Field_reconstructed = complex(zeros(Grid_num_point,Grid_num_point,'double'));
        for ii=1:num_iter
            Field_reconstructed = Field_reconstructed + Cin.Cavity_scan_all_field(:,:,ii) * exp(1i*Cin.Laser_in.k_prop* Length_scan(qq)*ii);
        end
        Dummy_E = Cin.Laser_in;
        Dummy_E.Field = Field_reconstructed;
        
        Power_scan(qq) = Calculate_power(Dummy_E);
        
        if (rem(qq,num_point_scan/100) == 0)
            fprintf('\b\b\b\b\b\b\b\b\b   %-3.0i %% ',100*qq/num_point_scan)
        end
    end
end


fprintf('Finished ... \n')


% figure(102)
% semilogy(Length_scan,Power_scan)
% title('Zoom on the resonance peak')

[~,ind_max] = max(Power_scan);
max_pos = Length_scan(ind_max);

%
% fprintf(' Microscopique resonance length: %-9.7d ... \n',max_pos);
%

Cout.Cavity_scan_RZ(:,1) = Length_scan;
Cout.Cavity_scan_RZ(:,2) = Power_scan;

% Find the additional round trip phase to put the cavity on resonance
Cout.Resonance_phase = exp(1i*Cin.Laser_in.k_prop* max_pos);

Cout.Cavity_scan_all_field = Cin.Cavity_scan_all_field;

end
