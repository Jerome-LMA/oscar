function Cout = calculate_rt_mat(Cin,varargin)
%calculate_rt_mat() Calculate the kernel for one round trip in the cavity
% C1 = calculate_rt_mat(C1), this function calculate the kernel for one
% round trip of the light in the cavity. From this kernel, one can derive
% the eigen modes and eigen vectors of the cavity.
% !! only use with small size of grid 64X64  with 4GB RAM, or 128X128 on
% more powerful machine

p  = inputParser;

% Check if the first argument is a cavity
p.addRequired('Cin', @(x)isa(x, 'CavityN' ));

% Check if a grid is given as a second argument
p.addOptional('Grid',[], @(x)isa(x, 'Grid' ));

% Check if we want to use the parralel computing toolbox or not
% By default off since it can take a huge amount of memory
p.addOptional('Use_PC',false, @(x)islogical(x));
p.parse(Cin,varargin{:})

if ~isempty(Cin.Cavity_EM_mat)
    disp('calculate_rt_mat(): Cavity kernel has already been calculated  ')
end

if isempty(p.Results.Grid)
    % No grid is given no resample is necessary
    resampled_grid = false;
    Gr = Cin.laser_in.Grid;
else
    resampled_grid = true;
    Gr = p.Results.Grid;
end

Cout = Cin;

if ~resampled_grid
    Num_point = Cin.laser_in.Grid.Num_point;
else
    if Cin.laser_in.Grid.Length ~= Gr.Length
        error('calculate_rt_mat(): The original grid and the new one have different lengths')
    end
    Num_point = Gr.Num_point;
end

tmp_mat_EM =  complex(zeros(Num_point^2));

for pp=1:Cin.Nb_mirror
    Cin.propagation_mat_array(pp).Use_DI = true;
    %Cin.propagation_mat_array(pp).Use_DI = false;
end


disp('Calculating the kernel:       ')

if license('test','distrib_computing_toolbox')  && p.Results.Use_PC         % check if the Parallel Computing Toolbox exists
    
    pool_obj = parpool;
    
    for mx=1:Num_point
        
        Tmp_RT_mat2 = complex(zeros(Num_point^2,Num_point));
        
        parfor py=1:Num_point
            
            if ~resampled_grid
                E_in = Cin.laser_in;
                E_in.Field = complex(zeros(Num_point));
                E_in.Field(mx,py) = 1;
                Field_Circ = E_in;
            else
                E_in = E_Field(Gr,'w0',450E-6); % dummy value
                E_in.Field = complex(zeros(Num_point));
                E_in.Field(mx,py) = 1;
                Field_Circ = Resample_E(E_in,Cin.laser_in.Grid);
            end
            
            for pp=1:Cin.Nb_mirror
                if pp ~= Cin.Nb_mirror % check we are not at the last iteration
                    Field_Circ = Propagate_E(Field_Circ,Cin.propagation_mat_array(pp));
                    Field_Circ = reflect_mirror(Field_Circ,Cin.I_array(pp+1),'Ref',1);
                else
                    Field_Circ = Propagate_E(Field_Circ,Cin.propagation_mat_array(pp));
                    Field_Circ = reflect_mirror(Field_Circ,Cin.I_array(1),'Ref',1);
                end
            end
            
            if ~resampled_grid
                Tmp_RT_mat2(:,py)  =  Field_Circ.Field(:);
            else
                Field_tmp = Resample_E(Field_Circ,Gr);
                Tmp_RT_mat2(:,py) = Field_tmp.Field(:);
            end
            
        end
        
        %tmp_mat_EM(:,(mx-1)*Num_point+1:mx*Num_point) = Tmp_RT_mat2;
        tmp_mat_EM(:,(0:Num_point-1)*Num_point + mx) = Tmp_RT_mat2;
        
    end
    
    delete(pool_obj);
    
else % if the PCT is not installed
    
    for mx=1:Num_point
        for py=1:Num_point
            
            if ~resampled_grid
                E_in = Cin.laser_in;
                E_in.Field = complex(zeros(Num_point));
                E_in.Field(mx,py) = 1;
                Field_Circ = E_in;
            else
                E_in = E_Field(Gr,'w0',450E-6); % dummy value
                E_in.Field = complex(zeros(Num_point));
                E_in.Field(mx,py) = 1;
                Field_Circ = Resample_E(E_in,Cin.laser_in.Grid);
            end
            
            for pp=1:Cin.Nb_mirror
                if pp ~= Cin.Nb_mirror % check we are not at the last iteration
                    Field_Circ = Propagate_E(Field_Circ,Cin.propagation_mat_array(pp));
                    Field_Circ = reflect_mirror(Field_Circ,Cin.I_array(pp+1),'Ref',1);
                else
                    Field_Circ = Propagate_E(Field_Circ,Cin.propagation_mat_array(pp));
                    Field_Circ = reflect_mirror(Field_Circ,Cin.I_array(1),'Ref',1);
                end
            end
            
            if ~resampled_grid
                tmp_mat_EM(:,(py-1)*Num_point+mx) = Field_Circ.Field(:);
            else
                Field_tmp = Resample_E(Field_Circ,Gr);
                tmp_mat_EM(:,(py-1)*Num_point+mx) = Field_tmp.Field(:);
                
            end
        end
        fprintf('\b\b\b\b\b\b\b %-3.0i %% ',round(100 * mx / Num_point))
    end
    fprintf('\b\b\b\b\b\b\b done! \n')
    
end

Cout.Cavity_EM_mat = tmp_mat_EM;




