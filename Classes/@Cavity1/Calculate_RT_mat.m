function Cout = Calculate_RT_mat(Cin,varargin)
%Calculate_RT_mat() Calculate the kernel for one round trip in the cavity
% C1 = Calculate_RT_mat(C1), this function calculate the kernel for one
% round trip of the light in the cavity. From this kernel, one can derive
% the eigen modes and eigen vectors of the cavity.
% !! only use with small size of grid 64X64  with 4GB RAM, or 128X128 on
% a more powerful machine

p  = inputParser;

% Check if the first argument is a cavity
p.addRequired('Cin', @(x)isa(x, 'Cavity1'));

% Check if we want to use the parralel computing toolbox or not
% By default off since it can take a huge amount of memory
p.addOptional('Use_PC',false, @(x)islogical(x));

p.parse(Cin,varargin{:})

if ~isempty(Cin.Cavity_EM_mat)
    disp('Calculate_RT_mat(): Cavity kernel has already been calculated  ')
end

Cout = Cin;

Num_point = Cin.Laser_in.Grid.Num_point;
tmp_mat_EM =  complex(zeros(Num_point^2));

Cin.Propagation_mat.Use_DI = true;

disp('Calculating the kernel:       ')

if license('test','distrib_computing_toolbox')  && p.Results.Use_PC        % check if the Parallel Computing Toolbox exists and we want to use
    
    pool_obj = parpool;
    
    for mx=1:Num_point
        
        Tmp_RT_mat2 = complex(zeros(Num_point^2,Num_point));
        
        parfor py=1:Num_point
            
            E_in = Cin.Laser_in;
            E_in.Field = complex(zeros(Num_point));
            E_in.Field(mx,py) = 1;
            
            Circ_field = Propagate_E(E_in,Cin.Propagation_mat);
            Circ_field = Reflect_mirror(Circ_field,Cin.I_end,'Ref',1);
            Circ_field = Propagate_E(Circ_field,Cin.Propagation_mat);
            Circ_field = Reflect_mirror(Circ_field,Cin.I_input,'Ref',1);
            
            Tmp_RT_mat2(:,py)  =  Circ_field.Field(:);
            
        end
        
        tmp_mat_EM(:,(0:Num_point-1)*Num_point + mx) = Tmp_RT_mat2;

    end
    
    delete(pool_obj);

else % if the PCT is not installed or we do not want to use it
    
    for mx=1:Num_point
        for py=1:Num_point
            
            E_in = Cin.Laser_in;
            E_in.Field = complex(zeros(Num_point));
            E_in.Field(mx,py) = 1;
            
            Circ_field = Propagate_E(E_in,Cin.Propagation_mat);
            Circ_field = Reflect_mirror(Circ_field,Cin.I_end,'Ref',1);
            Circ_field = Propagate_E(Circ_field,Cin.Propagation_mat);
            Circ_field = Reflect_mirror(Circ_field,Cin.I_input,'Ref',1);
            
            tmp_mat_EM(:,(py-1)*Num_point+mx) = Circ_field.Field(:);       
        end
            fprintf('\b\b\b\b\b\b\b %-3.0i %% ',round(100 * mx / Num_point))            
    end
 fprintf('\b\b\b\b\b\b\b done! \n')            

    
end

Cout.Cavity_EM_mat = tmp_mat_EM;

end

