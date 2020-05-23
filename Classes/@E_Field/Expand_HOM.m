function varargout  = Expand_HOM(Ein,max_mode_order,varargin)
% Expand_HOM, give the power of a beam in the HOM up to a certain order
% Expand_HOM(Ein,max_mode_order) expand the beam Ein into the HG basis
% higher order modes
% The type of screen output can be set by the variable 'display' and we can
% also override the basis of the mode for the caculation with 'basis'.
% The function can also return the scattering matrix into HOM with:
% A = Expand_HOM(Ein,max_mode);

p = inputParser;
p.FunctionName = 'Expand an E_field into higher order modes';

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check if the second argument is a scalar
p.addRequired('max_mode_order', @(x)isscalar(x) && x>0);

% Check if we display the results as a vector or pyramid
p.addOptional('display',[], @(x)strcmpi(x,'vector') | ...
    strcmpi(x,'pyramid'));

% Can fix arbitrary the basis, give the beam radius and the wavefront
% curvature of the basis mode (in meter)
p.addOptional('basis',[], @(x)isvector(x) && length(x)==2);


p.parse(Ein,max_mode_order,varargin{:})

Max_HOM = round(p.Results.max_mode_order);

% Expand on the basis of the perfect Gaussian beam, so find first the beam
% parameters:
if isempty(p.Results.basis)
[Fit_w Fit_R] = Fit_TEM00(Ein);
else
    Fit_w = p.Results.basis(1);
    Fit_R = p.Results.basis(2);
end

Matrix_scattering = zeros(Max_HOM +1);

for jj_mode_order = 0:Max_HOM
    for ii = 1:jj_mode_order + 1;
        
        mode_name = ['HG ' num2str(ii-1) ' ' num2str(jj_mode_order - ii+1) ];
        Coeff_overlap = Calculate_Overlap(E_Field(Ein.Grid,'w',Fit_w,'R',Fit_R,'mode',mode_name),Ein);
        
        Matrix_scattering(jj_mode_order+1,ii) = abs(Coeff_overlap)^2;
        
    end
end

for pp = 1:Max_HOM+1
    Power_per_mode(pp) = sum(Matrix_scattering(pp,1:pp));
end

% Do the display now
if ~isempty(p.Results.display)
    disp('Power content:')
    if strcmp(p.Results.display,'vector')
        for pp = 1:Max_HOM+1
            tmp_name = ['Mode order ' num2str(pp-1) '     '];
            tmp_str = num2str(Power_per_mode(pp), ' %2.3d,');
            disp([tmp_name tmp_str] )
        end
        
    elseif strcmp(p.Results.display,'pyramid')
        for pp = 1:Max_HOM+1
            tmp_name = ['Mode order ' num2str(pp-1) '     '];
            tmp_str = num2str(Matrix_scattering(pp,1:pp), ' %2.3d,');
            disp([tmp_name tmp_str] )
        end
        
    else
        disp('Expand_HOM(): wrong argument for the display  ')
        
    end
end

if nargout == 0
    semilogy(0:Max_HOM, Power_per_mode,'s','MarkerSize',10,'MarkerFaceColor','b')
    xlim([-0.4 Max_HOM+0.4])
    xlabel(' Mode order')
    ylabel(' Power content   ')
elseif nargout == 1
    varargout(1) = {Matrix_scattering};
elseif nargout == 2    
    varargout(1) = {Matrix_scattering};
    varargout(2) = {Power_per_mode};
else
    disp('Expand_HOM(): wrong number of output argument')
end

end

