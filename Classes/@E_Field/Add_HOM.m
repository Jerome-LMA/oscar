function varargout  = Add_HOM(Ein,max_mode_order,varargin)
% Add_HOM: add some higher order mode to a laser beam
% the input beam is a TEM00 and we will add the mode in the same base
p = inputParser;
p.FunctionName = 'Add some higher order modes';

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check if the second argument is a scalar
p.addRequired('max_mode_order', @(x)isscalar(x) && x>0);

p.parse(Ein,max_mode_order,varargin{:})
Max_HOM = round(p.Results.max_mode_order);

% higher order mode distribution
% (1 = all the mode with the same amplitude)
% (2 = lower mode have more amplitude, not implemented)
HOM_distri = 2;

[Fit_w,Fit_R] = Fit_TEM00(Ein);

Nb_mode = 0;
for ii = 0:max_mode_order
    Nb_mode = Nb_mode + (ii+1);
end

if HOM_distri == 1
    Ampli_HOM = 0.05 * ones(1,Max_HOM);
end

if HOM_distri == 2
    Ampli_HOM(1:max_mode_order) = 1./(1:max_mode_order);
    Ampli_HOM = Ampli_HOM / sum(Ampli_HOM);
    Ampli_HOM = Ampli_HOM.^2;
end

for jj_mode_order = 1:Max_HOM
    for ii = 1:jj_mode_order + 1;
        mode_name = ['HG ' num2str(ii-1) ' ' num2str(jj_mode_order - ii+1) ];
        Ein = Ein + Ampli_HOM(jj_mode_order) * E_Field(Ein.Grid,'w',Fit_w,'R',Fit_R,'mode',mode_name);
       % E_plot(E_Field(Ein.Grid,'w',Fit_w,'R',Fit_R,'mode',mode_name)); pause;
    end
end


if nargout == 0
    E_plot(Ein)
    return
elseif nargout == 1
    varargout(1) = {Ein};
else
    disp('Add_HOM(): wrong number of output argument')
end

end

