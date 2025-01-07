function I_out = Cut_Frequency_Interface(I_in,option,f_cut)
% Cut_frequency_Interface() filter in the spatial frequency domain an
% interface in a crude way (no windowing)
% option = 'LP' to low pass filter the data with a corner frequency f_cut
% option = 'HP' to high pass filter the data with a corner frequency f_cut
% option = 'BP' to band pass the data between the frequency f_cut(1) and
% f_cut(2)
% f_cut is given in m^-1

p = inputParser;

% Check if the first argument is an Interface
p.addRequired('Iin', @(x)isa(x, 'Interface'));

% Check if we display the amplitude or intensity
p.addRequired('option', @(x)(strcmpi(x,'LP') | ...
    strcmpi(x,'HP') | strcmpi(x,'BP')));

% Check if the scale is linear or logarithmic
p.addRequired('f_cut',@(x)isnumeric(x) && x>0);

p.parse(I_in,option,f_cut)

f_cut = p.Results.f_cut;
I_out = I_in;


% Cut the frequency
Grid_FFT2D = sqrt(I_in.Grid.D2_FFT_X.^2 + I_in.Grid.D2_FFT_Y.^2);

% Create a round filter filter
Filter_mir1 = ones(I_in.Grid.Num_point,I_in.Grid.Num_point);
Filter_mir2 = ones(I_in.Grid.Num_point,I_in.Grid.Num_point);

if strcmpi(p.Results.option,'LP')
    Filter_mir1(Grid_FFT2D > f_cut) = 0;
elseif strcmpi(p.Results.option,'HP')
    Filter_mir1(Grid_FFT2D < f_cut) = 0;
else
    if length(f_cut) == 2
        Filter_mir1(Grid_FFT2D < f_cut(1)) = 0;
        Filter_mir2(Grid_FFT2D > f_cut(2)) = 0;
    else
        error('Cut_frequency_Interface(): for band pass filtering 2 frequency are required: f_cut = [50 100]')
    end
end

I_out.surface = real(ifft2(ifftshift( Filter_mir1 .* Filter_mir2 .*fftshift(fft2(I_in.surface)) )));


end

