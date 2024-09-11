function PSD_reconstructed = FitFunctionPSD(freq_data, power_law)
% Return the parameterised PSD 1D function with several segments
% using the set of parameters P and the frequency axis xdata

% P contains the coefficients for the fit, if 2 coefficients:
% PSD = P(1) * 1 / f^(-P(2))    P(1) is the amplitude, P(2) is the power
% law

% if 4 coefficient: P(1) is the amplitude, P(2) the cutting frequency
% between the 2 segments, P(3) the first segment power law, P(4) the second
% segment power law
% for freq < P(2)  PSD = P(1) * 1 / f^(-P(3)),
% for freq > P(2)  PSD = Amp * 1 / f^(-P(4)) with of course the continuity in the PSD
%and so on, for an arbitrary number of segment

% The frequency must be ordered in increasing value

power_law(1) = power_law(1)*1E-21; % for the fit to be easier, ensure that the fitting coefficient are in the same range


if(rem(length(power_law),2) ~= 0)
    warning('P has to be an even number of parameters')
end

Nb_segment = length(power_law)/2;

Freq_cut1 = power_law(1);
Freq_cut_last = power_law(Nb_segment);

% if Freq_cut1 < freq_data(1)
%    Freq_cut1
%    freq_data(1)
%     error('FitFunctionPSD(): very unlikely error: contact the developer')
% end

if Freq_cut_last > freq_data(end)
    warning('FitFunctionPSD(): the power law of the PSD covers a wider spatial range that what is allowed by the grid size and resolution')
end


PSD_reconstructed = zeros(1,length(freq_data));
%PSD_reconstructed = 1E-30*ones(1,length(fdata));

if Nb_segment == 1
    PSD_reconstructed = power_law(1) * freq_data.^(-power_law(2));
    
elseif Nb_segment == 2
    
    jj = 1; % index to scan the frequency vector
    
    while freq_data(jj) <=  power_law(2)
        PSD_reconstructed(jj) = power_law(1) * freq_data(jj).^(-power_law(3));
        jj = jj + 1;
        if jj > length(freq_data)
            return
        end
    end
    
    Amp_2 = power_law(1) * power_law(2).^(power_law(4) - power_law(3)); % Calculate the starting amplitude of the second segment to have the continuity
    
    for jj = jj:length(freq_data)
        PSD_reconstructed(jj) = Amp_2 * freq_data(jj).^(-power_law(4));
    end
    
else
    jj = 1; % index to scan the frequency vector
    
    while freq_data(jj) <=  power_law(2)
        PSD_reconstructed(jj) = power_law(1) * freq_data(jj).^(-power_law(Nb_segment+1));
        jj = jj + 1;
        if jj > length(freq_data)
            return
        end
    end
    
    
    for ii = 2:Nb_segment-1
        Amp(ii) = PSD_reconstructed(jj-1) * power_law(ii).^(power_law(Nb_segment+ii));
        
        while freq_data(jj) <=  power_law(ii+1)
            PSD_reconstructed(jj) = Amp(ii) * freq_data(jj).^(-power_law(Nb_segment+ii));
            jj = jj + 1;
            if jj > length(freq_data)
                return
            end
        end
    end
    
    % for the last segment
    Amp(ii+1) = PSD_reconstructed(jj-1) * power_law(Nb_segment).^(power_law(end));
    
    for jj = jj:length(freq_data)
        PSD_reconstructed(jj) = Amp(ii+1) * freq_data(jj).^(-power_law(end));
    end
    
end

PSD_reconstructed = PSD_reconstructed';  % to match the arrangement during the fitting

end
