function PSD_reconstructed = FitFunctionPSD(P,fdata)
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
% P(1) = P(1)*1E-16; % for the fit to be easier, ensure that the fitting coefficient are in the same range
% and so on, for an arbitrary number of segment

if(rem(length(P),2) ~= 0)
    warning('P has to be an even number of parameters')
end

Nb_segment = length(P)/2;

PSD_reconstructed = zeros(1,length(fdata));
%PSD_reconstructed = 1E-30*ones(1,length(fdata));



if Nb_segment == 1
    PSD_reconstructed = P(1) * fdata.^(-P(2));
    
elseif Nb_segment == 2
    
    jj = 1; % index to scan the frequency vector
    
    while fdata(jj) <=  P(2)
        PSD_reconstructed(jj) = P(1) * fdata(jj).^(-P(3));
        jj = jj + 1;
    end
    
    Amp_2 = P(1) * P(2).^(P(4) - P(3)); % Calculate the starting amplitude of the second segment to have the continuity
    
    for jj = jj:length(fdata)
        PSD_reconstructed(jj) = Amp_2 * fdata(jj).^(-P(4));
        
    end
    
else
    jj = 1; % index to scan the frequency vector
    
    while fdata(jj) <=  P(2)
        PSD_reconstructed(jj) = P(1) * fdata(jj).^(-P(Nb_segment+1));
        jj = jj + 1;
    end
    
    
    for ii = 2:Nb_segment-1
        Amp(ii) = PSD_reconstructed(jj-1) * P(ii).^(P(Nb_segment+ii));
        
        while fdata(jj) <=  P(ii+1)
            PSD_reconstructed(jj) = Amp(ii) * fdata(jj).^(-P(Nb_segment+ii));
            jj = jj + 1;
        end
    end
    
    % for the last segment
    Amp(ii+1) = PSD_reconstructed(jj-1) * P(Nb_segment).^(P(end));
    
    for jj = jj:length(fdata)
        PSD_reconstructed(jj) = Amp(ii+1) * fdata(jj).^(-P(end));
    end
    
end

PSD_reconstructed = PSD_reconstructed';  % to match the arrangement during the fitting

end
