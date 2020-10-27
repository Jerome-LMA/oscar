function Cout = Cavity_lock_PDH(Cin)
%  Cout = Cavity_lock_PDH(Cin) Lock on the PDH signal
% SB must be present in the input field

% Calculate first the resonance length
% Must start outside the cavity

if Cin.Laser_start_on_input
    error(['Cavity_lock_PDH(' inputname(1) '): to calculate the PDH error signal, the beam must be defined outside the cavity, set Laser_start_on_input = false'])
end

if isempty(Cin.Laser_in.Field_SBl)
    error(['Cavity_lock_PDH(' inputname(1) '): to calculate the PDH error signal, a SB field must be present'])
end

Cout = Cin;

Scan.Nb_points = 5;
Scan.span = 0.001; % in rad around the resonance which maximise the circulating power

Scan_reso =  linspace(-Scan.span,Scan.span,Scan.Nb_points);
Scan_reso2 = Cin.Resonance_phase .* exp(1i*Scan_reso);
disp('PDF locking:       ')

for jj = 1:Scan.Nb_points
    Cin.Resonance_phase = Scan_reso2(jj);
    Cout_tmp = Calculate_fields(Cin,'accuracy',0.01);    
    Sig.p(jj) = Demodulate_SB(Cout_tmp.Field_ref);
    fprintf('\b\b\b\b\b\b %3.0d %%',round(jj /Scan.Nb_points * 100))
end

fprintf('\b\b\b\b\b\b done! \n')

% check one pos, one neg
if (sign(Sig.p(1)) * sign(Sig.p(end))) == 1
    disp('Warning: PDH locking() scanning range too small, result may be inacurate')
end

plot(Scan_reso,Sig.p)

p = polyfit(Scan_reso,Sig.p,1);

Offset = - p(2)/p(1);

Cout.Resonance_phase = Cout.Resonance_phase * exp(1i*Offset);
fprintf('PDH locking added phase offset [Rad]: %g \n',Offset)

end