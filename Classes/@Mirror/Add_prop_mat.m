function M2 = Add_prop_mat(M1,E1)
% Pre-calculate the propagation matrix inside the substrate.
% used to speed up the calculations
% E1 is taken outside the mirror, it is used to get the wavelength of the
% main laser beam

M2 = M1;

M2.Propagation_mat_sub = Prop_operator(E1,M1.length_substrate,'n',M1.n_substrate);
  
end

