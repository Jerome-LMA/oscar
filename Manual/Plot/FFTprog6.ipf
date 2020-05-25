#pragma rtGlobals=1		// Use modern global access method.

Macro createEField()

	//NewPath/O 	FFT_Folder,"C:Documents and Settings:demonstrator:My Documents:Jerome1:Temp_result:FFT_data"

	Variable Power,Radius
	Variable/G Num_point, Length_Grid, Lambda, k_prop, step_grid
	Variable/G Length_cav,RofC_ITM,RofC_ETM,R_ITM,R_ETM,Refrac_index,Waist
	Variable i,j,radius_temp_sq,amplitude_E
	Variable X_dist,Y_dist
	Variable/C dist_temp
	Variable RofC_ITM_defor,RofC_ETM_defor
	Variable/G Return_power
	Variable RofC_ITM_tmp,RofC_ETM_tmp
	
	Num_point = 128						//  Number of points for the field matrix
	Length_Grid = 0.30						//  Side length of the grid where the calculations will be done (in meter)
	Lambda = 1064e-9						//  Wavelength of the laser beam (m)
	k_prop = (2 * Pi)/Lambda 					//  Cte propagation
	step_grid = Length_Grid/Num_point			//  Step for the integration
	
	Length_cav = 80+ 1.91709e-09				// length of the cavity (m)
	RofC_ITM = 2e20							// radius of curvature ITM (m)
	RofC_ETM = 1E3							// radius of curvature ETM (m)
	R_ITM = 0.99816							// (power) reflectivity ITM 
	R_ETM = 0.99884						// (power) reflectivity ETM 
	Refrac_index = 1.75						// Optical refractive index of the TM substrate
	
	RofC_ITM_defor = 247					// Radius of curvature of the ITM induced deformation
	RofC_ETM_defor = 3.44e5					// Radius of curvature of the ETM induced deformation
	
	//RofC_ITM = 1/(1/RofC_ITM + 1/RofC_ITM_defor)
	//RofC_ETM = 1/(1/RofC_ETM - 1/RofC_ETM_defor)
		
	//Waist = 0.010							// Gingin Test 1 waist = 0.00855308
	Waist = 0.04 						//  Beam radius (i.e. waist) in meter
	Radius = 1e20							//  Curvature of the laser beam in meter
	//Radius = 100
	Power = 1								//  Power of the laser beam in watt
	
	PauseUpdate; Silent 1
	Make/O/C/N = (Num_point,Num_point)/D Field_Start
	
	SetScale/P x -Length_Grid/2+step_grid/2,step_grid,"", Field_Start
	SetScale/P y -Length_Grid/2+step_grid/2,step_grid,"", Field_Start
	
	amplitude_E = 1
	
	i = 0
	Do
		j = 0
		Do
			X_dist = - Length_Grid/2 + step_grid/2 + i * step_grid
			Y_dist = - Length_Grid/2 + step_grid/2 + j * step_grid
			
			radius_temp_sq = X_dist^2 + Y_dist ^2
			dist_temp = amplitude_E * Exp(-radius_temp_sq/Waist^2) * Exp(- sqrt(-1) *  k_prop * radius_temp_sq/(2 * Radius))
			Field_Start[i][j] = dist_temp
			j = j + 1	
		While ( j < Num_point)			 
		i = i+ 1	
	While ( i < Num_point)			
	
	
	Make/O/N = (Num_point)/D Grid_axis
	i = 0
	Do 
			Grid_axis[i] = -Length_Grid/2+step_grid/2 + i*step_grid
			i = i + 1	
	While ( i < Num_point)
	
	Make/O/N = (Num_point)/D Grid_axis_FFT
	i = 0
	Do 
			Grid_axis_FFT[i] = -1/(2*step_grid) + i * 1/(Num_point  * step_grid)
			i = i + 1	
	While ( i < Num_point)
	



End



Macro Calculate_power(Wave_name)
	String Wave_name

	Variable i,j
	 Variable power_temp
	 PauseUpdate; Silent 1
	
	power_temp = 0
	i = 0
	Do
		j = 0
		Do			
			power_temp = power_temp + magsqr($Wave_name[i][j]) * step_grid^2
			j = j + 1	
		While ( j < Num_point)			 
		i = i+ 1	
	While ( i < Num_point)
	
	Print "Laser beam power (W):"
	Print power_temp
	Return_power = power_temp
	
End

Macro Beam_parameter(Wave_name)
	String Wave_name

	Variable i,j,mil_temp
	Variable Phasetemp
	Variable/G Waist_beam_temp
	PauseUpdate; Silent 1
	
	Make/O/N = (Num_point,Num_point)/D Amplitude_beam
	
	SetScale/P x -Length_Grid/2+step_grid/2,step_grid,"", Amplitude_beam
	SetScale/P y -Length_Grid/2+step_grid/2,step_grid,"", Amplitude_beam

	i = 0
	Do
		j = 0
		Do				
			Amplitude_beam[i][j] = cabs($Wave_name[i][j])
			j = j + 1	
		While ( j < Num_point)			 
		i = i+ 1	
	While ( i < Num_point)
	
	mil_temp = Num_point/2
	
	Make/O/N = (Num_point)/D Cross_sec_amp
	Make/O/N = (Num_point)/D Cross_sec_phase
	
	j = 0
	Do				
		Cross_sec_amp[j] = Amplitude_beam[mil_temp][j]
		Phasetemp = imag(r2polar($Wave_name[mil_temp][j]))
		Cross_sec_phase[j] = Phasetemp
		j = j + 1	
	While ( j < Num_point)

	CurveFit/Q gauss Cross_sec_amp /X=Grid_axis /D
	Print "Waist of the beam:"
	Print W_coef[3]
	Waist_beam_temp = W_coef[3]
	
	Unwrap 6.28319, Cross_sec_phase 
	
	Make/O/N = (Num_point)/D Waveexp
	
	i = 0
	Do
		Waveexp[i] = exp(-2*Grid_axis[i]*Grid_axis[i]/(W_coef[3]*W_coef[3]))
		i = i +1
	While(i<Num_point)
	
	CurveFit/Q poly 3, Cross_sec_phase /X=Grid_axis /W=Waveexp /D  

	Print "Radius of the beam:"
	Print -k_prop/(2*W_coef[3])
	 	 
End


Proc Make_propagation(Wave_name,Distance)
	String Wave_name
	Variable Distance
	
	Variable i,j
	
	PauseUpdate; Silent 1
	
	Make/O/C/N = (Num_point,Num_point)/D Mat_propagation
	
	FFT $Wave_name
	
	i = 0
	Do
		j = 0
		Do				
			Mat_propagation[i][j] = Exp(sqrt(-1) * ( - k_prop*Distance + Pi *  Lambda * (Grid_axis_FFT[i]^2 + Grid_axis_FFT[j]^2) * Distance ) ) 
			j = j + 1	
		While ( j < Num_point)			 
		i = i+ 1	
	While ( i < Num_point)
	
	$Wave_name = $Wave_name * Mat_propagation
	IFFT /C $Wave_name
	
End

Proc Prop_Square()

	Variable i,j,radius_temp_sq,amplitude_E
	Variable X_dist,Y_dist
	PauseUpdate; Silent 1
	
	Make/O/C/N = (Num_point,Num_point)/D Field_Square
	
	SetScale/P x -Length_Grid/2+step_grid/2,step_grid,"", Field_Square
	SetScale/P y -Length_Grid/2+step_grid/2,step_grid,"", Field_Square
	
	i = 0
	Do
		j = 0
		Do
			X_dist = - Length_Grid/2 + step_grid/2 + i * step_grid
		       Y_dist = - Length_Grid/2 + step_grid/2 + j * step_grid
			
			radius_temp_sq = X_dist^2 + Y_dist ^2
		
			if ( (i > Num_point*3/8) &  (i < Num_point*5/8) )
				if ( (j > Num_point*3/8) &  (j < Num_point*5/8) )
					Field_Square[i][j] = 1
				endif
			endif
			j = j + 1	
		While ( j < Num_point)			 
		i = i+ 1	
	While ( i < Num_point)
	
	Make_propagation("Field_Square")

	SetScale/P x -Length_Grid/2+step_grid/2,step_grid,"", Field_test
	SetScale/P y -Length_Grid/2+step_grid/2,step_grid,"", Field_test

End

Proc Create_mirror()

	Variable i,j,radius_temp_sq,X_dist, Y_dist
	Variable RofC

	Make/O/C/N = (Num_point,Num_point)/D Mirror_ITM_trans, Mirror_ITM_cav,Mirror_ETM_cav,Mirror_ITM_out
	
	// Mirror_ITM_in grid for the beam entering the cavity
	// Mirror_ITM_cav grid for the beam reflecting in the cavity
	
	SetScale/P x -Length_Grid/2+step_grid/2,step_grid,"", Mirror_ITM_cav
	SetScale/P y -Length_Grid/2+step_grid/2,step_grid,"", Mirror_ITM_cav

	PauseUpdate; Silent 1

	i = 0
	Do
		j = 0
		Do
			X_dist = - Length_Grid/2 + step_grid/2 + i * step_grid
		      Y_dist = - Length_Grid/2 + step_grid/2 + j * step_grid
			
			radius_temp_sq = X_dist^2 + Y_dist ^2
			
			Mirror_ITM_trans[i][j] = -(RofC_ITM - sqrt(RofC_ITM^2 - radius_temp_sq)) 
			Mirror_ITM_cav[i][j] = -(RofC_ITM - sqrt(RofC_ITM^2 - radius_temp_sq)) * 2
			Mirror_ETM_cav[i][j] = -(RofC_ETM - sqrt(RofC_ETM^2 - radius_temp_sq)) * 2
			j = j + 1	
		While ( j < Num_point)			 
		i = i+ 1	
	While ( i < Num_point)
	
	Mirror_ITM_out = (-1) * Mirror_ITM_cav
	
End

Proc Create_mirror_ANSYS()

	Variable i,j
	Variable X_dist,Y_dist,radius_temp

	PauseUpdate; Silent 1

	Duplicate/O Grid_axis ITM_trans_extra,ITM_cav_extra,ETM_cav_extra

	Interpolate/T=2/N=100/E=2/I=3/A=20/J=2/Y=ITM_trans_extra/X=Grid_axis FFT_ITM_Trans /X=Wave_x
	Interpolate/T=2/N=100/E=2/I=3/A=20/J=2/Y=ITM_cav_extra/X=Grid_axis FFT_ITM_Cav /X=Wave_x
	Interpolate/T=2/N=100/E=2/I=3/A=20/J=2/Y=ETM_cav_extra/X=Grid_axis FFT_ETM_Cav /X=Wave_x

	Make/O/N = (Num_point,Num_point)/D Mirror_ITM_trans, Mirror_ITM_cav,Mirror_ETM_cav, Mirror_ITM_out

	SetScale/P x -Length_Grid/2+step_grid/2,step_grid,"", Mirror_ITM_trans
	SetScale/P y -Length_Grid/2+step_grid/2,step_grid,"", Mirror_ITM_trans
	
	i = 0
	Do
		j = 0
		Do
			X_dist = - Length_Grid/2 + step_grid/2 + i * step_grid
		      Y_dist = - Length_Grid/2 + step_grid/2 + j * step_grid
			
			radius_temp = sqrt(X_dist^2 + Y_dist ^2)
			
			if (radius_temp > Length_Grid/2)
				Mirror_ITM_trans[i][j] = 0
				Mirror_ITM_cav[i][j] = 0
				Mirror_ETM_cav[i][j] = 0
			else
				Mirror_ITM_trans[i][j] = interp(radius_temp, Grid_axis, ITM_trans_extra )
				Mirror_ITM_cav[i][j] = interp(radius_temp, Grid_axis, ITM_cav_extra )
				Mirror_ETM_cav[i][j] =interp(radius_temp, Grid_axis, ETM_cav_extra )
			endif
			j = j + 1	
		While ( j < Num_point)			 
		i = i+ 1	
	While ( i < Num_point)

	Mirror_ITM_out = (-1) * Mirror_ITM_cav

End

Proc Propa_mirror(Wave_Field,Wave_mirror,reflec)
	String Wave_Field,Wave_mirror
	Variable reflec

	Variable i,j
	Variable/C Dephas_temp

	PauseUpdate; Silent 1

	i = 0
	Do
		j = 0
		Do
			Dephas_temp = Exp ( -sqrt(-1) * k_prop * $Wave_mirror[i][j] )
			$Wave_Field[i][j] = $Wave_Field[i][j]  * Dephas_temp * reflec
			j = j + 1	
		While ( j < Num_point)			 
		i = i+ 1	
	While ( i < Num_point)
End

Macro Test_Sphere()
	
	Print Secs2Time(DateTime,0)
	
	Variable i,j
	Variable iter_max,r_ITM_amp,r_ETM_amp,t_ITM_amp,t_ETM_amp
	
	PauseUpdate; Silent 1
	iter_max = 10000

	createEField()

	printf "RofC ITM = %g and RofC ETM %g\r", RofC_ITM, RofC_ETM

	Make/O/C/N = (Num_point,Num_point)/D Field_Cav,Field_circ,Field_in,Field_out,Field_temp
	Make/O/N = (Num_point,Num_point)/D Field_dephas
	
	Make/O/N = (iter_max)/D Tab_beam_waist
	
	Field_in = Field_Start
	Beam_parameter("Field_in")
	
	Field_Cav = 0
	Field_out = 0
	
	Create_mirror()
	
	// Find the resonant length of the cavity
	Field_circ = Field_Start
	Make_propagation("Field_circ",Length_cav)
	Propa_mirror("Field_circ","Mirror_ETM_cav",-1)
	Make_propagation("Field_circ",Length_cav)
	Propa_mirror("Field_circ","Mirror_ITM_cav",1)
	
	Field_dephas  = imag(r2polar(Field_circ/Field_Start))
	
	Length_cav = Length_cav + Field_dephas[Num_point/2][Num_point/2]/k_prop/2
	
	Print Field_dephas[Num_point/2][Num_point/2]/k_prop/2
	
	r_ITM_amp = sqrt(R_ITM)
	r_ETM_amp = sqrt(R_ETM)
	t_ITM_amp = sqrt(1- R_ITM)		// transmission of the ITM
	t_ETM_amp = sqrt(1 - R_ETM)
	
	Field_temp = Field_in
	Propa_mirror("Field_temp","Mirror_ITM_trans",t_ITM_amp)
	Field_circ = Field_temp
	
	Field_temp = Field_in	
	Propa_mirror("Field_temp","Mirror_ITM_out",-r_ITM_amp)
	Field_out = Field_temp
	
	i = 0
	Do
		//Beam_parameter("Field_circ")
		//Tab_beam_waist[i] = Waist_beam_temp
		Field_Cav =  Field_Cav + Field_circ 
		Make_propagation("Field_circ",Length_cav )
		Propa_mirror("Field_circ","Mirror_ETM_cav",-r_ETM_amp)
		Make_propagation("Field_circ",Length_cav )
		Field_temp = Field_circ
		Propa_mirror("Field_temp","Mirror_ITM_trans",t_ITM_amp)
		Field_out = Field_out + Field_temp
		Propa_mirror("Field_circ","Mirror_ITM_cav",r_ITM_amp)
		i = i+ 1	
		print i
	While ( i < iter_max)

	Print "Power field inside the cavity"
	Calculate_power("Field_Cav")

	Print "Power field out"
	Calculate_power("Field_out")
	
	Print Secs2Time(DateTime,0)

End

Macro Test_ANSYS()
	
	Print Secs2Time(DateTime,0)
	
	Variable i,j
	Variable iter_max,r_ITM_amp,r_ETM_amp,t_ITM_amp,t_ETM_amp
	
	PauseUpdate; Silent 1
	iter_max = 800

	createEField()

	Make/O/C/N = (Num_point,Num_point)/D Field_Cav,Field_circ,Field_in,Field_out,Field_temp
	Make/O/N = (Num_point,Num_point)/D Field_dephas
	
	
	Field_in = Field_Start
	Beam_parameter("Field_in")
	
	Field_Cav = 0
	Field_out = 0
	
	Create_mirror_ANSYS()
	
	r_ITM_amp = sqrt(R_ITM)
	r_ETM_amp = sqrt(R_ETM)
	t_ITM_amp = sqrt(1 - R_ITM)
	t_ETM_amp = sqrt(1 - R_ETM)
	
	
	//Length_cav = 71.99999986048559 +0e-12
	
	Field_temp = Field_in
	Propa_mirror("Field_temp","Mirror_ITM_trans",t_ITM_amp)
	Field_circ = Field_temp
	
	Field_temp = Field_in	
	Propa_mirror("Field_temp","Mirror_ITM_out",-r_ITM_amp)
	Field_out = Field_temp
	
	i = 0
	Do
		Field_Cav =  Field_Cav + Field_circ 
		Make_propagation("Field_circ",Length_cav )
		Propa_mirror("Field_circ","Mirror_ETM_cav",-r_ETM_amp)
		Make_propagation("Field_circ",Length_cav )
		Field_temp = Field_circ
		Propa_mirror("Field_temp","Mirror_ITM_trans",t_ITM_amp)
		Field_out = Field_out + Field_temp
		Propa_mirror("Field_circ","Mirror_ITM_cav",r_ITM_amp)
		i = i+ 1	
		print i
	While ( i < iter_max)

	
	Calculate_power("Field_Cav")

	Calculate_power("Field_out")
	
	Print Secs2Time(DateTime,0)

End


Macro Analyze_output()

	Variable i,j,mil_temp,Power_field_out_factor,temp_phase
	Variable/C Conv_temp1,Conv_temp2,Conv_temp3,

	PauseUpdate; Silent 1
	mil_temp = Num_point/2
	Make/O/N = (Num_point)/D Cross_sec_output
	Make/O/N = (Num_point)/D Cross_sec_input
	Make/O/N = (Num_point)/D Fit_outputfield
	
	j = 0
	Do				
//		Cross_sec_output[j] = cabs(Field_out[mil_temp][j]) *sign(sin(imag(r2polar(Field_out[mil_temp][j]))))
		Cross_sec_output[j] = cabs(Field_out[mil_temp][j])
		Cross_sec_input[j] = cabs(Field_in[mil_temp][j])
//		print j
//		print sin(imag(r2polar(Field_out[mil_temp][j])))
		j = j + 1	
	While ( j < Num_point)
	
	Display Cross_sec_output vs Grid_axis	
		
	K0 = 0;K2 = 0;K3 =Waist 	;
	CurveFit/H="1011" gauss Cross_sec_output /X=Grid_axis /D=Fit_outputfield  /R
	
	Power_field_out_factor = Fit_outputfield[mil_temp]/Cross_sec_input[mil_temp]
	
	Beam_parameter("Field_out")
	
	WaveStats/Q/R=(30,94) Cross_sec_phase
	temp_phase = V_avg
		
	Make/O/C/N = (Num_point,Num_point)/D Field_out_defor
	Make/O/C/N = (Num_point,Num_point)/D Field_out_Fit
	
	//Field_out_Fit = Field_in * sqrt(Power_field_out_factor)
	Field_out_Fit = Field_in * Power_field_out_factor
	Field_out_defor = Field_out_Fit - (Field_out * exp(-sqrt(-1) * temp_phase))
	
	Conv_temp1 = 0
	Conv_temp2 = 0
	Conv_temp3 = 0
	
	i = 0
	Do
		j = 0
		Do
	//		radius_temp_sq = X_dist^2 + Y_dist ^2
	//		dist_temp = amplitude_E * Exp(-radius_temp_sq/Waist^2) * Exp(- sqrt(-1) *  k_prop * radius_temp_sq/(2 * Radius))
	//		Field_Start[i][j] = dist_temp
			j = j + 1	
		While ( j < Num_point)			 
		i = i+ 1	
	While ( i < Num_point)	
	
	i = 0
	Do
		j = 0
		Do			
			Conv_temp1 = Conv_temp1 + Field_in[i][j] * conj(Field_in[i][j]) * step_grid^2
			Conv_temp2 = Conv_temp2 + Field_out[i][j] * conj(Field_out[i][j]) * step_grid^2
			Conv_temp3 = Conv_temp3 + (Field_out_Fit[i][j]) * conj(Field_out_Fit[i][j]) * step_grid^2
			j = j + 1	
		While ( j < Num_point)			 
		i = i+ 1	
	While ( i < Num_point)
	
	printf "Power of the input TEM00 %g and the output power %g\r", Conv_temp1,Conv_temp2
	printf "Output power in the mode TEM00 %g \r", (Conv_temp3)

	//printf"Power scattered out of the input mode TEM00:"
	//Calculate_power("Field_out_defor")
	
End

Macro Find_resonance_length()
	
	Print Secs2Time(DateTime,0)
	
	Variable i,j
	Variable iter_max_field, Num_step_length,First_peak,Peak_pos
	
	PauseUpdate; Silent 1
	iter_max_field = 200
	Num_step_length = 40

	createEField()

	Make/O/C/N = (Num_point,Num_point)/D Field_Cav,Field_circ,Field_in,Field_out,Field_temp
	Make/O/N = (Num_point,Num_point)/D Field_dephas
	
	Make/O/N = (Num_step_length)/D Field_circ_temp
	Make/O/N = (Num_step_length)/D Length_temp
		
	Field_in = Field_Start
	Beam_parameter("Field_in")
	Create_mirror()	

	j = 0
	Do
		Length_temp[j] = j * Lambda/(Num_step_length * 2 )
		//Length_temp[j] = j * Lambda/(Num_step_length)
		
		Field_Cav = 0

		Field_temp = Field_in
		Propa_mirror("Field_temp","Mirror_ITM_trans",1)
		Field_circ = Field_temp
	
		i = 0
		Do
			Field_Cav =  Field_Cav + Field_circ 
			Make_propagation("Field_circ",Length_temp[j]+Length_cav)
			Propa_mirror("Field_circ","Mirror_ETM_cav",-1)
			Make_propagation("Field_circ",Length_temp[j]+Length_cav)
			Propa_mirror("Field_circ","Mirror_ITM_cav",1)		
			i = i+ 1	
		While ( i <iter_max_field )
		Calculate_power("Field_Cav")
		Field_circ_temp[j] = Return_power
		j = j + 1
		print j
	While ( j < Num_step_length)
	
	WaveStats/Q Field_circ_temp
	
	FindPeak/M=1 Field_circ_temp
	First_peak = V_PeakLoc
	Peak_pos = Length_temp[First_peak]
	Print Peak_pos
	
	Length_temp = 0
	
	iter_max_field = 800
	Num_step_length = 30
	
	
	j = 0
	Do		
		Length_temp[j] = (Peak_pos + ( j * (4e-9/ Num_step_length) -2e-9 ))
		Field_Cav = 0

		Field_temp = Field_in
		Propa_mirror("Field_temp","Mirror_ITM_trans",1)
		Field_circ = Field_temp
	
		i = 0
		Do
			Field_Cav =  Field_Cav + Field_circ 
			Make_propagation("Field_circ",Length_temp[j]+Length_cav)
			Propa_mirror("Field_circ","Mirror_ETM_cav",-1)
			Make_propagation("Field_circ",Length_temp[j]+Length_cav)
			Propa_mirror("Field_circ","Mirror_ITM_cav",1)		
			i = i+ 1	
		While ( i <iter_max_field )
		Calculate_power("Field_Cav")
		Field_circ_temp[j] = Return_power
		j = j + 1
		print j
	While ( j < Num_step_length)

	Print Secs2Time(DateTime,0)

End

Macro Find_resonance_length_Fine(DeltaL)
	Variable DeltaL
	
	Print Secs2Time(DateTime,0)
	
	Variable i,j
	Variable iter_max_field, Num_step_length,First_peak,Peak_pos
	
	PauseUpdate; Silent 1
	iter_max_field = 200
	Num_step_length = 20

	createEField()

	Make/O/C/N = (Num_point,Num_point)/D Field_Cav,Field_circ,Field_in,Field_out,Field_temp
	Make/O/N = (Num_point,Num_point)/D Field_dephas
	
	Make/O/N = (Num_step_length)/D Field_circ_temp
	Make/O/N = (Num_step_length)/D Length_temp
		
	Field_in = Field_Start
	Beam_parameter("Field_in")
	Create_mirror()	

	j = 0
	Do
		Length_temp[j] = (DeltaL + ( j * (1e-9/ Num_step_length) -0.5e-9 ))
		
		Field_Cav = 0

		Field_temp = Field_in
		Propa_mirror("Field_temp","Mirror_ITM_trans",1)
		Field_circ = Field_temp
	
		i = 0
		Do
			Field_Cav =  Field_Cav + Field_circ 
			Make_propagation("Field_circ",Length_temp[j]+Length_cav)
			Propa_mirror("Field_circ","Mirror_ETM_cav",-1)
			Make_propagation("Field_circ",Length_temp[j]+Length_cav)
			Propa_mirror("Field_circ","Mirror_ITM_cav",1)		
			i = i+ 1	
		While ( i <iter_max_field )
		Calculate_power("Field_Cav")
		Field_circ_temp[j] = Return_power
		j = j + 1
		print j
	While ( j < Num_step_length)
	
	WaveStats/Q Field_circ_temp
	
	FindPeak/M=20 Field_circ_temp
	First_peak = V_PeakLoc
	Peak_pos = Length_temp[First_peak]
	Print Peak_pos
	
	Length_temp = 0

	Print Secs2Time(DateTime,0)

End


Macro Test2()
	Variable i,j
	Variable/C comp
	
	Make/O/N = (Num_point,Num_point)/D Field_diff
	
	PauseUpdate; Silent 1

	i = 0
	Do
		j = 0
		Do			
			// Field_diff[i][j] = imag(r2polar(Field_out[i][j]))
			j = j + 1	
		While ( j < Num_point)			 
		i = i+ 1	
	While ( i < Num_point)

	 comp = (1+sqrt(-1))

	print magsqr(comp)
	
	print (comp * conj(comp))

End

Macro Test5()
	Variable i,j
	Variable X_dist,Y_dist,radius_temp
	
		PauseUpdate; Silent 1
	
	Make/O/N = (Num_point+1)/D Grid_FFT
	i = 0
	Do
		Grid_FFT[i] = -Length_Grid/2+ i*step_grid
		i = i+ 1	
	While ( i < Num_point+1)
	
      Make/O/N = (Num_point,Num_point)/D Mirror_grid
	
	Mirror_grid = 0
	
	i = 0
	Do
		j = 0
		Do
			X_dist = - Length_Grid/2 + step_grid/2 + i * step_grid
		      Y_dist = - Length_Grid/2 + step_grid/2 + j * step_grid
			
			radius_temp = sqrt(X_dist^2 + Y_dist ^2)
			
			if ( radius_temp <0.05)
				Mirror_grid[i][j] = 1
			endif
	
			j = j + 1	
		While ( j < Num_point)			 
		i = i+ 1	
	While ( i < Num_point)		
	
end