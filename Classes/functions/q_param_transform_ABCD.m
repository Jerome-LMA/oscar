function q_out = q_param_transform_ABCD(Mat_ABCD,q_in)
%q_param_from_ABCD() simple function to change to pass the Gaussian complex
%beam parameter q_in to the ABCD matrix Mat_ABCD to return q_out
%   Detailed explanation goes here

q_out =  (Mat_ABCD(1,1)*q_in + Mat_ABCD(1,2))/(Mat_ABCD(2,1)*q_in + Mat_ABCD(2,2));

end