function y = FitFunctionPSD(P,xdata)
% Return the idealised PSD 1D function with several segments
% using the set of parameters P on the frequency axis xdata

P(1) = P(1)*1E-16;
if(rem(length(P),2) ~= 0)
    warning('P has to be an even number of parameters')
end

Half = length(P)/2;

Frequency_x = zeros(1,Half);
x = zeros(1,Half);
A = zeros(1,Half);

for j = 1:Half
    x(j) = P(Half+j);
    Frequency_x(j) = P(j);
    if j == 1
        A(j) = P(1);
    else
        A(j) =  A(j-1)*(power(Frequency_x(j),x(j)-x(j-1)));
    end
end

y = zeros(1,length(xdata))';

% ---- First segment ---- %
index_vec = (xdata<Frequency_x(2));
y(index_vec) = A(1)*power(xdata(index_vec),-x(1));

% ---- Intermediate segments ---- %
for j = 2:Half-1
    index_vec = ((xdata>=Frequency_x(j))&(xdata<Frequency_x(j+1)));
    y(index_vec) = A(j)*(power(xdata(index_vec),-x(j)));
end

% ---- Last segment ---- %
index_vec = (xdata>=Frequency_x(end));
y(index_vec) = A(end)*(power(xdata(index_vec),-x(end)));
end