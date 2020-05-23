function I1 = Do_Virtual_map(G,Power_law)
%Do_Virtual_map() Create a synthetic map according to a parametrised PSD
%   Grid is the Grid object for the calculations
%   Power_law is a vector of length 8 with the power law, see the function
%   DiagonalPSDLaw2D_LMA2 for the definition
%   I1 is the created virtual map.

p = inputParser;
p.FunctionName = 'Create a virtual map';

% Check if the first argument is a grid
addRequired(p,'G',@(x)isa(x, 'Grid'));

% Check if the PSD law is given
addRequired(p,'Power_law',@(x)isnumeric(x))

p.parse(G,Power_law)

I1 = FakeMirror2D(G.Step,G.Num_point,Power_law);

%% Here the great original functions from F. Bondu with only light modifications/simplifications or comments


    function fit_spec = DiagonalPSDLaw2D_LMA2(freq,law)
        % spec=DiagonalPSDLaw2D(f)
        % linear spectral density dummy spectrum
        %
        % in m^2/m^{-2}
        % single sideband
        %
        
        % initiate the matrix for the right size
        fit_spec = freq;
        
        nb_param = length(law);
        
        switch nb_param
            
            case 2
                amp = law(1);
                expon = law(2);
                
                amp = amp*1E-16;
                
                fit_spec = amp * freq.^(expon-1);
                fit_spec(freq==0) = 0;
            
            case 5
                amp = law(1:2);
                expon = law(3:4);
                
                amp(1) = amp(1)*1E-16;
                amp(2) = amp(2)*1E-14;
               
                f1 = law(5);
                
                fi = find( freq <= f1 );
                fit_spec(fi) = amp(1) * freq(fi).^(expon(1)-1);
                              
                fi = find( freq > f1);
                fit_spec(fi) = amp(2) * freq(fi).^(expon(2)-1);
                
                fit_spec(freq==0) = 0;
                              
            case 8  
                amp = law(1:3);
                expon = law(4:6);
                
                amp(1) = amp(1)*1E-16;
                amp(2) = amp(2)*1E-23;
                amp(3) = amp(3)*1E-12;
                
                f1 = law(7);
                f2 = law(8);
                
                fi = find( freq <= f1 );
                fit_spec(fi) = amp(1) * freq(fi).^(expon(1)-1);
                
                fi = find( freq <= f2 & freq > f1 );
                fit_spec(fi) = amp(2) * freq(fi).^(expon(2)-1);
                
                fi = find( freq > f2 );
                fit_spec(fi) = amp(3) * freq(fi).^(expon(3)-1);
                
                fit_spec(freq==0) = 0;
            
            otherwise
                %nb_param
                error('DiagonalPSDLaw2D_LMA2(freq,law), the ''law'' variable must be a vector with 2, 5 or 8 elements')
                
        end
        %figure(300);loglog(freq(1,2:end),fit_spec(1,2:end))     
    end

    function mirgrid=FakeMirror2D(deltax,N,law)
        % mirgrid=FakeMirror2D(deltax,N,law);
        % input: deltax: grid step size
        %        N: number of points horizontally and vertically - even !!
        %        law: the parametrised 1D PSD
        % output: a square grid
        
        if (mod(N,2)==1)
            disp('FakeMirror2D(), number of point of the grid must be even' )
            mirgrid = 0;
            return
        end
        
        deltaf = 1/(N*deltax);
        
        % first define in frequency domain on a half plane
        
        % two squares as bricks to define a half space
        matU = SquareDiagonalLSD(deltaf,N/2,law);
        matV = SquareDiagonalLSD(deltaf,N/2,law);
        
        %figure(301); imagesc(abs(matU)); axis square
        
        % first quadrant : simplest one
        matA = matU;
        
        % fourth quadrant
        matD = conj(matU);
        matD = rot90(matD,2);
        matD(2:end,:) = matD(1:end-1,:);
        matD(:,2:end) = matD(:,1:end-1);
        matD(:,1) = 0; % Nyquist frequency
        matD(1,:) = 0; % Nyquist frequency
        
        % second quadrant: remove column with DC,
        %                add line with zero at Nyquist frequency
        %matB = matV;
        matB = fliplr(matV);
        matB(:,2:end) = matB(:,1:end-1);
        matB(:,1) = 0;
        matB(1,2:end) = conj(fliplr(matA(1,2:end)));
        
        % third quadrant:
        matC = flipud(conj(matV));
        matC(2:end,:) = matC(1:end-1,:);
        matC(1,:) = 0;
        matC(2:end,1) = conj(flipud(matA(2:end,1)));
        
        Freq_mirgrid = [matA matB ; matC matD];
        
        %mirgrid = fftshift(Freq_mirgrid);
        %figure(302); imagesc(abs(Freq_mirgrid)); axis square
        
        
        mirgrid = ifft2(Freq_mirgrid,'symmetric');
        
        %figure(303); imagesc(abs(mirgrid)); axis square
        
        calib = 1/(sqrt(2)*deltax/N);
        mirgrid = real(mirgrid*calib);
        
    end

    function tab = SquareDiagonalLSD(deltaf,N,law)
        % generates a square of random numbers
        % tab(1,1) = mean = 0
        % tab(u,v) is an imaginary number
        % |tab(u,v)|^2 is PSD(f) with f^2=u^2+v^2
        % uses DiagonalPowerSpectralDensityLaw2D: definition of radial PSD
        
        % x, y and frequency grids;
        [xgrid,ygrid] = meshgrid(0:deltaf:deltaf*(N-1),0:deltaf:deltaf*(N-1));
        fgrid = sqrt(xgrid.^2 + ygrid.^2);
        
        % tabs of gaussian random numbers
        [tabx,taby] = gaussian(N,0,1/sqrt(2));
        
        lsd = sqrt(DiagonalPSDLaw2D_LMA2(fgrid,law));
        
        %figure(304); imagesc(fgrid); axis square
        
        
        
        sqx = tabx.*lsd;
        sqy = taby.*lsd;
        
        tab = sqx+1i*sqy;
        
    end

    function [y1,y2]=gaussian(N,moy,sigma)
        %[y1,y2]=gaussian(N,moy,sigma)
        % input:  N, mean, standard deviation
        % output: 2 N*N tables of random numbers
        
        x1 = rand(N); % table of N*N random numbers
        x2 = rand(N); % an other one
        
        y1 = sqrt(-2*log(x1)).*cos(2*pi*x2)*sigma+moy;
        y2 = sqrt(-2*log(x1)).*sin(2*pi*x2)*sigma+moy;
        
    end




end


