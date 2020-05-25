function Eout = Transmit_Aperture(Ein,Diam,varargin)
% Transmit_Aperture() Propagate the E_field through an aperture circular or
% squared
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Transmit_Aperture(E_in,diam) pass the field E_in through a circular
%      aperture of diameter diam
%     Transmit_Aperture(E_in,diam,'Shape','square') pass the field E_in through  a
%     square aperture of side length diam
%     Transmit_Aperture(E_in,'Shape','batman') pass the field E_in through  a
%     batman sign, size must be between 0 and 1, 1 being for the sign on
%     the full grid.
%     To make your own image, took a black and white jpeg picture. load it
%     and convert it to mask
%      image = imread('where is my file');
%      image_gray =  image(:,:,1);
%      then save it as a matlab file called : 'sign.mat' and save it under
%      the folder '@E_Field'
%
%      All the aperture are centered on the grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p  = inputParser;

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check if the first argument is an E_Field
p.addRequired('Diam', @(x)isnumeric(x) && x>0);

% Check if the first argument is an E_Field
p.addParameter('Shape','round',@(x)strcmpi(x,'round') | ...
    strcmpi(x,'square') | strcmpi(x,'batman') );

p.parse(Ein,Diam,varargin{:})

Diam_mask = Diam;

if strcmp(p.Results.Shape,'round')
    Mask = zeros(Ein.Grid.Num_point,Ein.Grid.Num_point,'double');
    Mask(Ein.Grid.D2_r < (Diam_mask/2)) = 1;
    
elseif strcmp(p.Results.Shape,'square')
    Mask = zeros(Ein.Grid.Num_point,Ein.Grid.Num_point,'double');
    Mask((abs(Ein.Grid.D2_X) < (Diam_mask/2)) & (abs(Ein.Grid.D2_Y) < (Diam_mask/2))) = 1;
    
elseif strcmp(p.Results.Shape,'batman')
    %Mask = ones(Ein.Grid.Num_point,Ein.Grid.Num_point);
    load('sign.mat')   
    [nlin,ncol] = size(image_gray);
    
    % Add 0 to render image square, in the original image ncol is
    % even and is superior to nlin which is also even
    
    diff_size = ncol - nlin;
    
    image_gray = vertcat(image_gray,zeros(diff_size/2,ncol));
    image_gray = vertcat(zeros(diff_size/2,ncol),image_gray);
    
    %figure(1);imagesc(image_gray); axis square;
    
    % find the 4 corners of the image
    
    [C,~] = max(image_gray);
    caa = find(C >1,1,'first');
    cbb = find(C >1,1,'last');
    
    [C,~] = max(image_gray,[],2);
    raa = find(C >1,1,'first');
    rbb = find(C >1,1,'last');
    
    % dimension of the white image in pixel
    
    clength = cbb - caa;
    rlength  = rbb - raa;
    
    % scaling factor desired image
    
    dfactor = (ncol / clength) * Diam_mask;
    
    % Create the scale for the image mask
    vec_scale = linspace(-(dfactor)*Ein.Grid.Length/2,(dfactor)*Ein.Grid.Length/2,ncol);
    
    [mask_D2_X,mask_D2_Y] = meshgrid(vec_scale);
    
    Mask = interp2(mask_D2_X,mask_D2_Y,image_gray,Ein.Grid.D2_X,Ein.Grid.D2_Y,'nearest',0);
    
    % The mask is from 0 to 255 in amplitude bring it back from 0
    % to 1
    mask_thresold = 10;
    
    Mask(Mask < mask_thresold)  = 0;
    Mask(Mask >= mask_thresold) =1;
    
    Mask = double(Mask);
    
    %figure(2);imagesc(vec_scale,vec_scale,image_gray); axis square;
    %figure(3);imagesc(Mask); axis square;
    
end

Eout =  Ein.*Mask;

end





