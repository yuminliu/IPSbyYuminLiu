function PPG = graygrad(f, T)
%COLORGRAD Computes the vector gradient of an RGB image.
%   PPG = COLORGRAD(F, T) computes PPG, the per-plane composite gradient
%   obtained by summing the 2-D gradients of the individual color 
%   planes. Input T is a threshold in the range [0, 1]. If it is
%   included in the argument list, the values of PPG are
%   thresholded by letting PPG(x,y) = 0 for values <= T and PPGG(x,y) =
%   PPG(x,y) otherwise.If T is not
%   included in the argument list then T is set to 0. Both output
%   gradients are scaled to the range [0, 1].

%   Copyright 2002-2015 Yumin Liu, R. C. Gonzalez, R. E. Woods, and S. L. Eddins
%   From the book Digital Image Processing Using MATLAB, 2nd ed.,
%   Gatesmark Publishing, 2009.
%
%   Book web site: http://www.imageprocessingplace.com
%   Publisher web site: http://www.gatesmark.com/DIPUM2e.htm

% Compute the x and y derivatives of the three component images 
% using Sobel operators.
if (~ismatrix(f))
   error('Input image must be gray image.');
end

sh = fspecial('sobel');
sv = sh';
Gx = imfilter(double(f(:, :, 1)), sh, 'replicate');
Gy = imfilter(double(f(:, :, 1)), sv, 'replicate');

% Compute the per-plane gradients.
G = sqrt(Gx.^2 + Gy.^2);

% Form the composite by adding the individual results and
% scale to [0, 1].
PPG = mat2gray(G);

% Threshold the result.
if nargin == 2
   PPG = (PPG > T).*PPG;
end
