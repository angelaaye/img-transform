function [b] = bilinear_interp(I, pt)
% BILINEAR_INTERP Performs bilinear interpolation for a given image point.
%
%   Given the (x, y) location of a point in an input image, use the surrounding
%   4 pixels to conmpute the bilinearly-interpolated output pixel intensity.
%
%   Note that images are (usually) integer-valued functions (in 2D), therefore
%   the intensity value you return must be an integer (use round()).
%
%   This function is for a *single* image band only - for RGB images, you will 
%   need to call the function once for each colour channel.
%
%  Inputs:
%  -------
%   I   - Single-band (greyscale) intensity image, 8-bit (i.e., uint8).
%   pt  - 2x1 point in input image (x, y), with subpixel precision.
%
%  Outputs
%  -------
%   b  - Interpolated brightness or intensity value (whole number >= 0).

%--- FILL ME IN ---

% Insert code here.

%------------------

% Get neighbouring 4 points

x = pt(2);
y = pt(1);

x1 = floor(x);
x2 = floor(x)+1;
y1 = floor(y);
y2 = floor(y)+1;

if (x1 == x) && (y1 == y)
    b = round(double(I(x, y)));
else
    % Bilinear interpolation (formula from Wiki)
    b = 1/((x2-x1)*(y2-y1))*[x2-x x-x1]*[double(I(x1, y1)) double(I(x1, y2)); double(I(x2, y1)) double(I(x2, y2))]*[y2-y; y-y1];
    b = round(b);
end

end