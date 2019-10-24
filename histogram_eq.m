function [J] = histogram_eq(I)
% HISTOGRAM_EQ Histogram equalization for greyscale image.
%
%   Perform histogram equalization on the 8-bit greyscale intensity image I
%   to produce a contrast-enhanced image J. Full details of the algorithm are
%   provided in the Szeliski text.
%
%   Inputs:
%   -------
%    I  - Single-band (greyscale) intensity image, 8-bit (i.e., uint8).
%
%   Outputs:
%   --------
%    J  - Contrast-enhanced greyscale intensity image, 8-bit (i.e., uint8).

%--- FILL ME IN ---

% Insert code here.

%------------------

% Find number of pixels at each grey scale level [0, 255]
h = zeros(1,256);
[h, bins] = histcounts(double(I), (0:1:256));
[m, n] = size(I);
N = m*n; %number of pixels in total

% Compute cumulative distribution function
c = zeros(1,256);
for i = 1:256
    if i == 1
        c(1,i) = h(1,i);
    else
        c(1,i) = c(1,i-1) + h(1,i);
    end
end
    
% Find minimum intensity that has a non-zero value
ind = find(h);
c_min = bins(ind(1));

% Scale intensities to [0, 255]
c_round = zeros(1,256);
for i = 1:256
    c_round(1,i) = round((c(1,i) - c_min)/(N - c_min)*255);
end

% Compute contrast-enhanced greyscale image
for i = 1:m
    for j = 1:n
        intensity = I(i,j);
        J(i,j) = uint8(c_round(intensity+1));
    end
end

end