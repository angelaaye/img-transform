function [H, A] = dlt_homography(I1pts, I2pts)
% DLT_HOMOGRAPHY Perspective Homography between two images.
%
%   Given 4 points from 2 separate images, compute the perspective homography
%   (warp) between these points using the DLT algorithm.
%
%   Inputs:
%   -------
%    I1pts  - 2x4 array of points from Image 1 (each column is x, y).
%    I2pts  - 2x4 array of points from Image 2 (1-to-1 correspondence).
%
%   Outputs:
%   --------
%    H  - 3x3 perspective homography (matrix map) between image coordinates.
%    A  - 8x9 DLT matrix used to determine homography.

%--- FILL ME IN ---

% Insert code here.

%------------------

% Step 1: Compute similarity transform for I1pts

% take the four points, compute their center
x_val = sum(I1pts(1, :))/4;
y_val = sum(I1pts(2, :))/4;
center = [x_val; y_val];

% compute the total distance from all four points to the origin
% divide each vector by the current average distance (to get unit average)
% multiply by sqrt(2)
tot_dist = 0;
for i = 1:4
    tot_dist = tot_dist + sqrt((I1pts(1, i)-x_val)^2 + (I1pts(2, i)-y_val)^2);
end
avg_dist = tot_dist/4;

% obtain points after scaling, sqrt(2) average distance 
I1pts_scaled = I1pts*sqrt(2)/avg_dist;
%compute new center
x_val = sum(I1pts_scaled(1, :))/4;
y_val = sum(I1pts_scaled(2, :))/4;
center = [x_val; y_val];

% find similarity transform T [sR|t], s = scale factor, t = displacement
T1 = [sqrt(2)/avg_dist 0 -x_val; 0 sqrt(2)/avg_dist -y_val; 0 0 1];
I1pts = [I1pts; 1 1 1 1]; %homogenous points
P1 = T1*I1pts;

% Step 2: Do the same for I2pts
x_val = sum(I2pts(1, :))/4;
y_val = sum(I2pts(2, :))/4;
center = [x_val; y_val];

tot_dist = 0;
for i = 1:4
    tot_dist = tot_dist + sqrt((I2pts(1, i)-x_val)^2 + (I2pts(2, i)-y_val)^2);
end
avg_dist = tot_dist/4;

I2pts_scaled = I2pts*sqrt(2)/avg_dist;
x_val = sum(I2pts_scaled(1, :))/4;
y_val = sum(I2pts_scaled(2, :))/4;
center = [x_val; y_val];

T2 = [sqrt(2)/avg_dist 0 -x_val; 0 sqrt(2)/avg_dist -y_val; 0 0 1];
I2pts = [I2pts; 1 1 1 1];
P2 = T2*I2pts;

% Step 3: Find homography matrix
A = [];
for i = 1:4
    A = [A; -P1(1, i) -P1(2, i) -1 0 0 0 P2(1, i)*P1(1, i) P2(1, i)*P1(2, i) P2(1, i); 0 0 0 -P1(1, i) -P1(2, i) -1 P2(2, i)*P1(1, i) P2(2, i)*P1(2, i) P2(2, i)];
end
H_col = null(A);
H_mtx = transpose(reshape(H_col, 3, 3));
H = inv(T2)*H_mtx*T1;


end
