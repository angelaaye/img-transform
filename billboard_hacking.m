% Billboard hack script file.

% Bounding box in Y & D Square image.
bbox = [404, 490, 404, 490;
         38,  38, 354, 354];

% Point correspondences.
Iyd_pts = [416, 485, 488, 410;
            40,  61, 353, 349];
Ist_pts = [2, 219, 219,   2; 
           2,   2, 410, 410];

Iyd = imread('yonge_dundas_square.jpg');
Ist = imread('uoft_soldiers_tower_dark.png');

% Make a copy for hacking... your script MUST use the Ihack variable name.
Ihack = Iyd; 


% Let's do the histogram equalization first...
Jst = histogram_eq(Ist);

% Compute the perspective homography we need...
[H, ~] = dlt_homography(Iyd_pts, Ist_pts);
%new_pts = H*[Iyd_pts; 1 1 1 1];
%new_pts = new_pts./new_pts(3,:);

% Main 'for' loop to do the warp and insertion - 
% this could be vectorized to be faster if needed.
for i = min(Iyd_pts(1, :)):max(Iyd_pts(1, :))
    for j = min(Iyd_pts(2, :)):max(Iyd_pts(2, :))
        if (inpolygon(i, j, Iyd_pts(1, :), Iyd_pts(2, :)) == 1)
            new_pt = H*[i; j; 1];
            new_pt = new_pt/new_pt(3);
            intensity = bilinear_interp(Jst, [new_pt(1), new_pt(2)]);
            Ihack(j, i, :) = uint8(intensity);
        end
    end
end

% You may wish to make use of the inpolygon() function (see MATLAB help).
%Ihack = Ihack(404:490, 38:354);
figure;  imshow(Ihack);

%---------- Functions Go Below ----------

% ADD YOUR DLT_HOMOGRAPHY CODE BELOW
function [H, A] = dlt_homography(I1pts, I2pts)

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

% ADD YOUR BILINEAR_INTERP CODE BELOW
function [b] = bilinear_interp(I, pt)

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

% ADD YOUR HISTOGRAM_EQ CODE BELOW
function [J] = histogram_eq(I)

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
