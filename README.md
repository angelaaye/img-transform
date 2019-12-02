# img-transform
This repository uses perspective transformation to replace a portion of an existing image with an alternate image. 

## Files
`dlt_homography.m` computes the perspective homography between two images, given four point correspondences.<br/> 
`bilinear_interp.m` performs bilinear interpolation to produce a pixel intenxity value, given an image and a subpixel location (point)<br/>
`histogram_eq.m` performs discrete historgram equilization on the input image (which is 8-bit and greyscale)<br/>
`billboard_hack.m` uses the above functions to create the composite image, which is stored in colour<br/>
