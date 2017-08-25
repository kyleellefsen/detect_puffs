# Guide for setting parameters in flika's detect_puffs plugin

There are quite a few parameters to set when running the threshold_cluster function in the detect_puffs plugin. The function GUI has brief 
explanations for these parameters, and this guide attempts to extend those explanations. This guide assumes that you have already read the 
detect_puffs [tutorial](http://htmlpreview.github.io/?https://github.com/kyleellefsen/detect_puffs/blob/master/docs/How%20to%20detect%20subcellular%20Ca2%2B%20signals%20using%20Flika.html).

## Parameters

### data_window (Window)
This window is the F/F<sub>0</sub> window.

### normalized_window (Window)
This Window should have the baseline at 0. The standard deviation of the noise should be equal to 1. 
This window is used for gaussian fitting and kinetics.

### blurred_window (Window)
The (spatially) gaussian blurred version of the normalized_window. The movie in this window should be blurred in proportion to how 
spatially diffuse the signals you are trying to detect are. This window will be used for event detection.

### roi_width (int)
The width of the ROI in pixels. This parameter cannot be adjusted after it is passed into the `threshold_cluster()` function.

### paddingXY (int)
How many pixels do you want to pad the ROI by when fitting with gaussian. If the `roi_width` is set to 3 and `paddingXY` is set to 10,
then the region the gaussian fitting will be performed on will be 10+3+10 = 23 pixels wide. 

### paddingT_pre (int)
How many frames before the event detection should we look for the start of the event.

### paddingT_post (int)
How many frames after the event detection should we look for the end of the event.

### maxSigmaForGaussianFit (int)
When fitting with a gaussian, what is the upper bound for the sigma (width) or sigmax and sigmay (if fitting with a rotatable gaussian) parameter.

### rotatedfit (bool)
Set this to True to fit to a 2D rotatable gaussian. If False, each event will be fit with a single circularly symmetric sigma (width) parameter. A rotatable gaussian has additional fitting parameters: a long axis with a sigmax (width), a short axis with a sigmay, and a theta (angle).

### radius (float)
Puffs seperated by less than this amount (in pixels) will be grouped together in a group. If two puff centroids are 5 pixels apart and 
the radius is set to 6, these two puff centroids will be grouped together. If you have three puff centroids (A, B, and C) in a line and 
the distance between A and B is 5, the distance between B and C is 5, the distance between A and C is 10, all three will be grouped
together. This is because A and B will be grouped together and B and C will be grouped together. Because of this, depending on the 
geometry and density of your events, events which are seperated by much greater than `radius` might be grouped together.

The `radius` parameter cannot be changed after it is passed into the `threshold_cluster()` function. However, the user can manually 
regroup events into clusters. Take the ROI, move it over a cluster, and press down the 'u' key on your keyboard to ungroup all events
in that cluster, so that each event is the only event in its own cluster. To group clusters together, move the ROI so that several 
clusters are within the ROI, and press down the 'g' key. 

### maxPuffLen (int)
maximum possible duration of puff length in frames.  If puffs last longer, use a bigger value.  This affects 'chunking' when 
determining distances between pixels.  It's better to overestimate, although increasing too much will slow down analysis.

### maxPuffDiameter (int)
This affects 'chunking' when determining distances between pixels.    This value divided by two is also used as the radius for 
calculating densities

### blur_thresh (float)
Get rid of all the pixels you are sure don't belong to puffs. All pixels in the blurred_window movie below this threshold won't be considered to be part of an event.

### time_factor (float)
When looking for the nearest point of higher density, time and space don't need to be weighted equally.  This is how much time will be 
'shrunk' in order to calculate distances relative to pixels.  Practically, the larger this number, the more points separated in time 
will tend to be clustered together.  Use a higher time factor if you are sampling at high frame rates.

### load_flika_file (bool)
If this is set to True and there is a flika file with the same name and in the same location as the data_window file, all the processing
will be skipped and the flika file will be loaded.


## FAQs

### Q) Which movie are the ΔF/F<sub>0</sub> values calculated from?

A) They are taken from the data_window. The value for each frame is the average of the pixel intensities inside the ROI at that frame. 
The width of the ROI is set by the the `roi_width` parameter, and the center of the ROI is the event centroid as determined during the 
gaussian fitting step.

### Q) Is the radius parameter the distance from one cluster center to another?

A) No. The radius is used to determine if two events are part of the same cluster. Inter cluster distances will all end up being greater
than the value of the `radius` parameter, because any two individual events that are within that radius will be grouped together inside 
a cluster.

### Q) In the puff analyser, I sometimes get peaks identified as puffs that don’t really look as “proper” puffs in the ΔF vs time graph. How can I best locate this specific puff in my image stacks to make sure that what I detect is a real event?

A) When a puff is selected in the puff analyzer GUI, the ROI in the data_window will move to be centered around that puff's centroid. 
The centroid will also be highlighted. If you link the data_window to the blurred_window (by right-clicking the timeline at the bottom
of the blurred_window), and copy and paste the ROI into the blurred_window, you can hover over different times in the F/F<sub>0</sub>
trace window to visualize the event. This makes visual determination of if a detected event is a true positive much easier.

### Q) Sometimes I detect puffs which take very long to terminate. The problem I have with them is that in the ΔF vs time graph I cannot drag the dotted lines further left or right to capture the whole duration of the puff. Is there a way around this?

A) Yes there is! You can change the x-axis in the puff analyzer's ΔF vs time graph by adjusting the `paddingT_pre` and `paddingT_post`
variables. Unfortunately this requires rerunning the `threshold_cluster()` function.
