Current Version: Release Candidate 1.6.3
email phi3@umbc.edu  with any questions or concerns

How to Use:

Take a picture of the spider web.
	All strands should be in focus.
	The web should be within a single plane orthogonal to the camera's angle. (future versions may be able to operate at any angle)

Change the image_name to the picture of the spider web.
	

Adjust Settings. (1-5 should be near optimized already; this is mainly for 6, 7, and 8)
	1. base_distance = the base radius of pixels counts a point as being near another
	2. scaling_distance = increases #1's radius by an amount that scales with the image or ROI size
	3. base_radius = feeds into subclust
	4. scaling_radius = the 
	5. percent_points = the percentage of points kept as deemed most accurate enough to feed into subtractive clustering
	6. scan_width = one dimension of the plane that the web sits on
	7. scan_length = the other dimension of the plane that the web sits on
	8. z_axis = the distance the web rests from the camera

Run script.
	Select a rectangular region of interest. 
	Try to avoid obstructions.
	(future versions may include the option of a polygonal region of interest=- NOTE: polygonal ROI will likely slow down runs significantly)

Go to the directory where the script is held.

Get "results.txt" and import it into Polytec's software.

