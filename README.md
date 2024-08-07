# Automation problems in otolith image analysis
     
## Strategy 
1. Determine the **absolute scale** of an image using the ruler contained in the photo.

2. **Contour of the otolith**:
   - Identify the **set of points** which make up the contour.
   - **Model the contour** points (e.g. using splines) 
   - **Identify reference axes** and orient the otolith contour and image a standardized way.
   - Show how **contours vary with age** using overlayed contours.

3. Determine the **position of the nucleus** of an otolith. Helpful clues:
   - If we consider the image as an intensity map, image gradients will generally point towards or away from the nucleus.
   - From the position of the nucleus, radial scans of its surroundings are more similar than those of other positions.

4. Using annotated otolith images, find a **non-linear transform** that remaps the ring markers so that they are equally spaced, on average.
   - Explore how this non-linear transform varies from otolith to otolith.

5. Using a set of otolith images that have been properly rotated and resized, find what an **“average” otolith** looks like. 

6. **Ring identification**:
   - Using the otolith contour shape as a starting point, vary its size (and shape?) such that it maps to all rings contained in the image (i.e. model how the otolith changed size and shape during its growth).
   - How do the average intensities of the rings vary with the size of the rings (i.e. this is a way of getting a standardized radial intensity profile).
   - Use radial scans from the nucleus to the edge of the otolith and check if their averages remove unwanted noise (i.e. boost the signals from “true” rings).
   - How the “true” rings identified in the annotated otolith images compare with the observed rings in images. What criteria are used to dismiss “false” rings?

## Analytical problems 

### Similarity measure between two warped vectors 

Say we have two image intensity profiles for two lines starting at the nucleus to the edge of an otolith. Even if we standardize for distance, their features (e.g. annuli) do not generally line up. 
The problem is this : can we create a univariate function that warps the distance in such a way that makes the features of one profile line up with the other? The matching between two vectors proceeds on two fronts: a warping of the coordinates and tranforming the intensities (signal) from one vector so that it is on the same scale as the other. 
 - Cumulative density function of a beta distribution provides 

### Squaring the Otolith
Look up image intensties along a set of rays that goes from the nucleus to edge of the otolith, at regular intervals around the otolith. This produces an unwrapped otolith image. The goal at this point is this:
- Transform the intensities along each ray such that the otolith rings are aligned and the intensities are normalized to some reference vector.
- Stretch out the squared otolith in such a way that its rings are more or less evenly spaced.

### Modelling intensities
Rings are periodic events, model the variations using a repeating element, e.g. a periodic spline. 







