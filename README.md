# Radial Profiles of Specific Star Formation in Galaxy Pairs with MaNGA

The goal of this project is to see whether or not galaxy mergers incite new, merger-induced, star foramtion in paired galaxies. I will compute the radial profiles of the specific star formation rate (ssfr) in galaxy pairs and isolated control galaxies. Then I will take the difference between the profiles of the paired galaxies and the isolated control galaxies to see if paired galaxies feature higher or lower ssfr in comparison to control galaxies.

## Radial Profiles
MaNGA provides a 2D distribution of galaxy spectra for each of its observations. To prep the data for comparison, I will reduce the 2D distribution of ssfr in each observation to a 1D ssfr profile as a function of the galaxy's radius. The radius is given in terms of the effective radius, the radius which contains 50% of the total light for the galaxy, so that the profiles of galaxies of various sizes can be compared with each other. The galaxy profiles will also be deprojected to a circular profile using each galaxy's position angle. 

## Sample Selection
Next I will define each of our pair and control samples. Paired galaxies are selected to be within a line-of-sight velocity of 500 km/s, a projected separation of 50 kpc, and a stellar mass range of log(M/M_sun) = 9-11.5. Pairs are selected from within the MaNGA fields-of-view and from the NASA-Sloan atlas. Control galaxies are the MaNGA galaxies with no observed companion galaxy. Star forming galaxies are selected by BPT analysis.

## Analysis
I then take the median ssfr radial profiles for pairs and control galaxies and compare the profiles against one another.
