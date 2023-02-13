# Radial Profiles of Specific Star Formation in Galaxy Pairs with MaNGA

## The Goal

The goal of this project is to see whether or not galaxy interactions can enhance their rates of star formation. We will compute the specific star formation rate (sSFR) across the surfaces of the galaxies in our sample. Then we will compare the sSFRs in paired galaxies to a sample of control galaxies to test if galaxy interactions do indeed enhance the sSFR. 

## The Survey 

This work is based on the data from the SDSS's (Sloan Digital Sky Survey) MaNGA (Mapping Nearby Galaxies at Apache Point Observatory) survey. MaNGA is an integral field spectroscopic (IFS) survey which has observed 10,000 nearby galaxies. Spectroscopy takes the light from an object and splits it by wavelength. We can tell many things about an object by its spectra. You can tell what king of object you are looking at, the ages of stars, the chemical composition, et cetera.

Traditionally, for large spectrscopic surveys, light is collected through fiber optic cables to prevent spectra from overlapping with other light sources. A single fiber optic cable is placed on a single source, and a single spectrum is collected. IFS takes a whole bundle of fiber optic cables and places them onto a single target. This gives astronomers several spectra for a single target, each of which covers a different part of the object. With MaNGA, IFS allows us to simultaneously study a galaxy's center and disk with a single observation.

![Figure 1](https://github.com/jlsteffen/Star-Formation-in-Mergers/blob/main/images/manga_v3.jpg)

## MaNGA Object Catalog

## Sample Selection
Next we will define our pair and control samples. While two galaxies may appear close to each other on the sky, they are not necessarily close to each other in space. To verify that they are paired galaxies, the two galaxies need to be near each other on the sky, have a similar redshift (which tells us how far the galaxies are from us), and have a limited relative velocity between the two galaxies. Paired galaxies are identified within the MaNGA survey's footprint and from other surveys which share MaNGA's footprint. The control galaxies are then the MaNGA galaxies with no observed companion galaxy.

## Radial Profiles
MaNGA provides a 2D distribution of galaxy spectra for each of its observations. To prep the data for comparison, we will reduce the 2D distribution of ssfr in each observation to a 1D ssfr profile as a function of the galaxy's radius. The radius is given in terms of the effective radius, the radius which contains 50% of the total light for the galaxy, so that the profiles of galaxies of various sizes can be compared with each other. Typical galaxies are circular; however, if they are at an angle with the observer they will appear to be an ellipse on the sky. We deproject the geometry of the galaxies to account for this effect.

## Analysis
With the two samples constructed and radiak profiles created, we compare the sSFR data between the pair and control samples.

![Figure 2](https://github.com/jlsteffen/Star-Formation-in-Mergers/blob/main/images/ssfr_comb.jpg)

We see that the paired galaxies have higher rates of sSFR in their centers than the control galaxies but the two samples have similar rates of sSFR in their disks.  

This work has been published to Steffen et al. 2021, https://arxiv.org/pdf/2102.03398.pdf.
