# Radial Profiles of Specific Star Formation in Galaxy Pairs with MaNGA

## The Goal

The goal of this project is to see whether or not galaxy interactions can enhance their rates of star formation. We will compute the specific star formation rate (sSFR) across the surfaces of the galaxies in our sample. Then we will compare the sSFRs in paired galaxies to a sample of control galaxies to test if galaxy interactions do indeed enhance the sSFR. 

## Radial Profiles
MaNGA provides a 2D distribution of galaxy spectra for each of its observations. To prep the data for comparison, we will reduce the 2D distribution of ssfr in each observation to a 1D ssfr profile as a function of the galaxy's radius. The radius is given in terms of the effective radius, the radius which contains 50% of the total light for the galaxy, so that the profiles of galaxies of various sizes can be compared with each other. Typical galaxies are circular; however, if they are at an angle with the observer they will appear to be an ellipse on the sky. We deproject the geometry of the galaxies to account for this effect.

## Sample Selection
Next we will define our pair and control samples. While two galaxies may appear close to each other on the sky, they are not necessarily close to each other in space. To verify that they are paired galaxies, the two galaxies need to be near each other on the sky, have a similar redshift (which tells us how far the galaxies are from us), and have a limited relative velocity between the two galaxies. Paired galaxies are identified within the MaNGA survey's footprint and from other surveys which share MaNGA's footprint. The control galaxies are then the MaNGA galaxies with no observed companion galaxy.

## Analysis
With the two sample constructed, we compare the sSFR data between the pair and control samples.
![](https://github.com/jlsteffen/MergerSF/blob/main/ssfr_comb.pdf)

This work has been published to Steffen et al. 2021, https://arxiv.org/pdf/2102.03398.pdf.
