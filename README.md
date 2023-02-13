# Radial Profiles of Specific Star Formation in Galaxy Pairs with MaNGA

## The Goal
The goal of this project is to see whether or not galaxy interactions can enhance their rates of star formation. We will compute the specific star formation rate (sSFR) across the surfaces of the galaxies in our sample. Then we will compare the sSFRs in paired galaxies to a sample of control galaxies to test if galaxy interactions do indeed enhance the sSFR. 

## The Survey 
<img align="left" width="400" height="350" src="https://github.com/jlsteffen/Star-Formation-in-Mergers/blob/main/images/manga_v3.jpg"> 

This work is based on the data from the SDSS's (Sloan Digital Sky Survey) MaNGA (Mapping Nearby Galaxies at Apache Point Observatory) survey. MaNGA is an integral field spectroscopic (IFS) survey which has observed 10,000 nearby galaxies. Spectroscopy takes the light from an object and splits it by wavelength. We can tell many things about an object by its spectra. You can tell what king of object you are looking at, the ages of stars, the chemical composition, et cetera.

Traditionally, for large spectrscopic surveys, light is collected through fiber optic cables to prevent spectra from overlapping with other light sources. A single fiber optic cable is placed on a single source, and a single spectrum is collected. IFS takes a whole bundle of fiber optic cables and places them onto a single target. This gives astronomers several spectra for a single target, each of which covers a different part of the object. With MaNGA, IFS allows us to simultaneously study a galaxy's center and disk with a single observation.

## MaNGA Object Catalog
The First step of the project is to build a sample of merging galaxies in the MaNGA survey. While the MaNGA survey contains information for its 10,000 target galaxies (~6,400 galaxies at the time of this project), there are many other objects which may fall within the fields-of-views of the survey. Many of these objects will be stars and foreground/background galaxies; however, some of these objects may be galaxies which are gravitationally paired with the MaNGA target galaxy. In a previous project ([MaNGAObj](https://github.com/jlsteffen/MaNGAObj)), I built a catalog of all ancillary objects within the survey. I will use the galaxies identified in the catalog as a starting point for my sample. 

## Star Formation Rates
The star formation rate is the number of stars that are created in a galaxy or a region of a galaxy every year. While this sounds like a simple definition, directly measuring this for galaxies is difficult. Inidividual stars in a galaxy cannot be separately resolved from the billions of other stars so astronomers cannot simply count the number of new stars in a galaxy. We instead have to rely on tracers of new stars to infer star formation rates.

When stars are born, they are not created one at a time. They are created in groups in what are called stellar nurseries. The new population of stars are not identical and will have varying masses. The most massive stars, O and B type stars, are very lumoinous but will also have short lifespans (only a few million years). This means that if we detect tracers of OB stars, the galaxy or region must have recently created new stars. One tracer of OB stars with optical spectra is the H-alpha emission line shown in the figure below. This emission line is created by hydrogen gas which has been ionized by nearby OB stars. The luminosity of the line is directly related to the star formation rate. That is, galaxies with stronger H-alpha luminosites are creating more new stars. 

![Figure 2](https://github.com/jlsteffen/Star-Formation-in-Mergers/blob/main/images/spectra.jpg)

Unfortunetely, there is a caveat to using the H-alpha line to infer star formation rates. There are other sources that can ionize hydrogen gas. The main one is an active galactic nucleus (AGN). An AGN is where the region around a galaxy's central supermassive black hole emits high energy photons because the black hole is accreting large amounts of material. Another sources is from hot low-mass evolved stars which may also be able to ionize hydrogen gas. Fortunetely, these ionization sources can be separted from each other using the galaxy's spectra. 

We use the figure below to select galaxies how emission lines are generated by star formation. The left panel shows the color-magnitude diagram for our galaxies. The plot shows the galaxies' color on the y-axis, bluer galaxies are further down on the figure while redder galaxier are further up on the figure. Massive stars with short lifespans will be bluer while small stars with long lifespans will be redder. That means that red galaxies are dominated by older less massive stars while blue galaxies are dominated by young massive stars. Since we are looking for young stars, we require that our galaxies fall below the lower dashed line (the blue cloud) to be in our sample. 

The middle panel shows the BPT diagram (Baldwin et al. 1981) of our galaxies. The panel plots two emission line ratios on its x and y axis which creates this "seagull-like" distribution. Each wing represents a different source of ionization. The right branch represents galaxies whose dominant ionization source is an AGN while the left branch represents galaxies whose dominant ionization source are OB stars. We select galaxies which fall below the dashed line (the left branch) to be in our sample. 

The right panel shows the WHAN diagram (Cid Fernandes et al. 2011) of our galaxies. The panel plots an emission line ratio on the x axis against the equivalent width of the H-alpha line on the y axis. Cid Fernandes et al. 2011 has shown that the equivalent width of the H-alpha line is a good tracer for hot low mass evolved stars and that galaxies whose emission lines are dominated by these objects can be separated by imposing an H-alpha equivalent width cut at 3 Angstroms. Using these three methods, we create a sample of star forming galaxies. 

![Figure 3](https://github.com/jlsteffen/Star-Formation-in-Mergers/blob/main/images/spectra.jpg)

## Sample Selection
With a set of star forming galaxies identified in the survey, we need to determine which galaxies are actually gravitationally paired with the target galaxies. While two galaxies may appear close to each other on the sky, there may be a significant distance between the two objects along our line-of-sight. To verify that they are paired galaxies, the two galaxies need to be near each other on the sky, have a similar redshift (which tells us how far the galaxies are from us), and have a limited relative velocity between the two galaxies. Since we want to study how galaxy mergers influence star formation, we also create a sample of isolated control galaxies which have no nearby companions. In total, we have 169 paired galaxies and 1830 control galaxies to work with. 

## Radial Profiles
MaNGA provides a 2D distribution of galaxy spectra for each of its observations. To prep the data for comparison, we will reduce this 2D distribution to a 1D profile as a function of the galaxy's radius. Typical galaxies are circular; however, if they are at an angle with the observer they will appear to be an ellipse on the sky. The geometry of the galaxies needs to be deprojected to account for this effect. I show this process in the plot below. The left panel shows one of the survey's galaxies with elliptical profiles from various previous surveys overlaid. I calculate the inclination angle for the galaxies using the major-to-minor axis ratio of the shown ellipses. The inclination angle is then used to calculate each pixel's radius as shown in the middle panel. In the right panel we now show the star formation rate as a function of galaxy radius. The black squares represent individual pixels while the red line shows the average star formation rate within discrete radius bins. We create these star formation profiles for each galaxy in our pair and control samples. 

![Figure 4]([https://github.com/jlsteffen/Star-Formation-in-Mergers/blob/main/images/ssfr_comb.jpg](https://github.com/jlsteffen/Star-Formation-in-Mergers/blob/main/images/8332-12702.pdf)

## Analysis
With the two samples constructed and radial profiles created, we compare the star formation rate data between the pair and control samples.

![Figure 5](https://github.com/jlsteffen/Star-Formation-in-Mergers/blob/main/images/ssfr_comb.jpg)

We see that the paired galaxies have higher rates of sSFR in their centers than the control galaxies but the two samples have similar rates of sSFR in their disks.  

This work has been published in [Steffen et al. 2021](https://ui.adsabs.harvard.edu/abs/2021ApJ...909..120S/abstract).
