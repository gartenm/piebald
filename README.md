# piebald

Hi,

this is code accompanying a publication by Garten et al. “Contacting domains segregate a lipid transporter from a solute transporter in the malarial host-parasite interface”. Please feel encouraged to contact the authors with any question about the scripts.

The scripts in folder "projector" produce Mollweide projections (see https://en.wikipedia.org/wiki/Mollweide_projection) from 3d voxel data. Those projections can be analyzed with scripts in the folder analyze. Scripts found in the folder distance are made to perform distance measurements on segmented images. Scripts in the bootstrap folder perform a bootstrap analysis of membrane distance data. One image example for the projection is provided with the scripts. Data for the bootstrap can be found with the publication.

The projection scripts are written for matlab 2018b with gpu acceleration. If the gpu package is not available to you, change the gpuarrays to regular arrays.

The scripts require functions that are freely available elsewhere:
The bioformats plugin can be found on https://www.openmicroscopy.org/bio-formats/

Finally, a little more detailed description of the individual files. 

Projector:
batchMap
	use this to go directly from the stack to the projection
	it needs the bioformats plug in to load the images - https://docs.openmicroscopy.org/bio-formats/6.1.0/users/matlab/index.html
	makes projections from czi or tiff stacks
	chose the right version depending on how many channels are in the image and which channel should be used to find the center

findCenter
	finds the center of the object in the image stack	

piebaldmap
	calculates the angular projection (and a cylidrical projection)
	the lower the angular resulotion the faster
	uses gpu features (can be edited out)

plotPol3d
	draws the mollweider projection (and pol projections)
	uses gpu features (can be edited out)

Analyze:
mapcolocCostes
	calculates the colocalization coeff from the projections
Distance:
measureLineDistance
	looks for pixels in 2nd layer of image and looks for closest pixel on 3rd layer
	uses scaling information of the image to get pyhsical distances

Bootstrap:
AnalyzeDistro.m
	used to fit histrogram frequency data, log transforms cdf
	full 2 component mixture model
	called in bootdtrapping functions

AnalyzeDistroGold.m
	as AnalyzeDistro.m but mu2 is fixed

plotMeasureLineDistanceV2.m
	makes histrograms from data
	saves summary and histrogram data, makes frequency data from histrograms (Avector) to make CDFs from scaled data, B are unscaled datapoints

bootstrapPVMPPMdist.m
	for PVMPPM membrane distributions
	bootstrapps image wise
	scales for pixel sizes

bootstrapallPVMPPMdistFlexiblePVlength.m
	for PVMPPM membrane distributions to compare to the gold distribution
	bootstrapps segment wise (from the scale), segment size can be set

bootstrapgold.m
	bootstrapps gold distribution, fits all parameters of mixture model

bootstrapgold2.m
	bootstrapps gold distribution, 2nd mean is fixed in AnalyzeDistroGold.m

bootstrapallPVMPPMdist.m
	bootstrapps data point wise when no scaling is needed

bootstrapPVMPPMdistNoResampleCompareDist.m
	uses a mix of distributions to define significance of mixture parameter a

CI95.m
	gives the central 95% of a histogram

plotA.m
	makes plots for the gold label figure

plotAA.m
	makes the histrogram graphs for the PfNCR1 KD figure
