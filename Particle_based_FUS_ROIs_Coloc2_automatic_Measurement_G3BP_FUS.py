import os
from ij import IJ
from ij.io import DirectoryChooser
from java.io import File
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.gui import WaitForUserDialog as WFU
from ij import ImagePlus as IP

# directory chooser
choose = DirectoryChooser("Choose an image")
path = choose.getDirectory()

files = os.listdir(path) # create list of files in directory

list_green = [] # list for green images
list_red = [] # list for red images

count = 0

# making a reference for the RoiManager to work later with it
rm = RoiManager().getInstance()

for i in range(len(files)):

	if "sdc488" in files[i]:
		list_green.append(files[i]) # append green image path to img_green
		
	if "sdc640" in files[i]:
		list_red.append(files[i]) # append red image path to img_red

for i in range(len(list_green)):

	rm.reset()

	count += 1
	
	img_green = IJ.openImage(path+list_green[i]) # open green image
	img_green.show() # show image
	img_green_name = list_green[i] # variable to store image name

	img_red = IJ.openImage(path+list_red[i]) # open red image
	img_red.show() # show image
	img_red_name = list_red[i] # variable to store image name

	# set correct scale for 60x SDC
	IJ.run(img_green, "Set Scale...", "distance=1 known=0.217 pixel=1 unit=um global")
	IJ.run(img_red, "Set Scale...", "distance=1 known=0.217 pixel=1 unit=um global")

	# Filtering with gaussian blur
	IJ.run(img_green, "Gaussian Blur...", "sigma=1")
	IJ.run(img_red, "Gaussian Blur...", "sigma=1")


	# Subtract background using roling ball algorithm
	IJ.run(img_green, "Subtract Background...", "rolling=100 disable")
	IJ.run(img_red, "Subtract Background...", "rolling=100 disable")

	# clear outside of ROI
	IJ.run(img_green, "Clear Outside", "")
	IJ.run(img_red, "Clear Outside", "")

	# duplicate image for intensity analysis
	img_green_duplicate = img_green.duplicate()
	img_red_duplicate = img_red.duplicate()
	dup_name_green = IP.getTitle(img_green_duplicate)	
	dup_name_red = IP.getTitle(img_red_duplicate)	
	img_green_duplicate.show()
	img_red_duplicate.show()

	# run otsu auto threshold
	IJ.run(img_green, "Auto Threshold", "method=Otsu ignore_black")
	IJ.run(img_red, "Auto Threshold", "method=Otsu ignore_black")
	
	# Analyse particles
	IJ.run(img_green, "Analyze Particles...", "size=0.2-35 display exclude add clear")
	
	# clear results from results table
	IJ.run("Clear Results", "")

	ROI_count = rm.getCount()

	# if there are ROIs in the manager, start measuring procedure
	if ROI_count > 0:
		
		for ROI in range(ROI_count): 
			
			# iterate over ROIs and measure each one
			rm.select(img_green_duplicate, ROI)
			# Run the coloc2 plugin for colocalisation analysis
			IJ.run("Coloc 2", "channel_1=["+dup_name_green+"] channel_2=["+dup_name_red+"] roi_or_mask=[ROI(s) in channel 1] threshold_regression=Bisection display_images_in_result li_histogram_channel_1 li_histogram_channel_2 li_icq manders'_correlation 2d_intensity_histogram psf=1 costes_randomisations=5")
	
	IJ.run("Close All", "")

	# for processing of large amounts of cells
	# wait for user to close pdf result windows
	if count > 100:
		WFU("Please Close All PDF Windows").show()
		count = 0

	print(count)
	
# Save log results as a csv file in path and then close it
# This method will save ALL data from ALL images into ONE csv!!
IJ.selectWindow("Log")
IJ.saveAs("Text", path+"/Particle_based_Coloc2_FUS_ROIs.csv") 
IJ.run("Close")

	
IJ.run("Close All", "")
print("Finished task")