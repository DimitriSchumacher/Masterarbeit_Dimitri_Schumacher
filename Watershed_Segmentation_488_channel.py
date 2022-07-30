import os
from ij import IJ
from ij import WindowManager
from ij.io import DirectoryChooser
from java.io import File
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable

choose = DirectoryChooser("Choose an image")
path = choose.getDirectory()

files = os.listdir(path) # create list of files in directory

list_green = [] # list for green images

# making a reference for the RoiManager to work later with it
rm = RoiManager().getInstance()

for i in range(len(files)):

	if "sdc488" in files[i]:
		list_green.append(files[i]) # append green image path to img_green

for i in range(len(list_green)):
	img_green = IJ.openImage(path+list_green[i]) # open green image
	img_green.show() # show image
	
	# Set image scale
	IJ.run(img_green, "Set Scale...", "distance=1 known=0.217 pixel=1 unit=um global")

	# This is done to later transpose the ROIs onto the original image (duplicate)
	img_green_duplicate = img_green.duplicate()
	img_green_duplicate.show()
	
	# Function to get a list of the active image titles 
	# in this case original and duplicate
	titles = WindowManager.getImageTitles()
	
	# Filter image using Gausian-Blur
	IJ.run(img_green, "Gaussian Blur...", "sigma=2")
	
	# define histogram stats and characteristics
	stats = img_green.getStatistics()
	standard_Dev = stats.stdDev
	mean = stats.mean
	
	x = 0 # How many standard deviations away from the mean should the cutoff be

	# define cutoff point for thresholding, x std. Devs away from mean
	cutoff = mean + x * standard_Dev

	# set threshold according to the cutoff variable
	IJ.setRawThreshold(img_green, cutoff, 65535, "null")
	IJ.run(img_green, "Convert to Mask", "")

	# Fill holes of the nuclei   
	IJ.run(img_green, "Fill Holes", "")

	# Watershed
	IJ.run(img_green, "Watershed", "")
	
	# Analyse particles
	IJ.run(img_green, "Analyze Particles...", "size=128-1000 display exclude clear add")

	# select all ROIs and save them in path
	ROIset = path+"ROIs_"+list_green[i]+".zip"
	#rm.runCommand("Save", ROIset)

	ROI_count = rm.getCount()

	# iterate over all ROIs, crop and save them
	for ROI in range(ROI_count):
		ROI_number = str(ROI)
		ROI_duplicate = img_green_duplicate.duplicate()
		ROI_duplicate.show()
		rm.select(ROI_duplicate, ROI)
		image_crop = IJ.run(ROI_duplicate, "Duplicate...", " ")
		IJ.saveAs(image_crop, "Tiff", path+"Single_Cell_"+list_green[i][0: -4]+"_Cell_"+ROI_number+".tif")
		
	IJ.run("Close All", "")

print("Finished")