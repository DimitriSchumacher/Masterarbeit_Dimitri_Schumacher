import os
from ij import IJ
from ij import WindowManager
from ij.io import DirectoryChooser
from java.io import File
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
import csv

choose = DirectoryChooser("Choose an image")
path = choose.getDirectory()

files = os.listdir(path) # create list of files in directory

list_green = [] # list for green images
list_red = []  # list for red images

green_results = []

# making a reference for the RoiManager to work later with it
rm = RoiManager().getInstance()

# reference of results table to work later with it
rt = ResultsTable.getResultsTable()

# variable for heading of results table
headings = rt.getHeadings()

# setup of measure tool
IJ.run("Set Measurements...", "area mean min shape median integrated display redirect=None decimal=3")

for i in range(len(files)):

	if "sdc488" in files[i]:
		list_green.append(files[i]) # append green image (GFP, 488 nm) path to img_green
	elif "sdc640" in files[i]:
		list_red.append(files[i]) # append red image (G3BP, 640 nm) path to img_red

for i in range(len(list_green)):

	rm.reset()
	
	
	img_green = IJ.openImage(path+list_green[i]) # open green image
	img_green.show() # show image
	
	# Set image scale
	IJ.run(img_green, "Set Scale...", "distance=1 known=0.217 pixel=1 unit=um global")

	# Filter image using Gausian-Blur 
	IJ.run(img_green, "Gaussian Blur...", "sigma=1")

	# Subtract backround
	IJ.run(img_green, "Subtract Background...", "rolling=100 disable")

	# Clear image outside of ROI
	IJ.run(img_green, "Clear Outside", "")

	# duplicate image for intensity analysis
	img_green_duplicate = img_green.duplicate()
	
	# run otsu auto threshold
	IJ.run(img_green, "Auto Threshold", "method=Otsu ignore_black")
	
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
			IJ.run(img_green_duplicate, "Measure", "")

		# append results into green_results list
		for row in range(rt.size()): 
			green_results.append(rt.getRowAsString(row))

		headings = rt.getHeadings()
		
	IJ.run("Close All", "")

heading = "Number" # name of first column 

# create complete heading for the data
for i in headings:
	heading = heading + "\t" + i

# create and write csv file
file_green = open(path + "Particle_Measurements_488.csv", "w")
writer_green = csv.writer(file_green, delimiter="\t", lineterminator="\n")

file_green.write(heading + "\n")

for row in green_results:
	file_green.write(row.encode("latin1") + "\n")
	
file_green.close()

print("Finished")