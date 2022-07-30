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

# Results for percentage of cells with aggregate
# Binary count: if aggr. in cell = 1. If not = 0

results = []

# making a reference for the RoiManager to work later with it
rm = RoiManager().getInstance()

for i in range(len(files)):

	if "sdc488" in files[i]:
		list_green.append(files[i]) # append green image (GFP, 488 nm) path to img_green
	elif "sdc640" in files[i]:
		list_red.append(files[i]) # append red image (G3BP, 640 nm) path to img_red

for i in range(len(list_green)):

	FUS_aggr = 0
	G3BP_aggr = 0
	both_aggr = 0
	name = list_green[i]
		
	img_green = IJ.openImage(path+list_green[i]) # open green image
	img_green.show() # show image
	
	img_red = IJ.openImage(path+list_red[i]) # open green image
	img_red.show() # show image
	
	# Set image scale
	IJ.run(img_green, "Set Scale...", "distance=1 known=0.217 pixel=1 unit=um global")
	IJ.run(img_red, "Set Scale...", "distance=1 known=0.217 pixel=1 unit=um global")

	# Filter image using Gausian-Blur 
	IJ.run(img_green, "Gaussian Blur...", "sigma=1")
	IJ.run(img_red, "Gaussian Blur...", "sigma=1")

	# Subtract backround
	IJ.run(img_green, "Subtract Background...", "rolling=100 disable")
	IJ.run(img_red, "Subtract Background...", "rolling=100 disable")

	# Clear image outside of ROI
	IJ.run(img_green, "Clear Outside", "")
	IJ.run(img_red, "Clear Outside", "")

	# set otsu auto threshold
	IJ.run(img_green, "Auto Threshold", "method=Otsu")# ignore_black")
	IJ.run(img_red, "Auto Threshold", "method=Otsu")# ignore_black")
	
	# Analyse FUS particles
	IJ.run(img_green, "Analyze Particles...", "size=0.2-35 display exclude clear add")

	ROI_count = rm.getCount()

	# checks if there is FUS aggr. in the cell
	if ROI_count > 0:
		
		FUS_aggr = 1
		rm.reset()
		
	# Analyse G3BP particles
	IJ.run(img_red, "Analyze Particles...", "size=0.2-35 display exclude clear add")
	
	ROI_count = rm.getCount()

	# checks if there is G3BP aggr. in the cell
	if ROI_count > 0:
		
		G3BP_aggr = 1
		rm.reset()

	if G3BP_aggr == 1 and FUS_aggr == 1:

		both_aggr = 1


	results.append([name, FUS_aggr, G3BP_aggr, both_aggr])

	IJ.run("Close All", "")
	
# print(results) # for debugging

heading = ["Name","FUS aggr.","G3BP aggr.","Both aggr."]

# create and write csv file
results_file = open(path + "Results_Binary_Cell_Aggr_Counter.csv", "w")
writer_results = csv.writer(results_file, delimiter=",", lineterminator="\n")

writer_results.writerow(heading)

for row in results:
	writer_results.writerow(row)
	
results_file.close()

print("Finished")