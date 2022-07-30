import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import filedialog
from scipy.stats import ttest_ind
from scipy.stats import kstest
from scipy.stats import median_test
import csv

# File chooser dialog
# Please select the csv file containing the cell measurements
root = tk.Tk()
root.withdraw()
path = filedialog.askopenfilename()

df = pd.read_csv(path, sep="\t", encoding='latin1') # create dataframe from csv file in path

stats = np.array(df.values) # creates numpy array out of df

# Empty list to append the needed data for analysis
# data = [[name],[cell_type],[condition], [Area], [Mean Intensity],[Dataset],[Time],[intDen]]
data = [[],[],[],[],[],[],[],[]]

mutants = ["WT", "R514S", "R521C"]
conditions = ["Control", "10 $\mu$M", "1 mM"]
time = ["24 h", "48 h"]
dataset = ["0","1","2"]

cutoff_area = 35
cutoff_intensity = 500
cutoff_intDen = 0

def create_data():

    # function to create data for analysis from df dataframe
    
    # dataframe slices 
    name = df["Label"]
    mean = df["Mean"]
    area = df["Area"]
    intDen = df["IntDen"]
    
    for i in range(len(name)):
        
        if area[i] < cutoff_area and mean[i] > cutoff_intensity and intDen[i] > cutoff_intDen:
            
            data[0].append(name[i])
            
            data[3].append(area[i])
            data[4].append(mean[i])
            data[7].append(intDen[i])
                 
            if "wt" in name[i] or "WT" in name[i]:
                data[1].append("WT")
            elif "514" in name[i]:
                data[1].append("R514S")
            elif "521" in name[i]:
                data[1].append("R521C")
            
            # append condition
            if "ohne" in name[i] or "no_peptide" in name[i]:
                data[2].append("Control")
            elif "_10" in name[i]:
                data[2].append("10 $\mu$M")
            elif "1mM" in name[i]:
                data[2].append("1 mM")
                
            # append dataset
            if "210818" in name[i]: 
                data[5].append("0")
            elif "211015" in name[i] or "211016" in name[i]: 
                data[5].append("1")
            elif "211021" in name[i] or "211022" in name[i]:    
                data[5].append("2")   
            
            # append timepoint
            if "24h" in name[i]:
                data[6].append("24 h")
            elif "48h" in name[i]:
                data[6].append("48 h")
        
    
create_data()

# create dataframe from data list, transposing it beforehand to fit formating
df_data = pd.DataFrame(np.array(data).T.tolist())

labels0 = []

def create_labels():

    for mutant in mutants: 
        
        for t in time: 
            
            for condition in conditions:
            
                labels0.append(mutant+" "+condition + " " + t)

create_labels()

def stat_test(test, alpha, parameter):
    
    # function to do statistical test (t-test, mood's-test or KS-test)
    # input: alpha (e.g. 0.05) and test="t" for t-test, test="mood" for mood's test or test="KS" for KS-test
    # parameter: "intDen", "intensity" or "area", to test either one

    print("--"*20)
    print("Statistical Test: " + test + "-Test of the " + parameter)
    print("--"*20)
    
    if parameter == "intensity":
        
        value = 4
        
    elif parameter == "area":
        
        value = 3
        
    elif parameter == "intDen":
        
        value = 7
        
    for mutant in mutants: 

        a = df_data[(df_data[1] == mutant)]
        
        for condition in conditions: 

            b = a[(a[2] == condition)]
            
            for t in time: 
                
                b2 = b[(b[6] == t)]
                
                for mutant2 in mutants:
                        
                    c = df_data[(df_data[1] == mutant2)]
                    
                    if mutant == mutant2:
                        
                        for condition2 in conditions: 
                            
                            d = c[(c[2] == condition2)]
                            
                            for t2 in time: 
                                
                                d2 = d[(d[6] == t2)]
                                
                                #print(d2)
                                
                                if test == "t": 
                                    
                                    res = ttest_ind(list(map(float, b2[value].values)), list(map(float, d2[value].values)), equal_var=False).pvalue
                                    
                                elif test == "KS":
                                    
                                    res = kstest(list(map(float, b2[value].values)), list(map(float, d2[value].values))).pvalue
                                    
                                elif test == "mood": 
                                                    
                                    res = median_test(list(map(float, b2[value].values)), list(map(float, d2[value].values)))[1]
                                    
                                if res < alpha: 
                                    
                                    print(mutant + " " + condition + " " + t + " vs " + mutant2 + " " + condition2 + " " + t2 + " " + test + "-Test")
                                    print(res)
                                        
                    elif mutant != mutant2:
                        continue

    print("--"*20)
    
#stat_test(test="t", alpha=0.05, parameter="area")

def print_stats():
    
    # function to print basic stats    
    
    print("--"*20)    
    print("Stats Printer")
    print("--"*20)    
    
    number_of_particles = []
    
    for mutant in mutants: 

        a = df_data[(df_data[1] == mutant)]
        
        for t in time: 
            
            b = a[(a[6] == t)]
        
            for condition in conditions: 
    
                c = b[(b[2] == condition)]
                
                particle_amount = len(c)
                number_of_particles.append(particle_amount)
    
    print("The stats are printed in the following order: ")
    print(labels0)
    print("--"*20)  
    print("Number of particles for each condition: ")
    print(number_of_particles)
    print("--"*20)    

print_stats()


norm_data_intensity = []
norm_data_area = []
norm_data_labels = []
norm_data_intDen = []

def normalize_data():    
    
    # Function to normalize values of intensity, intDen and area
    # Values are normalized to the mean of the 48 h control
    
    # means to be used as normalization values (48 h control mean)
    norm_values_intensity = []
    norm_values_area = []
    norm_values_intDen = []
    
    for mutant in mutants: 
        
        a0 = df_data[(df_data[1] == mutant)]
        
        b0 = a0[(a0[2] == "Control")]
        
        c0 = b0[(b0[6] == "48 h")]
        
        mean_intensity_ctrl_48h_mean = np.mean(list(map(float, c0[4].values)))
        norm_values_intensity.append(mean_intensity_ctrl_48h_mean)
        
        area_ctrl_48h_mean = np.mean(list(map(float, c0[3].values)))
        norm_values_area.append(area_ctrl_48h_mean)
        
        intDen_ctrl_48h_mean = np.mean(list(map(float, c0[7].values)))
        norm_values_intDen.append(intDen_ctrl_48h_mean)

    count = -1

    for mutant in mutants: 
        
        count += 1

        a = df_data[(df_data[1] == mutant)]
        
        for t in time: 
            
            b = a[(a[6] == t)]
        
            for condition in conditions: 
    
                c = b[(b[2] == condition)]
            
                c2 = list(map(float, c[4].values))
                c3 = list(map(float, c[3].values))
                c4 = list(map(float, c[7].values))
                
                division_intensities = np.divide(c2, norm_values_intensity[count])
                division_area = np.divide(c3, norm_values_area[count])
                division_intDen = np.divide(c4, norm_values_intDen[count])
                
                norm_data_intensity.append(division_intensities)
                norm_data_area.append(division_area)
                norm_data_intDen.append(division_intDen)
                
                label_data = str(mutant + " " + condition + " " + t)            
                norm_data_labels.append(label_data)


normalize_data()

norm_data_intensity2 = []
norm_data_area2 = []
norm_data_intDen2 = []

norm_data_intensity2_summarized = []
norm_data_area2_summarized = []
norm_data_intDen2_summarized = []

def normalize_data2():    
    
    # Function to normalize values of intensity, intDen and area
    # The data is normalized for each dataset sepparately
    # Values are normalized to the mean of the 48 h control
    
    # means to be used as normalization values (48 h control mean)
    norm_values_intensity = []
    norm_values_area = []
    norm_values_intDen = []
    
    for mutant in mutants: 
        
        for datas in dataset:
            
            a0 = df_data[(df_data[1] == mutant)]
            
            b0 = a0[(a0[2] == "Control")]
            
            c0 = b0[(b0[6] == "48 h")]
            
            d0 = c0[(c0[5] == datas)]
            
            mean_intensity_ctrl_48h_mean = np.mean(list(map(float, d0[4].values)))
            norm_values_intensity.append(mean_intensity_ctrl_48h_mean)
            
            area_ctrl_48h_mean = np.mean(list(map(float, d0[3].values)))
            norm_values_area.append(area_ctrl_48h_mean)
            
            intDen_ctrl_48h_mean = np.mean(list(map(float, d0[7].values)))
            norm_values_intDen.append(intDen_ctrl_48h_mean)

    count = -1

    for mutant in mutants: 

        a = df_data[(df_data[1] == mutant)]
        
        for datas in dataset: 
            
            count += 1
            
            b0 = a[(a[5] == datas)]
        
            for t in time: 
                
                b = b0[(b0[6] == t)]
            
                for condition in conditions: 
                    
                    c = b[(b[2] == condition)]
                
                    c2 = list(map(float, c[4].values))
                    c3 = list(map(float, c[3].values))
                    c4 = list(map(float, c[7].values))
                    
                    division_intensities = np.divide(c2, norm_values_intensity[count])
                    division_area = np.divide(c3, norm_values_area[count])
                    division_intDen = np.divide(c4, norm_values_intDen[count])
                    
                    norm_data_intensity2.append(division_intensities)
                    norm_data_area2.append(division_area)
                    norm_data_intDen2.append(division_intDen)
    
    for i in range(len(mutants)):
        
        for i2 in range(len(dataset)*len(time)):
            
            concatenated_intensity = norm_data_intensity2[i*18+i2].tolist() + norm_data_intensity2[i*18+i2+6].tolist() + norm_data_intensity2[i*18+i2+12].tolist()
            norm_data_intensity2_summarized.append(concatenated_intensity)
            
            concatenated_area = norm_data_area2[i*18+i2].tolist() + norm_data_area2[i*18+i2+6].tolist() + norm_data_area2[i*18+i2+12].tolist()
            norm_data_area2_summarized.append(concatenated_area)
            
            concatenated_intDen = norm_data_intDen2[i*18+i2].tolist() + norm_data_intDen2[i*18+i2+6].tolist() + norm_data_intDen2[i*18+i2+12].tolist()
            norm_data_intDen2_summarized.append(concatenated_intDen)
            
    
normalize_data2()

def norm_stat_test(test, alpha, parameter):
    
    # function to do statistical test (t-test, mood's-test or KS-test) on normalized values to each dataset (normalize_data2 function)
    # input: alpha (e.g. 0.05) and test="t" for t-test, test="mood" for mood's test or test="KS" for KS-test
    # parameter: "intensity", "intDen" or "area", to test either one or the other
    
    print("--"*20)
    print("Statistical test of normalized data: " + test + "-test, " + parameter)
    print("--"*20)
    
    if parameter == "intensity":
        
        for i in range(len(conditions)*len(time)):
            
            for i2 in range(len(conditions)):
                
                x1 = 0
                x2 = 0
                
                if i2 == 0:
                    
                    x1 = 1
                    x2 = 2
                
                if test == "t":
                    
                    res =  ttest_ind(norm_data_intensity2_summarized[i*3+x1], norm_data_intensity2_summarized[i*3+i2+x2], equal_var=False).pvalue
                
                elif test == "mood":
                    
                    res = median_test(norm_data_intensity2_summarized[i*3+x1], norm_data_intensity2_summarized[i*3+i2+x2])[1]
                    
                elif test == "KS":
                    
                    res = kstest(norm_data_intensity2_summarized[i*3+x1], norm_data_intensity2_summarized[i*3+i2+x2]).pvalue
                    
                if res < alpha:
        
                    print(norm_data_labels[i*3+x1] + " vs. " + norm_data_labels[i*3+i2+x2])
                    print(res)
                    
    elif parameter == "area":
        
        for i in range(len(conditions)*len(time)):
            
            for i2 in range(len(conditions)):
                
                x1 = 0
                x2 = 0
                
                if i2 == 0:
                    
                    x1 = 1
                    x2 = 2
                
                if test == "t":
                    
                    res =  ttest_ind(norm_data_area2_summarized[i*3+x1], norm_data_area2_summarized[i*3+i2+x2], equal_var=False).pvalue
                
                elif test == "mood":
                    
                    res = median_test(norm_data_area2_summarized[i*3+x1], norm_data_area2_summarized[i*3+i2+x2])[1]
                    
                elif test == "KS":
                    
                    res = kstest(norm_data_area2_summarized[i*3+x1], norm_data_area2_summarized[i*3+i2+x2]).pvalue
                    
                if res < alpha:
        
                    print(norm_data_labels[i*3+x1] + " vs. " + norm_data_labels[i*3+i2+x2])
                    print(res)
    
    elif parameter == "intDen":
        
        for i in range(len(conditions)*len(time)):
            
            for i2 in range(len(conditions)):
                
                x1 = 0
                x2 = 0
                
                if i2 == 0:
                    
                    x1 = 1
                    x2 = 2
                
                if test == "t":
                    
                    res =  ttest_ind(norm_data_intDen2_summarized[i*3+x1], norm_data_intDen2_summarized[i*3+i2+x2], equal_var=False).pvalue
                
                elif test == "mood":
                    
                    res = median_test(norm_data_intDen2_summarized[i*3+x1], norm_data_intDen2_summarized[i*3+i2+x2])[1]
                    
                elif test == "KS":
                    
                    res = kstest(norm_data_intDen2_summarized[i*3+x1], norm_data_intDen2_summarized[i*3+i2+x2]).pvalue
                    
                if res < alpha:
        
                    print(norm_data_labels[i*3+x1] + " vs. " + norm_data_labels[i*3+i2+x2])
                    print(res)
    
    print("--"*20)
    
norm_stat_test(test="t", alpha=0.05, parameter="area")
                        
def create_stats_csv():
    
    # creates csv files with the mean and median values for: mean intensity, area, intDen
    # function for normalize_data function
    
    header = ["Labels", "IntDen Mean", "IntDen Median", "Area Mean", "Area Median", "Mean Intensity Mean", "Mean Intensity Median"]
    csv_data = []
    
    norm_data_intDen_mean = []
    norm_data_area_mean = []
    norm_data_intensity_mean = []
    
    norm_data_intDen_median = []
    norm_data_area_median = []
    norm_data_intensity_median = []
    
    for x in range(len(labels0)):
        
        norm_data_intDen_mean.append(np.mean(norm_data_intDen[x]))
        norm_data_area_mean.append(np.mean(norm_data_area[x]))
        norm_data_intensity_mean.append(np.mean(norm_data_intensity[x]))
        
        norm_data_intDen_median.append(np.median(norm_data_intDen[x]))
        norm_data_area_median.append(np.median(norm_data_area[x]))
        norm_data_intensity_median.append(np.median(norm_data_intensity[x]))
    
    for x in range(len(labels0)):
        
        row_csv = [labels0[x], norm_data_intDen_mean[x], norm_data_intDen_median[x], norm_data_area_mean[x], norm_data_area_median[x], norm_data_intensity_mean[x], norm_data_intensity_median[x]]
        
        csv_data.append(row_csv)
        
    f = open(path[0:-4] + "_Stats_Summary.csv", "w", newline="")
    
    writer = csv.writer(f)
    
    writer.writerow(header)
    
    for row in csv_data:
        writer.writerow(row)
    
    f.close()
                  
#create_stats_csv()   
    
def create_stats_csv2():
    
    # creates csv files with the mean and median values for: mean intensity, area, intDen
    # function for normalize_data2 function
    
    header = ["Labels", "IntDen Mean", "IntDen Median", "Area Mean", "Area Median", "Mean Intensity Mean", "Mean Intensity Median"]
    csv_data = []
    
    norm_data_intDen_mean = []
    norm_data_area_mean = []
    norm_data_intensity_mean = []
    
    norm_data_intDen_median = []
    norm_data_area_median = []
    norm_data_intensity_median = []
    
    for x in range(len(labels0)):
        
        norm_data_intDen_mean.append(np.mean(norm_data_intDen2_summarized[x]))
        norm_data_area_mean.append(np.mean(norm_data_area2_summarized[x]))
        norm_data_intensity_mean.append(np.mean(norm_data_intensity2_summarized[x]))
        
        norm_data_intDen_median.append(np.median(norm_data_intDen2_summarized[x]))
        norm_data_area_median.append(np.median(norm_data_area2_summarized[x]))
        norm_data_intensity_median.append(np.median(norm_data_intensity2_summarized[x]))
    
    for x in range(len(labels0)):
        
        row_csv = [labels0[x], norm_data_intDen_mean[x], norm_data_intDen_median[x], norm_data_area_mean[x], norm_data_area_median[x], norm_data_intensity_mean[x], norm_data_intensity_median[x]]
        
        csv_data.append(row_csv)
        
    f = open(path[0:-4] + "_Stats_Summary_ea_dataset.csv", "w", newline="")
    
    writer = csv.writer(f)
    
    writer.writerow(header)
    
    for row in csv_data:
        writer.writerow(row)
    
    f.close()
                  
#create_stats_csv2()   

def plot_violinplots_intensity():
    
    pos = [1,2,3,4,5,6]
    count = -1
    
    fig, axs = plt.subplots(1,3, sharey=True)
    fig.suptitle("QSYSQ: Mean Particle Intensity")
    axs[0].set_ylabel("Mean Intensity [a.u.]")
    #axs[0].set_ylim(-1.2,1.5)
        
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=80)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels0[0+x*6:6+x*6])
    
    
    for mutant in mutants: 
        
        count += 1
        count2 = -1
        
        a = df_data[(df_data[1] == mutant)]

        for condition in conditions: 
            
            count2 += 1
            
            b = a[(a[2] == condition)]
   
            c = b[(b[6] == "24 h")]
            c2 = b[(b[6] == "48 h")]
            
            axs[count].violinplot(list(map(float, c[4].values)), positions=[pos[count2]])#, showmedians=True)
            axs[count].violinplot(list(map(float, c2[4].values)), positions=[pos[count2+3]])#, showmedians=True)


def plot_boxplots_intensity():
    
    pos = [1,2,3,4,5,6]
    count = -1
    
    fig, axs = plt.subplots (1,3, sharey=True)    
    fig.suptitle("QSYSQ: Mean Particle Intensity")
    axs[0].set_ylabel("Mean Intensity [a.u.]")
    axs[0].set_ylim(-1000,30000)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=80)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels0[0+x*6:6+x*6])
    

    for mutant in mutants: 
        
        count += 1
        count2 = -1
        
        a = df_data[(df_data[1] == mutant)]

        for condition in conditions: 
            
            count2 += 1
            
            b = a[(a[2] == condition)]
   
            c = b[(b[6] == "24 h")]
            c2 = b[(b[6] == "48 h")]
            
            axs[count].boxplot(list(map(float, c[4].values)), positions=[pos[count2]], showfliers=False, showmeans=True, manage_ticks=False)
            axs[count].boxplot(list(map(float, c2[4].values)), positions=[pos[count2+3]], showfliers=False, showmeans=True, manage_ticks=False)

def plot_violinplots_area():
    
    pos = [1,2,3,4,5,6]
    count = -1
    
    fig, axs = plt.subplots(1,3, sharey=True)
    fig.suptitle("QSYSQ: Particle Area")
    axs[0].set_ylabel("Area [$\mu m^2$]")
    #axs[0].set_ylim(-1.2,1.5)
        
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=80)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels0[0+x*6:6+x*6])
    
    for mutant in mutants: 
        
        count += 1
        count2 = -1
        
        a = df_data[(df_data[1] == mutant)]

        for condition in conditions: 
            
            count2 += 1
            
            b = a[(a[2] == condition)]
   
            c = b[(b[6] == "24 h")]
            c2 = b[(b[6] == "48 h")]
            
            axs[count].violinplot(list(map(float, c[3].values)), positions=[pos[count2]])#, showmedians=True)
            axs[count].violinplot(list(map(float, c2[3].values)), positions=[pos[count2+3]])#, showmedians=True)


def plot_boxplots_area():
    
    pos = [1,2,3,4,5,6]
    count = -1
    
    fig, axs = plt.subplots (1,3, sharey=True)    
    fig.suptitle("QSYSQ: Particle Area")
    axs[0].set_ylabel("Area [$\mu m^2$]")
    #axs[0].set_ylim(-1.2,1.5)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=80)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels0[0+x*6:6+x*6])
    
    for mutant in mutants: 
        
        count += 1
        count2 = -1
        
        a = df_data[(df_data[1] == mutant)]

        for condition in conditions: 
            
            count2 += 1
            
            b = a[(a[2] == condition)]
   
            c = b[(b[6] == "24 h")]
            c2 = b[(b[6] == "48 h")]
            
            axs[count].boxplot(list(map(float, c[3].values)), positions=[pos[count2]], showfliers=False, showmeans=True, manage_ticks=False)
            axs[count].boxplot(list(map(float, c2[3].values)), positions=[pos[count2+3]], showfliers=False, showmeans=True, manage_ticks=False)

def plot_dataset_boxplots_intensity(to_plot):
    # select dataset to be plotted through to_plot variable ("0","1","2")
    
    pos = [1,2,3,4,5,6]
    count = -1
    
    fig, axs = plt.subplots (1,3, sharey=True)    
    fig.suptitle("QSYSQ: Mean Particle Intensity, Dataset " + to_plot)
    axs[0].set_ylabel("Mean Intensity [a.u.]")
    #axs[0].set_ylim(-1.2,1.5)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=80)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels0[0+x*6:6+x*6])
    
    a0 = df_data[(df_data[5] == to_plot)]
    
    for mutant in mutants: 
        
        count += 1
        count2 = -1
        
        a = a0[(a0[1] == mutant)]

        for condition in conditions: 
            
            count2 += 1
            
            b = a[(a[2] == condition)]
   
            c = b[(b[6] == "24 h")]
            c2 = b[(b[6] == "48 h")]
            
            axs[count].boxplot(list(map(float, c[4].values)), positions=[pos[count2]], showfliers=False, showmeans=True, manage_ticks=False)
            axs[count].boxplot(list(map(float, c2[4].values)), positions=[pos[count2+3]], showfliers=False, showmeans=True, manage_ticks=False)


def plot_normalized_intensities():
    
    pos = [1,2,3,4,5,6]
    count = -1
    
    fig, axs = plt.subplots (1,3, sharey=True)    
    fig.suptitle("QSYSQ: Normalized Mean Particle Intensity")
    axs[0].set_ylabel("Normalized Mean Intensity")
    axs[0].set_ylim(-0.1,4)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=80)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels0[0+x*6:6+x*6])
    
    for i in range(len(mutants)):
        
        for i2 in range(len(time)*len(conditions)):
        
            count += 1
        
            axs[i].boxplot(norm_data_intensity[count], positions=[pos[i2]], showfliers=False, showmeans=True, manage_ticks=False)
        
def plot_normalized_area():
    
    pos = [1,2,3,4,5,6]
    count = -1
    
    fig, axs = plt.subplots (1,3, sharey=True)    
    fig.suptitle("QSYSQ: Normalized Area")
    axs[0].set_ylabel("Normalized Area")
    axs[0].set_ylim(-0.1,4)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=80)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels0[0+x*6:6+x*6])
    
    for i in range(len(mutants)):
        
        for i2 in range(len(time)*len(conditions)):
        
            count += 1
        
            axs[i].boxplot(norm_data_area[count], positions=[pos[i2]], showfliers=False, showmeans=True, manage_ticks=False)
 
def plot_normalized_intensities_no_wt():
    
    mutants_no_wt = ["R514S", "R521C"]
    pos = [1,2,3,4,5,6]
    count = -1
    
    fig, axs = plt.subplots (1,2, sharey=True)    
    fig.suptitle("QSYSQ: Normalized Mean Particle Intensity")
    axs[0].set_ylabel("Normalized Mean Intensity")
    axs[0].set_ylim(-0.1,4)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=80)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels0[6+x*6:12+x*6])
    
    for i in range(len(mutants_no_wt)):
        
        for i2 in range(len(time)*len(conditions)):
        
            count += 1
        
            axs[i].boxplot(norm_data_intensity[count+6], positions=[pos[i2]], showfliers=False, showmeans=True, manage_ticks=False)
        
def plot_normalized_area_no_wt():
    
    mutants_no_wt = ["R514S", "R521C"]
    pos = [1,2,3,4,5,6]
    count = -1
    
    fig, axs = plt.subplots (1,2, sharey=True)    
    fig.suptitle("QSYSQ: Normalized Area")
    axs[0].set_ylabel("Normalized Area")
    axs[0].set_ylim(-0.1,4)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=80)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels0[6+x*6:12+x*6])
    
    for i in range(len(mutants_no_wt)):
        
        for i2 in range(len(time)*len(conditions)):
        
            count += 1
        
            axs[i].boxplot(norm_data_area[count+6], positions=[pos[i2]], showfliers=False, showmeans=True, manage_ticks=False)
            
            
def plot_boxplots_intDen():
    
    pos = [1,2,3,4,5,6]
    count = -1
    
    fig, axs = plt.subplots (1,3, sharey=True)    
    fig.suptitle("QSYSQ: Particle Integrated Density")
    axs[0].set_ylabel("IntDen [a.u.]")
    #axs[0].set_ylim(-1.2,1.5)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=80)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels0[0+x*6:6+x*6])
    
    for mutant in mutants: 
        
        count += 1
        count2 = -1
        
        a = df_data[(df_data[1] == mutant)]

        for condition in conditions: 
            
            count2 += 1
            
            b = a[(a[2] == condition)]
   
            c = b[(b[6] == "24 h")]
            c2 = b[(b[6] == "48 h")]
            
            axs[count].boxplot(list(map(float, c[7].values)), positions=[pos[count2]], showfliers=False, showmeans=True, manage_ticks=False)
            axs[count].boxplot(list(map(float, c2[7].values)), positions=[pos[count2+3]], showfliers=False, showmeans=True, manage_ticks=False)
 
def plot_normalized_intDen_no_wt():
    
    mutants_no_wt = ["R514S", "R521C"]
    pos = [1,2,3,4,5,6]
    count = -1
    
    fig, axs = plt.subplots (1,2, sharey=True)    
    fig.suptitle("QSYSQ: Normalized intDen")
    axs[0].set_ylabel("Normalized intDen")
    axs[0].set_ylim(-0.1,4)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=80)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels0[6+x*6:12+x*6])
    
    for i in range(len(mutants_no_wt)):
        
        for i2 in range(len(time)*len(conditions)):
        
            count += 1
        
            axs[i].boxplot(norm_data_intDen[count+6], positions=[pos[i2]], showfliers=False, showmeans=True, manage_ticks=False)

def plot_normalized_intDen():
    
    pos = [1,2,3,4,5,6]
    count = -1
    
    fig, axs = plt.subplots (1,3, sharey=True)    
    fig.suptitle("QSYSQ: Normalized intDen")
    axs[0].set_ylabel("Normalized intDen")
    axs[0].set_ylim(-0.1,4)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=80)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels0[0+x*6:6+x*6])
    
    for i in range(len(mutants)):
        
        for i2 in range(len(time)*len(conditions)):
        
            count += 1
        
            axs[i].boxplot(norm_data_intDen[count], positions=[pos[i2]], showfliers=False, showmeans=True, manage_ticks=False)


def scatter_plots():
    
    mutants_no_wt = ["R514S", "R521C"]
    
    count = -1
    count_2 = -1
    
    fig, axs = plt.subplots (4,3, sharex=True, sharey=True)    
    fig.suptitle("QSYSQ: Mean Intensity vs. Area")
    axs[1,0].set_ylabel("Area")
    axs[-1,1].set_xlabel("Mean Intensity")
    axs[-1,1].set_xlim(0,800)
    axs[1,0].set_ylim(0,5)
    
    for mutant in mutants_no_wt: 
        
        a = df_data[(df_data[1] == mutant)]
        
        for t in time: 
            
            count += 1
            count_2 = -1
            
            b = a[(a[6] == t)]
        
            for condition in conditions:
                
                count_2 += 1
    
                c = b[(b[2] == condition)]
                
                area = list(map(float, c[3].values))
                intensity = list(map(float, c[4].values))
                
                axs[count,count_2].scatter(intensity, area, s=0.25)

def plot_normalized_intensities2():
    
    # plots normalized data to mean of the 48h control
    # normalized for each dataset sepparately
    
    pos = [1,2,3,4,5,6]
    count = -1
    
    fig, axs = plt.subplots (1,3, sharey=True)    
    fig.suptitle("QSYSQ: Normalized Mean Particle Intensity (ea Dataset)")
    axs[0].set_ylabel("Normalized Mean Intensity")
    axs[0].set_ylim(-0.1,4)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=80)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels0[0+x*6:6+x*6])
    
    for i in range(len(mutants)):
        
        for i2 in range(len(time)*len(conditions)):
        
            count += 1
        
            axs[i].boxplot(norm_data_intensity2_summarized[count], positions=[pos[i2]], showfliers=False, showmeans=True, manage_ticks=False)

def plot_normalized_intDen2():
    
    pos = [1,2,3,4,5,6]
    count = -1
    
    fig, axs = plt.subplots (1,3, sharey=True)    
    fig.suptitle("QSYSQ: Normalized intDen (ea Dataset)")
    axs[0].set_ylabel("Normalized intDen")
    axs[0].set_ylim(-0.1,4)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=80)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels0[0+x*6:6+x*6])
    
    for i in range(len(mutants)):
        
        for i2 in range(len(time)*len(conditions)):
        
            count += 1
        
            axs[i].boxplot(norm_data_intDen2_summarized[count], positions=[pos[i2]], showfliers=False, showmeans=True, manage_ticks=False)

def plot_normalized_area2():
    
    pos = [1,2,3,4,5,6]
    count = -1
    
    fig, axs = plt.subplots (1,3, sharey=True)    
    fig.suptitle("QSYSQ: Normalized Area (ea Dataset)")
    axs[0].set_ylabel("Normalized Area")
    axs[0].set_ylim(-0.1,4)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=80)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels0[0+x*6:6+x*6])
    
    for i in range(len(mutants)):
        
        for i2 in range(len(time)*len(conditions)):
        
            count += 1
        
            axs[i].boxplot(norm_data_area2_summarized[count], positions=[pos[i2]], showfliers=False, showmeans=True, manage_ticks=False)

plot_normalized_intensities2()

#plot_normalized_intensities()

plot_normalized_area2()

#plot_normalized_area()

#plot_normalized_intDen2()

#plot_normalized_intDen()

#scatter_plots()
    
#plot_normalized_intDen()
 
#plot_normalized_intDen_no_wt()

#plot_boxplots_intDen()                          

#plot_normalized_area_no_wt()

#plot_normalized_intensities_no_wt()

#plot_normalized_area()

#plot_normalized_intensities()
            
plot_boxplots_intensity()               
#plot_violinplots_intensity()
plot_boxplots_area()               
#plot_violinplots_area()
            
#plot_dataset_boxplots_intensity(to_plot="2")

print("Finished")