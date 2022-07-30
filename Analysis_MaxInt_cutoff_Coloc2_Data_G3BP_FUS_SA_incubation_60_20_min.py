import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import filedialog
from scipy.stats import ttest_ind
from scipy.stats import kstest
from scipy.stats import median_test

# File chooser dialog
# Please select the csv file containing the coloc2 measurements
root = tk.Tk()
root.withdraw()
path = filedialog.askopenfilename()

df = pd.read_csv(path) # create dataframe from csv file in path

coloc2_stats = df.values # returns values out of df in array

# distance (in indices) to important variables, to calculate index to append values in lists
d_p1 = 0 # distance (indices) to Pearson's R coefficient with no threshold
d_name = 0 # distance to the name of the individual coloc2 analysis
d_MaxInt = 0 # distance to the channel 1 max intensity

# Empty list(s) to append the needed data for analysis, specially Pearson's R
# data = [[name],[cell_type], [condition],[Pearson's no thresh.], [Pearson's below], [Pearson's above]]
data = [[],[],[],[],[],[]]

cutoff1 = 0 # intDen cutoff to determine if a cell expressed enough FUS to be analysed
cutoff2 = 0 # intDen cutoff for antibody labelled G3BP

# loop to find d_p1 index
for distance in range(len(coloc2_stats)):    
    if "Pearson's R value (no threshold)" in coloc2_stats[distance]:
        d_p1 = distance
        break
# loop to find d_name index
for distance in range(len(coloc2_stats)):     
    if "Coloc_Job_Name" in coloc2_stats[distance]:
        d_name = distance
        break
# loop to find d_MaxInt index
for distance in range(len(coloc2_stats)):     
    if "Channel 1 Max" in coloc2_stats[distance]:
        d_MaxInt = distance
        break
            
# iterate over csv file, sort info into data list
for i in range(len(coloc2_stats)):
        
    if "Coloc_Job_Name" in coloc2_stats[i]:
        
        channel1_MaxInt = float(coloc2_stats[i+d_MaxInt-d_name, 1])
        channel2_MaxInt = float(coloc2_stats[i+d_MaxInt-d_name+1, 1])
        
        if channel1_MaxInt > cutoff1 and channel2_MaxInt > cutoff2:
        
            # append name
            data[0].append(coloc2_stats[i,1])
        
            name = coloc2_stats[i,1]
            
            # append Pearson's R (no threshold)
            data[3].append(float(coloc2_stats[i+d_p1-d_name, 1]))
            # append Pearson's R (below threshold)
            data[4].append(float(coloc2_stats[i+d_p1-d_name+1, 1]))
            # append Pearson's R (above threshold)
            data[5].append(float(coloc2_stats[i+d_p1-d_name+2, 1]))
        
            # append cell_type
            if "wt" in name or "WT" in name:
            
                data[1].append("WT")
            
            elif "514" in name:
            
                data[1].append("R514S")
            
            elif "521" in name:
            
               data[1].append("R521C")
        
            # append condition
            if "60min" in name:
            
                data[2].append("60min")
            
            elif "20min" in name:
            
                data[2].append("20min")
                
            elif "no_SA" in name:
                
                data[2].append("Control")
        
# create dataframe from data list, transposing it beforehand to fit formating
df_data0 = pd.DataFrame(np.array(data).T.tolist())
df_data = df_data0.replace("nan", 0)

labels_p1 = [] # Labels of the conditions
mean_p1 = [] # Mean Pearson's R (no threshold) for all conditions
std_p1 = [] # std Dev for Pearson's R (no threshold) for all conditions
    
    
def data_stats_calculation():
    
    mutations = ["WT", "R514S", "R521C"]
    conditions = ["Control","20min", "60min"]
    
    print("--"*20)
    
    # iterate over mutations and conditions to calculate mean and std Dev 
    for mutant in mutations:   
        
        a = df_data[(df_data[1] == mutant)]  
        
        for condition in conditions:
            
            b = a[a[2] == condition].values 
            
            # calculate Pearson's mean (b[0:,3]), but 
            # change the values from str to float beforehand
            mean = np.mean(list(map(float, b[0:,3])))
            std = np.std(list(map(float, b[0:,3])))
            
            mean_p1.append(mean)
            std_p1.append(std)
            
            labels_p1.append(mutant +" "+ condition)
            
            print("Number of cells in " + mutant + " " + condition + ":")
            print(len(b))
            
    print("--"*20)
    
data_stats_calculation()


def print_variables():
    
    print("--"*20)
    print("Mean Pearson's R (no threshold):")
    print(mean_p1)
    print("Std Deviation of Pearson's R (no threshold):")
    print(std_p1)
    print("--"*20)
 
def stat_test(test, alpha):
    
    # function to do statistical test (t-test, mood's-test or KS-test)
    # input: alpha (e.g. 0.05) and test="t" for t-test, test="mood" for mood's test or test="KS" for KS-test
    
    mutations = ["WT", "R514S", "R521C"]
    conditions = ["Control","20min", "60min"]

    print("--"*20)
    print("Statistical Test: " + test + "-Test")
    print("--"*20)
        
    for mutant in mutations: 

        a = df_data[(df_data[1] == mutant)]
        
        for condition in conditions: 

            b = a[(a[2] == condition)]
                
            for mutant2 in mutations:
                    
                c = df_data[(df_data[1] == mutant2)]
                
                if mutant == mutant2:
                    
                    for condition2 in conditions: 
                        
                        d = c[(c[2] == condition2)]
                        
                        if test == "t": 
                            
                            res = ttest_ind(list(map(float, b[3].values)), list(map(float, d[3].values)), equal_var=False).pvalue
                            
                        elif test == "KS":
                            
                            res = kstest(list(map(float, b[3].values)), list(map(float, d[3].values))).pvalue
                            
                        elif test == "mood": 
                                            
                            res = median_test(list(map(float, b[3].values)), list(map(float, d[3].values)))[1]
                            
                        if res < alpha: 
                            
                            print(mutant + " " + condition + " vs " + mutant2 + " " + condition2 + " " + test + "-Test")
                            print(res)
                                
                elif mutant != mutant2:
                    continue

    print("--"*20)
    
stat_test(test="t", alpha=0.05)

def plot_all_pearsons_together():
    
    fig, ax = plt.subplots()
    fig.suptitle("All Pearson's R Values")
    ax.bar(range(len(data[0])), data[3])
    ax.set_ylabel("Pearson's R")
    ax.set_ylim(-1,1)
    #ax.yaxis.grid(True)
    plt.show()

def plot_mean_pearson_individual():    
    
    fig, axs = plt.subplots(1, 3, sharey=True)
    fig.suptitle("Mean Pearson's R Values (No Threshold)")
    axs[0].bar(labels_p1[0:2], mean_p1[0:2], yerr=std_p1[0:2], capsize=8, color = ["b", "r"])
    axs[0].set_ylabel("Mean Pearson's R (No Threshold)")
    axs[0].set_ylim(-1,1)
    axs[1].bar(labels_p1[2:4], mean_p1[2:4], yerr=std_p1[2:4], capsize=8, color = ["b", "r"])
    axs[2].bar(labels_p1[4:6], mean_p1[4:6], yerr=std_p1[4:6], capsize=8, color = ["b", "r"]) 
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=70)
        axs[x].yaxis.grid(True)
    plt.show()
    
def plot_bar_scatter():
    
    mutations = ["WT", "R514S", "R521C"]
    conditions = ["Control", "20min", "60min"]
    c = [] # empty list for later use in plotting individual datapoints over bars
    
    position = [1,2,3] # position for bars in bar diagram
    w = 0.5 # width
    
    # iterate over mutants and conditions, to select Pearson's values corresponding to the mutant+condition
    # this is done to plot individual datapoints over the bars
    for mutant in mutations:   
        
        a = df_data[(df_data[1] == mutant)]  
        
        for condition in conditions:
            
            b = a[a[2] == condition].values
            c.append(b[0:,3])

    
    fig, axs = plt.subplots(1, 3, sharey=True)
    fig.suptitle("Mean Pearson's R Values (No Threshold)")
    axs[0].set_ylabel("Mean Pearson's R (No Threshold)")
    axs[0].set_ylim(-1,1)
    axs[0].bar(position, mean_p1[0:3], yerr=std_p1[0:3], width=w, capsize=8, color=(0,0,0,0), tick_label=labels_p1[0:3], edgecolor=["b", "orange", "r"] )
    axs[0].scatter(position[0] + np.random.random(c[0].size) * w - w / 2, np.array(list(map(float, c[0]))), s=2, color="b")
    axs[0].scatter(position[1] + np.random.random(c[1].size) * w - w / 2, np.array(list(map(float, c[1]))), s=2, color="orange")
    axs[0].scatter(position[2] + np.random.random(c[2].size) * w - w / 2, np.array(list(map(float, c[2]))), s=2, color="r")
    axs[1].bar(position, mean_p1[3:6], yerr=std_p1[3:6], width=w, capsize=8, color=(0,0,0,0), tick_label=labels_p1[3:6], edgecolor=["b", "orange","r"])
    axs[1].scatter(position[0] + np.random.random(c[3].size) * w - w / 2, np.array(list(map(float, c[3]))), s=2, color="b")
    axs[1].scatter(position[1] + np.random.random(c[4].size) * w - w / 2, np.array(list(map(float, c[4]))), s=2, color="orange")
    axs[1].scatter(position[2] + np.random.random(c[5].size) * w - w / 2, np.array(list(map(float, c[5]))), s=2, color="r")
    axs[2].bar(position, mean_p1[6:9], yerr=std_p1[6:9], width=w, capsize=8, color=(0,0,0,0), tick_label=labels_p1[6:9], edgecolor=["b", "orange","r"]) 
    axs[2].scatter(position[0] + np.random.random(c[6].size) * w - w / 2, np.array(list(map(float, c[6]))), s=2, color="b")
    axs[2].scatter(position[1] + np.random.random(c[7].size) * w - w / 2, np.array(list(map(float, c[7]))), s=2, color="orange")
    axs[2].scatter(position[2] + np.random.random(c[8].size) * w - w / 2, np.array(list(map(float, c[8]))), s=2, color="r")
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=70)
        axs[x].yaxis.grid(True)
    plt.show()

def plot_violinplots():
    
    pos = [1,2,3]
    count = -1
    
    mutations = ["WT", "R514S", "R521C"]
    conditions = ["Control", "20min", "60min"]
    
    fig, axs = plt.subplots(1,3, sharey=True)
    fig.suptitle("Pearson's R Values Distribution (No Threshold)")
    axs[0].set_ylabel("Pearson's R (No Threshold)")
    axs[0].set_ylim(-1.2,1.5)
        
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=70)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels_p1[0+x*3:4+x*3])
    
    for mutant in mutations: 
        
        count += 1
        count2 = -1
        
        a = df_data[(df_data[1] == mutant)]
        
        for condition in conditions: 
            
            count2 += 1
            
            b = a[(a[2] == condition)]
            
            axs[count].violinplot(list(map(float, b[3].values)), positions=[pos[count2]], showmeans=True)

    
def plot_boxplots():
    
    pos = [1,2,3]
    count = -1
    
    mutations = ["WT", "R514S", "R521C"]
    conditions = ["Control", "20min", "60min"]
    
    fig, axs = plt.subplots (1,3, sharey=True)    
    fig.suptitle("Mean Pearson's R Values (No Threshold)")
    axs[0].set_ylabel("Mean Pearson's R (No Threshold)")
    axs[0].set_ylim(-1.2,1.5)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=70)
        axs[x].set_xticklabels(labels_p1[0+x*3:4+x*3])
    
    for mutant in mutations: 
        
        count += 1
        count2 = -1
        
        a = df_data[(df_data[1] == mutant)]
        
        for condition in conditions: 
            
            count2 += 1
            
            b = a[(a[2] == condition)]
            
            axs[count].boxplot(list(map(float, b[3].values)), positions=[pos[count2]], showfliers=False, showmeans=True)

    
print_variables()
#plot_all_pearsons_together()
#plot_mean_pearson_individual()
#plot_bar_scatter()
plot_violinplots()
#plot_boxplots()

print("Finished")