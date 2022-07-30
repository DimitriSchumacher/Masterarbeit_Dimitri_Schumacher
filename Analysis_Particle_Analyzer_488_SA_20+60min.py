import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import filedialog
from scipy.stats import ttest_ind
from scipy.stats import kstest
from scipy.stats import median_test

# File chooser dialog
# Please select the csv file containing the cell measurements
root = tk.Tk()
root.withdraw()
path = filedialog.askopenfilename()

df = pd.read_csv(path, sep="\t") # create dataframe from csv file in path

stats = np.array(df.values) # creates numpy array out of df

# Empty list to append the needed data for analysis
# data = [[name],[cell_type],[condition], [Area], [Mean Intensity]]
data = [[],[],[],[],[]]

mutants = ["WT", "R514S", "R521C"]
conditions = ["Control", "20 min", "60 min"]

def create_data():

    # function to create data for analysis from df dataframe
    
    # dataframe slices 
    name = df["Label"]
    mean = df["Mean"]
    area = df["Area"]
    
    for i in range(len(name)):
        
        data[0].append(name[i])
        data[3].append(area[i])
        data[4].append(mean[i])
        
            
        if "wt" in name[i] or "WT" in name[i]:
            
            data[1].append("WT")
            
        elif "514" in name[i]:
            
            data[1].append("R514S")
            
        elif "521" in name[i]:
            
            data[1].append("R521C")

        if "20min" in name[i]:
            
            data[2].append("20 min")
            
        elif "60min" in name[i]:
            
            data[2].append("60 min")
        
        elif "no_SA" in name[i]:
            
            data[2].append("Control")
    
    
create_data()

# create dataframe from data list, transposing it beforehand to fit formating

df_data0 = pd.DataFrame(np.array(data).T.tolist())

df_data = df_data0.replace("nan", 0)

labels = []

def create_labels():
    
    for mutant in mutants: 
           
        for condition in conditions:
            
            labels.append(mutant+" "+condition)

create_labels()

def stat_test(test, alpha, parameter):
    
    # function to do statistical test (t-test, mood's-test or KS-test)
    # input: alpha (e.g. 0.05) and test="t" for t-test, test="mood" for mood's test or test="KS" for KS-test
    
    print("--"*20)
    print("Statistical Test: " + test + "-Test, " + parameter)
    print("--"*20)
        
    for mutant in mutants: 

        a = df_data[(df_data[1] == mutant)]
        
        for condition in conditions: 

            b = a[(a[2] == condition)]
                
            for mutant2 in mutants:
                    
                c = df_data[(df_data[1] == mutant2)]
                    
                if mutant == mutant2:
                    
                    for condition2 in conditions: 
                        
                        d = c[(c[2] == condition2)]
                        
                        if parameter == "intensity":
                        
                            if test == "t": 
                                
                                res = ttest_ind(list(map(float, b[4].values)), list(map(float, d[4].values)), equal_var=False).pvalue
                                
                            elif test == "KS":
                                
                                res = kstest(list(map(float, b[4].values)), list(map(float, d[4].values))).pvalue
                            
                            elif test == "mood": 
                                                
                                res = median_test(list(map(float, b[4].values)), list(map(float, d[4].values)))[1]
                                
                            if res < alpha: 
                                
                                print(mutant + " " + condition + " vs " + mutant2 + " " + condition2 + " " + test + "-Test")
                                print(res)
                                
                        elif parameter == "area":
                            
                            if test == "t": 
                                
                                res = ttest_ind(list(map(float, b[3].values)), list(map(float, d[3].values)), equal_var=False).pvalue
                                
                            elif test == "KS":
                                
                                res = kstest(list(map(float, b[3].values)), list(map(float, d[3].values))).pvalue
                            
                            elif test == "mood": 
                                                
                                res = median_test(list(map(float, b[3].values)), list(map(float, d[3].values)))[1]
                                
                            if res < alpha: 
                                
                                print(mutant + " " + condition + " vs " + mutant2 + " " + condition2 + " " + test + "-Test, " + parameter)
                                print(res)
                                
                elif mutant != mutant2:
                    continue

    print("--"*20)
    
stat_test(test="t", alpha=0.05, parameter="intensity")

particle_number = []
mean_mean_intensity = []
mean_area = []

def stats_printer(): 
    
    print("--"*20)
    print("Stats Printer:")
    print("--"*20)
    
    for mutant in mutants: 

        a = df_data[(df_data[1] == mutant)]
        
        for condition in conditions: 

            b = a[(a[2] == condition)]
            
            length = len(b)
            intensity_mean = np.mean(list(map(float, b[4].values)))
            area_mean = np.mean(list(map(float, b[3].values)))
            
            particle_number.append(length)
            mean_mean_intensity.append(intensity_mean)
            mean_area.append(area_mean)
    
    print("Number of particles in each condition: ")
    print(particle_number)
    print("Mean mean intensity of the particles for each condition: ")
    print(mean_mean_intensity)
    print("Mean area of the particles for each condition: ")
    print(mean_area)
    
stats_printer()

def plot_boxplots_intensity():
    
    count = -1
    pos = [1,2,3]
    
    fig, axs = plt.subplots(1, 3, sharey=True)
    fig.suptitle("Boxplots: FUS Particle Mean Intensity")
    axs[0].set_ylabel("Mean Intensity [a.u.]")
    axs[0].set_ylim(-1000,25000)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=70)
        axs[x].set_xticklabels(labels[0+x*3:4+x*3])
    
    for mutant in mutants: 
        
        count += 1
        count2 = -1
        
        a = df_data[(df_data[1] == mutant)]
        
        for condition in conditions: 
            
            count2 += 1
            
            b = a[(a[2] == condition)]
            
            axs[count].boxplot(list(map(float, b[4].values)), positions=[pos[count2]], showfliers=False, showmeans=True)

def plot_violinplots_intensity():
    
    count = -1
    pos = [1,2,3]
    
    fig, axs = plt.subplots(1,3, sharey=True)
    fig.suptitle("Violinplots: FUS Particle Mean Intensity")
    axs[0].set_ylabel("Mean Intensity [a.u.]")
    #axs[0].set_ylim(-1000,30000)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=70)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels[0+x*3:4+x*3])
        
    for mutant in mutants: 
        
        count += 1
        count2 = -1
        
        a = df_data[(df_data[1] == mutant)]
        
        for condition in conditions: 
            
            count2 += 1
            
            b = a[(a[2] == condition)]
            
            axs[count].violinplot(list(map(float, b[4].values)), positions=[pos[count2]])

    
def plot_boxplots_area():
    
    count = -1
    pos = [1,2,3]
    
    fig, axs = plt.subplots(1, 3, sharey=True)
    fig.suptitle("Boxplots: FUS Particle Area")
    axs[0].set_ylabel("Area [$\mu m^2$]")
    axs[0].set_ylim(0,10)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=70)
        axs[x].set_xticklabels(labels[0+x*3:4+x*3])
    
    for mutant in mutants: 
        
        count += 1
        count2 = -1
        
        a = df_data[(df_data[1] == mutant)]
        
        for condition in conditions: 
            
            count2 += 1
            
            b = a[(a[2] == condition)]
            
            axs[count].boxplot(list(map(float, b[3].values)), positions=[pos[count2]], showfliers=False, showmeans=True)

def plot_violinplots_area():
    
    count = -1
    pos = [1,2,3]
    
    fig, axs = plt.subplots(1,3, sharey=True)
    fig.suptitle("Violinplots: FUS Particle Area")
    axs[0].set_ylabel("Area [$\mu m^2$]")
    #axs[0].set_ylim(0,5000)
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=70)
        axs[x].set_xticks(pos)
        axs[x].set_xticklabels(labels[0+x*3:4+x*3])
        
    for mutant in mutants: 
        
        count += 1
        count2 = -1
        
        a = df_data[(df_data[1] == mutant)]
        
        for condition in conditions: 
            
            count2 += 1
            
            b = a[(a[2] == condition)]
            
            axs[count].violinplot(list(map(float, b[3].values)), positions=[pos[count2]])
    

plot_boxplots_intensity()  
#plot_violinplots_intensity()
plot_boxplots_area()
#plot_violinplots_area()

print("Finished")