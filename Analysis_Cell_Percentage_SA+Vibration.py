import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import filedialog
from scipy.stats import ttest_ind


# File chooser dialog
# Please select the csv file containing the cell measurements
root = tk.Tk()
root.withdraw()
path = filedialog.askopenfilename()

df = pd.read_csv(path) # create dataframe from csv file in path

stats = np.array(df.values) # creates numpy array out of df

# Empty list to append the needed data for analysis
# data = [[name],[cell_type],[condition],[FUS-aggr.],[G3BP-Aggr.],[Both-aggr.],[Dataset]]
data = [[],[],[],[],[],[],[]]

mutants = ["WT", "R514S", "R521C"]
conditions = ["SA Control", "SA + Vib."]
channels = ["sdc488", "sdc640"]
dataset = ["0","1","2"]

def create_data():
    
    name = df["Name"]
    FUS = df["FUS aggr."]
    G3BP = df["G3BP aggr."]
    Both = df["Both aggr."]
    
    for i in range(len(df)):
        
        data[0].append(name[i])
        data[3].append(FUS[i])
        data[4].append(G3BP[i])
        data[5].append(Both[i])

        if "wt" in name[i] or "WT" in name[i]:
            data[1].append("WT")
        elif "514" in name[i]:
            data[1].append("R514S")
        elif "521" in name[i]:
            data[1].append("R521C")
        
        # append condition
        if "RT" in name[i]:
            data[2].append("SA Control")
        elif "1h_Vibration" in name[i]:
            data[2].append("SA + Vib.")
        elif "NO_Vibration" in name[i]:
            data[2].append("SA Control")
        
            
        # append Dataset
        if "210803" in name[i]:
            data[6].append("0")
        elif "Sample_1" in name[i]: 
            data[6].append("1") 
        elif "Sample_2" in name[i]: 
            data[6].append("2")  
            
create_data()

# create dataframe from data list, transposing it beforehand to fit formating
df_data = pd.DataFrame(np.array(data).T.tolist())

total_cells = []
total_FUS = []
total_G3BP = []
total_both = []
labels = []
labels2 = []

def count_all():
    
    # function to count cells with aggregate for each condition, mutant and dataset
    for mutant in mutants: 

        a = df_data[(df_data[1] == mutant)]
        
        for condition in conditions: 

            b = a[(a[2] == condition)]
            
            labels.append(mutant + " " + condition)
            
            for datas in dataset: 
                
                labels2.append(mutant + " " + condition + " Dataset " + datas)
                
                c = b[(b[6] == datas)]
                
                total_cells.append(float(len(c)))
                total_FUS.append(sum(list(map(float, c[3].values))))
                total_G3BP.append(sum(list(map(float, c[4].values))))
                total_both.append(sum(list(map(float, c[5].values))))
    
    #print(labels)        
    #print(total_cells)
    #print(total_FUS)
    #print(total_G3BP)
    #print(total_both)
            
    
count_all()

# individual percentages for each condition
percentages_FUS = []
percentages_G3BP = []
percentages_both = []

# mean percentages
mean_percent_FUS = []
mean_percent_G3BP = []
mean_percent_both = []

# std deviation 
std_FUS = []
std_G3BP = []
std_both = []

def calculate_stats():

        # function to calculate percentages of cells with aggregates
        for i in range(len(total_cells)):
            
            # temporary solution to bad R521C 20 min data (0 cells)
            if total_cells[i] == 0:
               
                x_FUS = 0
                x_G3BP = 0
                x_both = 0
                
                percentages_FUS.append(x_FUS)
                percentages_G3BP.append(x_G3BP)
                percentages_both.append(x_both)
                continue
            
            # calculate and append percentages
            x_FUS = total_FUS[i]/total_cells[i]*100
            x_G3BP = total_G3BP[i]/total_cells[i]*100
            x_both = total_both[i]/total_cells[i]*100
            
            percentages_FUS.append(x_FUS)
            percentages_G3BP.append(x_G3BP)
            percentages_both.append(x_both)
        
        count = -len(dataset)
        
        iterations = int(len(percentages_FUS)/len(dataset))
        
        for x in range(iterations):

            count += len(dataset)
            
            mean_FUS = np.mean(percentages_FUS[count:count+len(dataset)])
            mean_G3BP = np.mean(percentages_G3BP[count:count+len(dataset)])
            mean_both = np.mean(percentages_both[count:count+len(dataset)])
            
            dev_FUS = np.std(percentages_FUS[count:count+len(dataset)])
            dev_G3BP = np.std(percentages_G3BP[count:count+len(dataset)])
            dev_both = np.std(percentages_both[count:count+len(dataset)])
            
            mean_percent_FUS.append(mean_FUS)
            mean_percent_G3BP.append(mean_G3BP)
            mean_percent_both.append(mean_both)
            
            std_FUS.append(dev_FUS)
            std_G3BP.append(dev_G3BP)
            std_both.append(dev_both)
        
calculate_stats()

def print_stats():
    
    print("--"*20)
    print("Stats Printer")
    print("--"*20)
    print("Mean Percent FUS:")
    print(mean_percent_FUS)
    print("Std Deviation:")
    print(std_FUS)
    print("Total amount of cells:")
    print(total_cells)
    
print_stats()

def stat_test(alpha):
    
    # function to do statistical t-test 
    # input: alpha (e.g. 0.05) 

    print("--"*20)
    print("Statistical Test: t-Test")
    print("--"*20)
    
    FUS_data = []
    
    count = -len(conditions)
    iterations = int(len(percentages_FUS)/len(dataset))
    
    for i in range(iterations):
        
        count += len(conditions)
        
        FUS_data.append(percentages_FUS[count:count+len(conditions)])
        
    for i in range(iterations):
        
        x = FUS_data[i]
        
        for i2 in range(len(conditions)):
            
            index2 = i2+i
            
            if index2 >= len(FUS_data):
                
                index2 = len(FUS_data)-1
            
            y = FUS_data[index2]
            
            res = ttest_ind(x, y, equal_var=False).pvalue
            
            if res < alpha:
                
                print(labels[i] + " vs. " + labels[index2])
                print(res)
    
    print("--"*20)
    
stat_test(alpha=0.05)


def plot_bars_FUS():
    
    count = -len(conditions)
    
    position = [1,2] # position for bars in bar diagram
    w = 0.5 # width
    
    fig, axs = plt.subplots(1, 3, sharey=True)
    fig.suptitle("Percentage of Cells with FUS Aggregate")
    axs[0].set_ylabel("Percentage of Cells with Aggregate")
    axs[0].set_ylim(0,100)
    
    for i in range(len(mutants)):
        
        count += len(conditions)

        axs[i].bar(position, mean_percent_FUS[count:count+len(conditions)], yerr=std_FUS[count:count+len(conditions)], width=w, capsize=8, tick_label=labels[count:count+len(conditions)], color = ["orange", "r"]) 
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=70)
        axs[x].yaxis.grid(True)
    plt.show()
    

def plot_bars_G3BP():
    
    count = -len(conditions)
    
    position = [1,2] # position for bars in bar diagram
    w = 0.5 # width
    
    fig, axs = plt.subplots(1, 3, sharey=True)
    fig.suptitle("Percentage of Cells with G3BP Aggregate")
    axs[0].set_ylabel("Percentage of Cells with Aggregate")
    axs[0].set_ylim(0,100)
    
    for i in range(len(mutants)):
        
        count += len(conditions)
        
        
        axs[i].bar(position, mean_percent_G3BP[count:count+len(conditions)], yerr=std_G3BP[count:count+len(conditions)], width=w, capsize=8, tick_label=labels[count:count+len(conditions)], color = ["orange", "r"]) 
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=70)
        axs[x].yaxis.grid(True)
    plt.show()


def plot_bars_both():
    
    count = -len(conditions)
    
    position = [1,2] # position for bars in bar diagram
    w = 0.5 # width
    
    fig, axs = plt.subplots(1, 3, sharey=True)
    fig.suptitle("Percentage of Cells with FUS and G3BP Aggregate")
    axs[0].set_ylabel("Percentage of Cells with Aggregate")
    axs[0].set_ylim(0,100)
    
    for i in range(len(mutants)):
        
        count += len(conditions)

        axs[i].bar(position, mean_percent_both[count:count+len(conditions)], yerr=std_both[count:count+len(conditions)], width=w, capsize=8, tick_label=labels[count:count+len(conditions)], color = ["orange", "r"]) 
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=70)
        axs[x].yaxis.grid(True)
    plt.show()


plot_bars_FUS()
plot_bars_G3BP()
plot_bars_both()


print("Finished")