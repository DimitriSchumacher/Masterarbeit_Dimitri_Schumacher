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

df = pd.read_csv(path, encoding='latin1') # create dataframe from csv file in path

stats = np.array(df.values) # creates numpy array out of df

# Empty list to append the needed data for analysis
# data = [[name],[cell_type],[condition],[FUS-aggr.],[Dataset],[time]]
data = [[],[],[],[],[],[]]

mutants = ["WT", "R514S", "R521C"]
conditions = ["Control", "10 $\mu$M", "1 mM"]
time = ["24 h", "48 h"]
dataset = ["0","1","2"]

def create_data():
    
    name = df["Name"]
    FUS = df["FUS aggr."]
    
    for i in range(len(df)):
        
        data[0].append(name[i])
        data[3].append(FUS[i])

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
        if "220223" in name[i] or "220224" in name[i]: 
            data[4].append("0")
        elif "Sample_1" in name[i]: 
            data[4].append("1")
        elif "Sample_2" in name[i]:    
            data[4].append("2")   
        
        # append timepoint
        if "24h" in name[i]:
            data[5].append("24 h")
        elif "48h" in name[i]:
            data[5].append("48 h")
            
create_data()

# create dataframe from data list, transposing it beforehand to fit formating
df_data = pd.DataFrame(np.array(data).T.tolist())

total_cells = []
total_FUS = []

labels = []
labels2 = []

def count_all():
    
    # function to count cells with aggregate for each condition, mutant and dataset
    for mutant in mutants: 

        a = df_data[(df_data[1] == mutant)]
        
        for t in time: 

            b = a[(a[5] == t)]
            
            for condition in conditions:
                
                c = b[(b[2] == condition)]
                
                labels.append(mutant + " " + condition + " " + t)
            
                for datas in dataset: 
                
                    labels2.append(mutant + " " + condition + " " + t + " Dataset " + datas)
                
                    d = c[(c[4] == datas)]
                
                    total_cells.append(float(len(d)))
                    total_FUS.append(sum(list(map(float, d[3].values))))
    
    #print(labels)        
    #print(total_cells)
    #print(total_FUS)
            
    
count_all()

# individual percentages for each condition
percentages_FUS = []

# mean percentages
mean_percent_FUS = []

# std deviation 
std_FUS = []

def calculate_stats():

        # function to calculate percentages of cells with aggregates
        for i in range(len(total_cells)):
            
            # calculate and append percentages
            x_FUS = total_FUS[i]/total_cells[i]*100

            percentages_FUS.append(x_FUS)
        
        count = -len(dataset)
        
        iterations = int(len(percentages_FUS)/len(dataset))
        
        for x in range(iterations):

            count += len(dataset)
            
            mean_FUS = np.mean(percentages_FUS[count:count+len(dataset)])
            dev_FUS = np.std(percentages_FUS[count:count+len(dataset)])
            mean_percent_FUS.append(mean_FUS)
            std_FUS.append(dev_FUS)

        #print(mean_percent_FUS)
        
calculate_stats()

def stat_test(alpha):
    
    # function to do statistical t-test 
    # input: alpha (e.g. 0.05) 

    print("--"*20)
    print("Statistical Test: t-Test")
    print("--"*20)
    
    FUS_data = []
    
    count = -len(dataset)*len(time)
    iterations = int(len(percentages_FUS)/len(dataset))
    
    for i in range(iterations):
        
        count += len(dataset)*len(time)
        
        FUS_data.append(percentages_FUS[count:count+len(dataset)*len(time)])
        
    for i1 in range(iterations):
        
        x = FUS_data[i1]
        
        for i2 in range(len(dataset)*len(time)):
            
            index2 = i2+i1
            
            if index2 >= len(FUS_data):
                
                index2 = len(FUS_data)-1
            
            y = FUS_data[index2]
            
            #print([x,y])
            
            res = ttest_ind(x, y, equal_var=False).pvalue
            
            if res < alpha:
                
                print(labels[i1] + " vs. " + labels[index2])
                print(res)
    
    print("--"*20)
    
stat_test(alpha=0.05)

def plot_bars_FUS():
    
    count = -(len(conditions)*len(time))
    
    position = [1,2,3,4,5,6] # position for bars in bar diagram
    w = 0.5 # width
    
    fig, axs = plt.subplots(1, 3, sharey=True)
    fig.suptitle("GGYGG: Percentage of Cells with FUS Aggregate")
    axs[0].set_ylabel("Percentage of Cells with Aggregate")
    axs[0].set_ylim(0,100)
    
    for i in range(len(mutants)):
        
        count += (len(conditions)*len(time))

        axs[i].bar(position, mean_percent_FUS[count:count+(len(conditions)*len(time))], yerr=std_FUS[count:count+(len(conditions)*len(time))], width=w, capsize=8, tick_label=labels[count:count+(len(conditions)*len(time))], color = ["orange", "orange", "orange", "r", "r","r"]) 
    
    for x in range(len(axs)):
        axs[x].tick_params(axis="x", rotation=70)
        axs[x].yaxis.grid(True)
    plt.show()

plot_bars_FUS()



print("Finished")