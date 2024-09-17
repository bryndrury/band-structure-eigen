import matplotlib.pyplot as plt
import pandas as pd
import time
from matplotlib.animation import FuncAnimation

def plot_data():
    df = pd.read_csv('eigenvalues.txt', sep=' ', header=None)
    bands = 6

    # First column is the k-values
    kvals = df[0]
    # All other columns are the eigenvalues
    eigenvalues = df.drop(0, axis=1)
    energies = eigenvalues.values
    
    plt.clf()

    # Plot the eigenvalues
    for n in range(bands):
        plt.scatter(kvals, energies[:,n], c="black", s=100/energies.shape[0])
        
    plt.hlines(11.6, min(kvals), max(kvals), colors='r', linestyles='dashed')

    plt.draw()
    
def update(event):
    plot_data()
    
    
time.sleep(2)
plot_data()

ani = FuncAnimation(plt.gcf(), update, interval=500)

plt.show()