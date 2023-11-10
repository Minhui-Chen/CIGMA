import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Create a 5x5 matrix of random numbers between 0 and 1
matrix = np.random.rand(10, 5)

# Define the size and color of the circles
circle_size = 0.2
#colors = np.random.rand(matrix.shape[0], matrix.shape[1]) 
#colors = 0.9 * (colors - np.min(colors)) / (np.max(colors) - np.min(colors)) + 0.1
colors = 1 - np.array([
    [0.6, 0.2, 0.6, 0.5, 0.5],
    [0.7, 0.1, 0.5, 0.4, 0.4],
    [0.1, 0.3, 0.4, 0.2, 0.1],
    [0.5, 0.1, 0.3, 0.3, 0.1],
    [0.3, 0.1, 0.4, 0.1, 0.1],
    [0.3, 0.6, 0.1, 0.4, 0.4],
    [0.1, 0.1, 0.1, 0.1, 0.3],
    [0.3, 0.7, 0.5, 0.4, 0.1],
    [0.2, 0.6, 0.1, 0.1, 0.6],
    [0.3, 0.5, 0.6, 0.7, 0.3]
    ])

# Create a figure and axes
fig, ax = plt.subplots(dpi=600)

# Loop through each element of the matrix and plot a circle with fixed size and random color
for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
        color = colors[i, j]
        if i in [0,1]:
            edge_color = 'red'
        elif i in [2,3,4]:
            edge_color = 'blue'
        elif i in [5,6,7]:
            edge_color = 'green'
        else:
            edge_color = 'black'
        circle = plt.Circle((j, i), circle_size, color=str(color))
        #circle = plt.Circle((j, i), circle_size, color=str(color), ec=edge_color, lw=0.5)
        ax.add_artist(circle)

# Set the x and y limits and turn off the axis ticks
ax.set_xlim([-0.9, matrix.shape[1] - 0.5])
ax.set_ylim([-0.5, matrix.shape[0] - 0.5])
ax.set_xticks([])
ax.set_yticks([])

# Set the aspect ratio of the figure to "equal"
ax.set_aspect("equal")

# Define the labels and positions for the vertical text
labels = ['CT2', 'CT1', 'CT2', 'CT1']
positions = [(-0.5, 0.5), (-0.5, 3), (-0.5, 6), (-0.5, 8.5)]
inds = ['Ind2', 'Ind1']
ind_positions = [(-0.5, 2.0), (-0.5, 7)]

# Add the vertical text labels
#ax.text( -1 - 0.6, 4.5, 'Cells', fontsize=14, ha='center', va='center', rotation=90, transform=ax.transData)
for i in range(matrix.shape[0]):
    ax.text(-.5, i, f'C{10-i}', fontsize=8, ha='center', va='center', transform=ax.transData)

for i in range(len(labels)):
    ax.text(positions[i][0] - 0.55, positions[i][1], labels[i], fontsize=10, ha='center', va='center', rotation=90, transform=ax.transData)

for i in range(len(inds)):
    ax.text(ind_positions[i][0] - 1.2, ind_positions[i][1], inds[i], fontsize=14, ha='center', va='center', rotation=90, transform=ax.transData)

ax.text( 2, 10.2, 'Genes', fontsize=14, ha='center', va='center', transform=ax.transData)
for j in range(matrix.shape[1]):
    ax.text(j, 9.5, f'G{j+1}', fontsize=8, ha='center', va='center', transform=ax.transData)

# Add the lines to indicate which rows belong to each group
line_x = -.85
ct1_line = Line2D([line_x, line_x], [0, 1], color='red', linewidth=2, transform=ax.transData)
ct2_line = Line2D([line_x, line_x], [2, 4], color='blue', linewidth=2, transform=ax.transData)
ct3_line = Line2D([line_x, line_x], [5, 7], color='red', linewidth=2, transform=ax.transData)
ct4_line = Line2D([line_x, line_x], [8, 9], color='blue', linewidth=2, transform=ax.transData)
ax.add_line(ct1_line)
ax.add_line(ct2_line)
ax.add_line(ct3_line)
ax.add_line(ct4_line)

# Remove the box around the plot
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)

# Save the plot to a file
plt.savefig('circle_matrix.png')
