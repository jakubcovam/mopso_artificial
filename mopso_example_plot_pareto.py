# --------------------------------------------
# Import libraries
# --------------------------------------------
import matplotlib.pyplot as plt
import pandas as pd
# --------------------------------------------

# --------------------------------------------
# Read file, prepare variables, and plot
# --------------------------------------------
# Read data from CSV file from C++ project
data = pd.read_csv('ParetoFront/all_example_solutions.csv')

# Separate Pareto front and dominated solutions
pareto = data[data['Type'] == 'Pareto']
dominated = data[data['Type'] == 'Dominated']

# Create a figure with two subplots side by side
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
# --------------------------------------------

# --------------------------------------------
# Plot Pareto Set in Decision Space (x1 vs. x2)
# --------------------------------------------
# Plot all solutions in decision space (optional)
# axes[0].scatter(data['x1'], data['x2'], label='All Solutions', color='gray', alpha=0.3)

# Plot Pareto front solutions in decision space
axes[0].scatter(pareto['x1'], pareto['x2'], label='Pareto Set', color='red')

# Set labels and title
axes[0].set_xlabel('x1')
axes[0].set_ylabel('x2')
axes[0].set_title('Pareto Set in Decision Space')
axes[0].legend(loc='best')
axes[0].grid(True)
# --------------------------------------------

# --------------------------------------------
# Plot Pareto Front in Objective Space (f1 vs. f2)
# --------------------------------------------
# Plot dominated solutions
axes[1].scatter(dominated['f1'], dominated['f2'], label='Dominated Solutions', color='gray', alpha=0.5)

# Plot Pareto front points
axes[1].scatter(pareto['f1'], pareto['f2'], label='Pareto Front Points', color='red')

# Draw a line through Pareto front solutions
# Sort Pareto front solutions based on one objective (f1)
pareto_sorted = pareto.sort_values(by='f1')
# Plot the Pareto front line with a label for the legend
axes[1].plot(pareto_sorted['f1'], pareto_sorted['f2'], color='blue', linestyle='--', linewidth=1.5, label='Pareto Front Line')

# Set labels and title
axes[1].set_xlabel('f1')
axes[1].set_ylabel('f2')
axes[1].set_title('Pareto Front in Objective Space')
axes[1].legend(loc='best')
axes[1].grid(True)
# --------------------------------------------

# --------------------------------------------
# Adjust layout, save and show
# --------------------------------------------
# Adjust layout to prevent overlap
plt.tight_layout()

# Save the plot as PDF
# plt.savefig('example_mopso.pdf')

# Save the plot as PNG
plt.savefig('example_mopso.png', dpi=300)

# Show the combined plot
plt.show()
# --------------------------------------------
