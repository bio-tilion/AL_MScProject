import random
import pandas as pd
import matplotlib.pyplot as plt

def chinese_restaurant_process(num_customers, alpha, seed=False):
    if seed:
        random.seed(42)
        
    if num_customers <= 0:
        return []
    
    # initial tables
    table_assignments = [1]
    next_open_table = 2

    # generate table assignments for all other customers
    for i in range(1, num_customers):
        if random.random() < alpha / (alpha + i):
            # sit at new table
            table_assignments.append(next_open_table)
            next_open_table += 1
        else:
            # sit at an existing table
            # equal weight is given to open tables
            which_table = random.choice(table_assignments)
            table_assignments.append(which_table)
    
    return table_assignments

n = 100
alphas = [10,5,4,3,2,1]

fig, axs = plt.subplots(2, 3, figsize=(9, 6))
axs = axs.flatten()

for i, alpha in enumerate(alphas):
    assignment = chinese_restaurant_process(n, alpha, seed=True)
    assignment = pd.Series(assignment)
    assignment.value_counts().sort_index().plot.bar(
        title=f"n={n}, alpha={alpha}",
        ax=axs[i],
        color="gray",
    )

plt.tight_layout()
plt.savefig("report/src/chinese_restaurant.png")
