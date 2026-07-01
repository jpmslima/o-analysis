import pandas as pd

# Load the csv
df = pd.read_csv('fubar-results.csv')

# The columns are 'Site', 'α', 'β', 'α-β', 'Prob[α<β]'
# Calculate omega = beta / alpha
# We should handle potential division by zero
df['omega'] = df['β'] / df['α']

# Save back to csv
df.to_csv('fubar-results.csv', index=False)
print("Added omega column successfully.")
