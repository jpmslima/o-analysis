import re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def main():
    # 1. Parse Clades.txt
    model_to_clade = {}
    current_clade = None
    with open('Clades.txt', 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('Clade'):
                # Store the clade name without the trailing colon
                current_clade = line.rstrip(':')
            else:
                parts = line.split()
                if len(parts) > 0:
                    model_id = parts[0]
                    model_to_clade[model_id] = current_clade

    # 2. Parse Cavities.csv to get the number of cavities for each model
    model_cavities = {}
    with open('Cavities.csv', 'r') as f:
        for line in f:
            line = line.strip()
            # Match lines like "9 cavities found for Pan troglodytes.cif #1"
            match = re.match(r'^(\d+)\s+cavities found for .*?#(\d+)$', line)
            if match:
                num_cavities = int(match.group(1))
                model_id = match.group(2)
                model_cavities[model_id] = num_cavities

    # 3. Join data
    data = []
    for model_id, num_cavities in model_cavities.items():
        clade = model_to_clade.get(model_id, "Unknown Clade")
        data.append({'Model_ID': model_id, 'Clade': clade, 'Cavities': num_cavities})

    df = pd.DataFrame(data)

    # Custom sort by Clade number
    def get_clade_num(c):
        match = re.search(r'Clade (\d+)', c)
        if match:
            return int(match.group(1))
        return 999
    
    df['Clade_Num'] = df['Clade'].apply(get_clade_num)
    df = df.sort_values('Clade_Num')

    # Get ordered unique clades
    unique_clades = df['Clade'].unique()
    
    # 4. Plotting
    plt.figure(figsize=(12, 7))
    
    # Group data for matplotlib boxplot
    clade_data = [df[df['Clade'] == clade]['Cavities'].values for clade in unique_clades]
    
    # Create the boxplot
    bp = plt.boxplot(clade_data, labels=unique_clades, patch_artist=True, 
                     medianprops=dict(color='black', linewidth=1.5),
                     boxprops=dict(linewidth=1.2),
                     whiskerprops=dict(linewidth=1.2),
                     capprops=dict(linewidth=1.2))
    
    # Apply colors
    colors = plt.cm.Pastel1(np.linspace(0, 1, len(unique_clades)))
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
        
    # Overlay actual points (swarm plot equivalent)
    for i, clade in enumerate(unique_clades):
        y = df[df['Clade'] == clade]['Cavities'].values
        # Add some jitter to x-axis so points don't overlap as much
        x = np.random.normal(i + 1, 0.05, size=len(y))
        plt.scatter(x, y, alpha=0.8, color='dimgrey', edgecolor='white', linewidth=0.5, s=50, zorder=3)

    plt.title('Distribution of Cavities Number across Clades', fontsize=16)
    plt.xlabel('Clade', fontsize=14)
    plt.ylabel('Number of Cavities', fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=11)
    plt.yticks(fontsize=11)
    plt.grid(axis='y', linestyle='--', alpha=0.7, zorder=0)
    plt.tight_layout()

    # Save
    plt.savefig('cavities_boxplot.png', dpi=300)
    plt.savefig('cavities_boxplot.pdf')
    print("Boxplot generated and saved as cavities_boxplot.png and cavities_boxplot.pdf")

if __name__ == "__main__":
    main()
