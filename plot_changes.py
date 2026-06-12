import csv
import matplotlib.pyplot as plt
import numpy as np

def main():
    filename = 'dna.aligned-changes-list.csv'
    
    with open(filename, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        lines = list(reader)
        
    # Find the row with the column headers (taxa names)
    header_row = None
    data_start_row = None
    
    for i, row in enumerate(lines):
        # The taxa line comes after the line with just 'Taxa'
        if len(row) > 1 and any(cell.strip() == 'Marmota monax' for cell in row):
            header_row = i
        if len(row) > 0 and row[0].strip() == 'Site':
            data_start_row = i + 1
            break
            
    if header_row is None or data_start_row is None:
        print("Could not find header or data rows in the CSV.")
        return
        
    headers = lines[header_row]
    
    # Initialize counts for each valid column
    counts = {i: 0 for i in range(1, len(headers)) if headers[i].strip() != ''}
    
    # Count non-empty cells in the data section
    for row in lines[data_start_row:]:
        for i in counts.keys():
            if i < len(row):
                cell = row[i].strip()
                # A change is considered as any value in the cell
                if cell != '':
                    counts[i] += 1
                    
    # Prepare data for plotting
    taxa_counts = {}
    for i, count in counts.items():
        name = headers[i].strip()
        # Filter to include only named taxa (excluding internal numbered nodes like "30")
        if any(c.isalpha() for c in name):
            taxa_counts[name] = count
            
    # Sort taxa by number of changes (descending)
    sorted_taxa = sorted(taxa_counts.items(), key=lambda x: x[1], reverse=True)
    names = [x[0] for x in sorted_taxa]
    values = [x[1] for x in sorted_taxa]
    
    # Plotting
    plt.figure(figsize=(14, 8))
    
    # Create the bar plot
    bars = plt.bar(names, values, color='cornflowerblue', edgecolor='black', linewidth=1)
    
    # Add data labels on top of each bar
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                 f'{int(height)}',
                 ha='center', va='bottom', fontsize=9)
    
    plt.title('Total Number of Changes per Taxa', fontsize=16)
    plt.xlabel('Taxa', fontsize=14)
    plt.ylabel('Number of Changes', fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=11)
    
    # Set y limits slightly higher to make room for text labels
    actual_max = max(values) if values else 0
    plt.ylim(0, actual_max * 1.1)
    
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('taxa_changes_plot.png', dpi=300)
    plt.savefig('taxa_changes_plot.pdf')
    
    print("Successfully generated taxa changes plot.")

if __name__ == '__main__':
    main()
