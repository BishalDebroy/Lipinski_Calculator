import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, PandasTools
import pandas as pd

def main():
    if len(sys.argv) != 3:
        print("Usage: python calculate_lipinski.py input.smi output.csv")
        return
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Read SMILES file
    df = pd.read_csv(input_file, sep='\t', header=None, names=['SMILES', 'Name'])
    
    # Calculate properties
    results = []
    for smi, name in zip(df['SMILES'], df['Name']):
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                hbd = Descriptors.NumHDonors(mol)
                hba = Descriptors.NumHAcceptors(mol)
                
                # Check Lipinski rules
                rule1 = mw <= 500
                rule2 = logp <= 5
                rule3 = hbd <= 5
                rule4 = hba <= 10
                passes = rule1 and rule2 and rule3 and rule4
                
                results.append({
                    'Name': name,
                    'SMILES': smi,
                    'MW': round(mw, 2),
                    'LogP': round(logp, 2),
                    'HBD': hbd,
                    'HBA': hba,
                    'Passes_Lipinski': passes
                })
        except:
            continue
    
    # Create results DataFrame
    results_df = pd.DataFrame(results)
    
    # Save to CSV
    results_df.to_csv(output_file, index=False)
    print(f"Processed {len(results_df)} molecules. Results saved to {output_file}")

if __name__ == "__main__":
    main()