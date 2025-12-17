import sys
import os


INPUT_PDB = "6m0j_fixed.pdb"  
CHAIN_TO_MUTATE = "B"         
RESIDUE_NUMBER = 486         
OUTPUT_FILENAME = f"mutant_{CHAIN_TO_MUTATE}{RESIDUE_NUMBER}_ALA.pdb"

def create_alanine_mutant(input_file, output_file, chain_id, res_num):
    print(f"Creating mutant: Chain {chain_id}, Residue {res_num} -> ALA...")
    
   
    allowed_atoms = [" N  ", " CA ", " C  ", " O  ", " CB ", " OXT"]
    
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                line_chain = line[21]
                try:
                    line_res = int(line[22:26])
                except ValueError:
                    f_out.write(line)
                    continue

         
                if line_chain == chain_id and line_res == res_num:
                    atom_name = line[12:16]
                    
                    if atom_name in allowed_atoms:
                        
                        new_line = line[:17] + "ALA " + line[20:]
                        f_out.write(new_line)
                    else:
                        
                        continue
                else:
                    
                    f_out.write(line)
            else:
                f_out.write(line)
                
    print(f"Success! File created: {output_file}")

# Executar
if __name__ == "__main__":
    if not os.path.exists(INPUT_PDB):
        print(f"Error: Can't find {INPUT_PDB}")
        # Tenta procurar na pasta pai se não estiver aqui
        if os.path.exists("../" + INPUT_PDB):
            INPUT_PDB = "../" + INPUT_PDB
            print(f"Found on previous folder: {INPUT_PDB}")
        else:
            sys.exit(1)

    create_alanine_mutant(INPUT_PDB, OUTPUT_FILENAME, CHAIN_TO_MUTATE, RESIDUE_NUMBER)