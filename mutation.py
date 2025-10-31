import matplotlib.pyplot as plt 
CODON_TABLE = {
    'ATA':'I','ATC':'I','ATT':'I','ATG':'M',
    'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'AAC':'N','AAT':'N','AAA':'K','AAG':'K',
    'AGC':'S','AGT':'S','AGA':'R','AGG':'R',
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L',
    'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
    'CGA':'R','CGC':'R','CGG':'R','CGT':'R',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
    'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'GAC':'D','GAT':'D','GAA':'E','GAG':'E',
    'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
    'TTC':'F','TTT':'F','TTA':'L','TTG':'L',
    'TAC':'Y','TAT':'Y','TAA':'Stop','TAG':'Stop',
    'TGC':'C','TGT':'C','TGA':'Stop','TGG':'W'
}

def translate_codon(codon):
    codon = codon.upper()
    if len(codon) != 3:
        return "invalid codon"
    return CODON_TABLE.get(codon, 'X') 

def count_syn_nonsyn(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("SEQUENCE MUST BE OF EQUAL LENGTH")
    
    if len(s1) % 3 !=0:
        raise ValueError("SEQUENCE LENGTH MUST BE A MULTIPLE OF 3")
    
    syn = 0
    nonsyn = 0
    changes = []
    
    for i in range(0, len(s1), 3):
        cod1 = s1[i:i+3]
        cod2 = s2[i:i+3]
        if cod1 == cod2:
            continue 
        aa1 = translate_codon(cod1)
        aa2 = translate_codon(cod2)
        if aa1 == aa2:
            syn += 1 
            changes.append((i//3 + 1, cod1, cod2, aa1, aa2, "Synonymous"))
        else:
            nonsyn += 1
            changes.append((i//3 + 1, cod1, cod2, aa1, aa2, "non-synonymous"))
            
    return syn, nonsyn, changes 

s1 = input("Enter the DNA sequence 1:").upper()
s2 = input("Enter the DNA sequence 2:").upper()

syn_count, nonsyn_count, change_list = count_syn_nonsyn(s1, s2)

print("Synonymous changes: ", syn_count)
print("Non-synonymous  changes: ", nonsyn_count)
print("Details of  changes: ")
for idx, c1, c2, aa1, aa2, kind in change_list:
        print(f"Codon #{idx}: {c1} → {c2} ; {aa1} → {aa2} : {kind}")
        
#BAR CHART(SYNONYMOUS AND NON-SYNONYMOUS COUNTS)    
labels = ['Synonymous', 'Non-synonymous']
counts = [syn_count, nonsyn_count]

plt.bar(labels, counts, color = ['blue', 'red'])
plt.title('SYNONYMOUS vs NON-SYNONYMOUS MUTATIONS')
plt.ylabel('Count')
plt.xlabel('Mutation type')
plt.show()




        
        
        
    
        
        
        