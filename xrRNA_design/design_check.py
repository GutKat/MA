import infrared as ir
import infrared.rna as rna
import RNA
import utils as ut
from collections import Counter
import csv
import random
from datetime import datetime



def remove_gaps(s, gaps):
    s = list(s)
    for gap in gaps:
        s[gap] = '-'
    
    return ''.join(s).replace('-', '')


check = True
target_gc = 0.5
target_energy = -30
target_len = 100


base_ss =   '((((((..((((((((......(((((....))))).(((((...............))))).))))))))..))))))............................'
pk1 =       '...................(((.........................................................................))).........'
pk2 =       '.................................................((((...............................................))))...'
iupac_seq = "NNNNXXXXNNNNXXXGGCAGCRCNNXXNNXXXXNNGACNNXXNNNNXXXGGUCNNXXXXNNGACXXXNNNNXXXXNNNNNNNNNNNNNNXXXXXXUGYXXGACCXXX"


csv_protocol = '/scr/aldea/kgutenbrunner/xrRNA_design/TBFV_design/seqs/designed_seqs.csv'
add_to_csv = True
seq = 'GGCUUC--GCAUCGUGGCAGCACCGGCUUGAGCCGGACCCGCUUUU---GGUCUU--GCGGGACACGGUGC--GAGGCCAUCCUUUUUUUUACUUUGCACGACC---' 
model = 'mc, 100 000'
date = datetime.now().strftime("%m/%d/%Y")
print(date)

seq_ = seq
gaps = [i for i in range(len(seq)) if seq[i] == '-']

base_ss = remove_gaps(base_ss, gaps)
pk1 = remove_gaps(pk1, gaps)
pk2 = remove_gaps(pk2, gaps)

seq = seq.replace('-', '').strip()

fc = RNA.fold_compound(seq)
fc.pf()
freq = fc.pr_structure(base_ss)

fc = RNA.fold_compound(seq)
(ss, mfe) = fc.mfe()
print('suboptimal structures:')
for s in fc.subopt(100):
    print(f"{s.structure} {s.energy:6.2f}")

if ss != base_ss:
    check = False

print(f'\nstructures:\n\t{base_ss}\n\t{pk1}\n\t{pk2}\n')

gc_content = (seq.count("G") + seq.count("C")) / len(seq)
print(f'\nseq:\t{seq}\nss:\t{ss}\nlength:\t{len(seq)}\nMFE:\t{round(mfe,3)}\nfreq:\t{round(freq, 3)}\ngc:\t{round(gc_content, 3)}\nfolds in structure: {check}\n')

if add_to_csv:
    with open(csv_protocol) as f:
        seq_ids = [row.split(';')[0] for row in f]
    seq_id = int(seq_ids[-1]) + 1
    with open(csv_protocol, 'a', newline='\n') as csvfile:
        row_writer = csv.writer(csvfile, delimiter=';',
                                quotechar='"')
        row_writer.writerow([seq_id, seq_, ss, len(seq), round(mfe,4), round(freq, 3), round(gc_content, 3), model, date])


# #check when we add 10 random nts (or just A's) to 5' and 3' whether xrRNA structure remains the same
# nts = 15
# long_seq = ''.join(random.choices(['A', "C", 'G', 'U'], k = nts)) + seq +''.join(random.choices(['A', "C", 'G', 'U'], k = nts))
# long_seq = 'A' * 10 + seq + 'A' * 10 

# print(f'seqs with added nts:\t{long_seq}')
# fc = RNA.fold_compound(long_seq)
# (long_ss, mfe) = fc.mfe()

# print(f'new ss:\t\t\t{long_ss}')
# print(f'ss w\o add nts:\t\t{long_ss[10:-10]}')
# print(f'r structures the same?\t{long_ss[10:-10] == ss}')
