import infrared as ir
import infrared.rna as rna
import RNA
import math
import utils as ut
from collections import Counter
import csv
import random
from datetime import datetime

def add_nt(seq, n_nts):
    return ''.join(random.choices(['A', "C", 'G', 'U'], k=n_nts)) + seq + ''.join(random.choices(['A', "C", 'G', 'U'], k=n_nts))

def remove_gaps(s, gaps):
    s = list(s)
    for gap in gaps:
        s[gap] = '-'
    return ''.join(s).replace('-', '')

def margin_left(left_text, right_text, padding):
    print(f"{left_text: <{padding}}{right_text}")

def print_line():
    print(f"\n{'':->150}\n")


target_len = [52, 70]
target_gc =  0.58
target_energy = -30


iupac_cons =    'WGUCAGSCCXXXXNNNXXXXXXXXGMYRCNXXXXXXXXXXXNNNXXXXXNXXXXXXXXNGUGCWGBCUGXXXXXXXXXXXNNNNN'
base_ss =       '...((((((((((.......))))).((((((((....................))))))))..)))))................'
pk1 = 		    '((............................................................)).....................'
pk2 = 		    '.........................................((((((((..........................))))))))..'
seq =           'AGUCAGGCCGUC-UUAU----GACGCCACCCCG-GC-----CCA-----G-----CGGGGUGCUGCCUG-----------UGGAA'
model = 'mc, 100 000,  cosntraints for beta and gamma, opt = target frequency - (ensemble defect + sigmoid(pk2_energy))'
date = datetime.now().strftime("%m/%d/%Y")

csv_protocol = '/scr/aldea/kgutenbrunner/working/xrRNA_design/MBFV_design/seqs/designed_seqs.csv'
add_to_csv = True

seq_ = seq
gaps = [i for i in range(len(seq)) if seq[i] == '-']

base_ss = remove_gaps(base_ss, gaps)
pk1 = remove_gaps(pk1, gaps)
pk2 = remove_gaps(pk2, gaps)

seq = seq.replace('-', '').strip()
fc = RNA.fold_compound(seq)
fc.pf()
freq = fc.pr_structure(base_ss)
ef = fc.ensemble_defect(base_ss)
pk2_e = fc.eval_structure(pk2)
obj_f = freq - (math.tanh(pk2_e*0.25) * 0.1)

(ss, mfe) = fc.mfe()
print('suboptimal structures:')
for s in fc.subopt(100):
    print(f"\t{s.structure}\t{s.energy:6.2f}")


print(f'\nstructures:\n{seq}\n{base_ss}\n{pk1}\n{pk2}\n')

gc_content = (seq.count("G") + seq.count("C")) / len(seq)
margin_left('seq:', seq, 15)
margin_left('ss:', ss, 15)
margin_left('length:', len(seq), 15)
margin_left('MFE:', round(mfe,3), 15)
margin_left('PK2 energy:', round(pk2_e,3), 15)
margin_left('freq:', round(freq, 3), 15)
margin_left('ensemble def:', round(ef, 3), 15)
margin_left('gc:', round(gc_content, 3), 15)
print_line()

if add_to_csv:
    with open(csv_protocol) as f:
        seq_ids = [row.split(';')[0] for row in f]
    seq_id = int(seq_ids[-1]) + 1
    with open(csv_protocol, 'a', newline='\n') as csvfile:
        row_writer = csv.writer(csvfile, delimiter=';',
                                quotechar='"')
        row_writer.writerow([seq_id, seq_, ss, len(seq), round(mfe,4), round(freq, 3), round(ef, 3), round(gc_content, 3), model, date])


# #check when we add 10 random nts (or just A's) to 5' and 3' whether xrRNA structure remains the same
nts = 10
long_seq = ''.join(random.choices(['A', "C", 'G', 'U'], k = nts)) + seq +''.join(random.choices(['A', "C", 'G', 'U'], k = nts))
#long_seq = 'A' * 10 + seq + 'A' * 10
long_seq_p = ''.join(' ' * nts) +  ''.join('*' * len(seq)) + ''.join(' ' * nts)


fc = RNA.fold_compound(long_seq)
(long_ss, mfe) = fc.mfe()
print('Adding nts up and downstream of designed sequence\n')
margin_left('designed seqs:', long_seq_p, 30)
margin_left('seqs with added nts:', long_seq, 30)
margin_left('new ss:', long_ss, 30)
margin_left('ss w\\o added nts:', ''.join(' ' * nts) + long_ss[nts:-nts] + ''.join(' ' * nts), 30)
margin_left('True structure:', ''.join(' ' * nts) + ss + ''.join(' ' * nts), 30)
