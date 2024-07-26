import design
import pandas as pd

data = []
steps = [50000, 100000, 250000]
each_step = 5

objective_function = {0: 'frequency', 1:'ensemble_defect', 2:'both'}

for step in steps:
    for key in objective_function.keys():
        for _ in range(each_step):
            print(objective_function[key], step)
            sample, freq, ed, pk2_e = design.testing(step,key)
            data.append([objective_function[key], step, sample, round(freq, 3), round(ed, 3), round(pk2_e, 3)])
            print(objective_function[key], step, sample, round(freq, 3), round(ed, 3), round(pk2_e, 3))

df = pd.DataFrame(data,
columns=['objective_function', 'steps', 'sequence', 'frequency', 'ensemble_defect', 'pk2_energy'])
print(df)
df.to_csv('/scr/aldea/kgutenbrunner/working/xrRNA_design/MBFV_design/testing_obj.csv')