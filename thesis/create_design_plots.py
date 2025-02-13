
import matplotlib.pyplot as plt
# import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from varnaapi import Structure
import varnaapi



def find_structures_TBFV(ss):
    pk1 = []
    pk1_start = False
    pk2 = []
    pk2_start = False
    alpha = []
    alpha_stack = []
    alpha_start = False
    beta = []
    beta_stack = []
    beta_start = False
    gamma_start = False
    gamma = []
    gamma_stack = []
    for i, nt in enumerate(ss):
        if nt == '<' or nt == '>':
            pk1.append(i+1)
            pk1_start = True
        elif nt == '[' or nt == ']':
            pk2.append(i+1)
            pk2_start = True
        elif nt == '(' and not pk1_start and not gamma_start:
            alpha.append(i+1)
        elif nt == '(' and pk1_start:
            beta_stack.append(i+1)
            beta.append(i+1)
        elif nt == ')' and beta_stack:
            beta_stack.pop()
            beta.append(i+1)
            gamma_start = True
            pk1_start = False
        elif nt == '(' and gamma_start:
            gamma.append(i+1)
        elif (nt == '(' and gamma_start) or (nt == ')' and gamma_start):
            gamma.append(i+1)
    n = len(alpha)
    alpha += gamma[-n:]
    gamma = gamma[:-n]
    return alpha, beta, gamma, pk1, pk2
    

def find_structures_MBFV(ss):
    pk1 = []
    pk2 = []
    alpha = []
    beta = []
    beta_stack = []
    gamma_start = False
    gamma = []
    for i, nt in enumerate(ss):
        if nt == '<' or nt == '>':
            pk1.append(i+1)
        elif nt == '[' or nt == ']':
            pk2.append(i+1)
        elif nt == '(' and len(alpha)<5:
            alpha.append(i+1)
        elif nt == '(' and not gamma_start:
            beta_stack.append(i+1)
            beta.append(i+1)
        elif nt == ')' and beta_stack:
            beta_stack.pop()
            beta.append(i+1)
            gamma_start = True
        elif (nt == '(' and gamma_start) or (nt == ')' and gamma_start):
            gamma.append(i+1)
    alpha += gamma[-5:]
    gamma = gamma[:-5]
    return alpha, beta, gamma, pk1, pk2



def plot_color_varnaapi_MBFV(seq, ss, ring):
    alpha, beta, gamma, pk1, pk2 = find_structures_MBFV(ss)
    v = varnaapi.Structure(structure=ss)
    v.update(bpStyle='simple', bp='#333232', drawBackbone=True,  fillBases=True, spaceBetweenBases=1, resolution=1, )

    colors ={
        'alpha':"#fbff00", #fbff00
        'beta_stem':'#ff1303',
        #'beta_loop':'#ffaea8',
        'gamma_stem':'#000dff',
        #'gamma_loop':'#999ef7',
        'ring':'cyan',
        'PK1': '#d272f7', #d272f7
        'PK2':'#ff9e03',
    }

    # alpha
    style = varnaapi.param.BasesStyle(fill=colors['alpha'])
    v.add_bases_style(style, alpha)
    a, b, c, d =  structure_parts_finder(alpha)
    v.add_highlight_region(a, b, radius=10, fill=colors['alpha'], outline=colors['alpha'])
    v.add_highlight_region(c, d, radius=10, fill=colors['alpha'], outline=colors['alpha'])

    # beta stem
    style = varnaapi.param.BasesStyle(fill=colors['beta_stem'])
    a, b, c, d =  structure_parts_finder(beta)
    v.add_bases_style(style, beta)
    v.add_highlight_region(a, d, radius=10, fill=colors['beta_stem'], outline=colors['beta_stem'])
    # # beta loop
    style = varnaapi.param.BasesStyle(fill=colors['beta_loop'])
    bps = [*range(b, c+1)]
    v.add_bases_style(style, bps)
    v.add_highlight_region(b, c, radius=10, fill=colors['beta_loop'], outline=colors['beta_loop'])


    # gamma stem
    style = varnaapi.param.BasesStyle(fill=colors['gamma_stem'])
    a, b, c, d =  structure_parts_finder(gamma)
    v.add_bases_style(style, gamma)
    v.add_highlight_region(a, d, radius=10, fill=colors['gamma_stem'], outline=colors['gamma_stem'])
    # gamma loop
    style = varnaapi.param.BasesStyle(fill=colors['gamma_loop'])
    bps = [*range(b, c+1)]
    v.add_bases_style(style, bps)
    v.add_highlight_region(b, c,cradius=10, fill=colors['gamma_loop'], outline=colors['gamma_loop'])

    #v.add_highlight_region(*ring, radius=10, fill=colors['ring'], outline=colors['ring'])



    # #  PK1
    style = varnaapi.param.BasesStyle(fill=colors['PK1'])
    v.add_bases_style(style, pk1)
    v.add_highlight_region(pk1[0], pk1[1], radius=10, fill=colors['PK1'], outline=colors['PK1'])
    v.add_highlight_region(pk1[2], pk1[3], radius=10, fill=colors['PK1'], outline=colors['PK1'])
    v.add_aux_BP(pk1[0], pk1[3], thickness=3,  color=colors['PK1'])
    v.add_aux_BP(pk1[1], pk1[2], thickness=3,  color=colors['PK1'])

    #  PK2
    style = varnaapi.param.BasesStyle(fill=colors['PK2'])
    a, b, c, d =  structure_parts_finder(pk2)
    v.add_bases_style(style, pk2)
    v.add_highlight_region(c, d, radius=10, fill=colors['PK2'], outline=colors['PK2'])

    for i in range(int(len(pk2)/2)):
        v.add_aux_BP(pk2[i], pk2[-(i+1)], thickness=3, color=colors['PK2'])


    v.update(baseNum='white')
    v.update(periodNum=100)

    v.show()
#v.savefig('/scr/aldea/kgutenbrunner/github/MA/thesis/images/MBFV_part_structures.png')


def plot_color_varnaapi_TBFV(seq, ss, ring):
    alpha, beta, gamma, pk1, pk2 = find_structures_TBFV(ss)

    v = varnaapi.Structure(structure=ss)
    v.update(bpStyle='simple', bp='#333232', drawBackbone=True,  fillBases=True, spaceBetweenBases=1, resolution=1, )

    colors ={
        'alpha':"#fbff00", #fbff00
        'beta_stem':'#ff1303',
        'beta_loop':'#ffaea8',
        'gamma_stem':'#000dff',
        'gamma_loop':'#999ef7',
        'ring':'cyan',
        'pk1': '#d272f7', #d272f7
        'pk2':'#ff9e03',
    }

    # alpha
    style = varnaapi.param.BasesStyle(fill=colors['alpha'])
    v.add_bases_style(style, alpha)
    a, b, c, d =  structure_parts_finder(alpha)
    v.add_bases_style(style, [*range(a,b), *range(c,d)])
    v.add_highlight_region(a, b, radius=10, fill=colors['alpha'], outline=colors['alpha'])
    v.add_highlight_region(c, d, radius=10, fill=colors['alpha'], outline=colors['alpha'])

    # beta stem
    style = varnaapi.param.BasesStyle(fill=colors['beta_stem'])
    a, b, c, d =  structure_parts_finder(beta)
    v.add_bases_style(style, beta)
    v.add_highlight_region(a, d, radius=10, fill=colors['beta_stem'], outline=colors['beta_stem'])


    # # beta loop
    style = varnaapi.param.BasesStyle(fill=colors['beta_loop'])
    bps = [*range(b, c+1)]
    v.add_bases_style(style, bps)
    v.add_highlight_region(b, c, radius=10, fill=colors['beta_loop'], outline=colors['beta_loop'])


    # gamma stem
    style = varnaapi.param.BasesStyle(fill=colors['gamma_stem'])
    a, b, c, d =  structure_parts_finder(gamma)
    v.add_bases_style(style, gamma)
    v.add_highlight_region(a, d, radius=10, fill=colors['gamma_stem'], outline=colors['gamma_stem'])


    # gamma loop
    style = varnaapi.param.BasesStyle(fill=colors['gamma_loop'])
    bps = [*range(b, c+1)]
    v.add_bases_style(style, bps)
    v.add_highlight_region(b, c,cradius=10, fill=colors['gamma_loop'], outline=colors['gamma_loop'])

    #v.add_highlight_region(*ring, radius=10, fill=colors['ring'], outline=colors['ring'])


    # #  PK1
    style = varnaapi.param.BasesStyle(fill=colors['PK1'])
    a, b, c, d =  structure_parts_finder(pk1)
    v.add_bases_style(style, pk1)
    v.add_highlight_region(a, b, radius=10, fill=colors['PK1'], outline=colors['PK1'])
    v.add_highlight_region(c, d, radius=10, fill=colors['PK1'], outline=colors['PK1'])

    for i in range(n):
        v.add_aux_BP(pk1[i], pk1[-(i+1)], thickness=3, color=colors['PK1'])


    #  PK2
    style = varnaapi.param.BasesStyle(fill=colors['PK2'])
    a, b, c, d =  structure_parts_finder(pk2)
    v.add_bases_style(style, pk2)
    v.add_highlight_region(c, d, radius=10, fill=colors['PK2'], outline=colors['PK2'])

    for i in range(int(len(pk2)/2)):
        v.add_aux_BP(pk2[i], pk2[-(i+1)], thickness=3, color=colors['PK2'])


    v.update(baseNum='white')
    v.update(periodNum=100)

    v.show()

def structure_parts_finder(structure):
    n = int(len(structure)/2)
    return structure[0], structure[n-1], structure[n], structure[-1]



def plot_colored_seq(seq, ss, FV, save_plot = None):
    plt.rcParams['text.usetex'] = True

    if FV=='TBFV':
        alpha, beta, gamma, pk1, pk2 = find_structures_TBFV(ss)
        alpha_1, alpha_2, alpha_3, alpha_4 = structure_parts_finder(alpha)
        alpha = [*range(alpha_1, alpha_2+1)] + [*range(alpha_3, alpha_4+1)]

    else:
        alpha, beta, gamma, pk1, pk2 = find_structures_MBFV(ss)

    _, beta_l1, beta_l2, _ =  structure_parts_finder(beta)
    beta_hl = [*range(beta_l1, beta_l2)]
    _, gamma_l1, gamma_l2, _ =  structure_parts_finder(gamma)
    gamma_hl = [*range(gamma_l1, gamma_l2)]
    gamma_hl = [x for x in gamma_hl if x not in pk2]


    colors ={
            'alpha':"#fbff00",
            'beta_stem':'#ff1303',
            #'beta_loop':'#ffaea8',
            'gamma_stem':'#000dff',
            #'gamma_loop':'#999ef7',
            #'ring':'cyan',
            'pk1': '#d272f7', 
            'pk2':'#ff9e03',
            None: 'gray'
            }


    coloring_dic = {'alpha': alpha, 'beta_stem':beta, 'beta_loop':beta_hl,'gamma_stem':gamma, 'gamma_loop':gamma_hl, 'pk1':pk1, 'pk2':pk2}

    def find_color(value):
        for key, values in coloring_dic.items():
            if value in values:
                return colors[key]
        return colors[None]  # Return None if the value is not found


    #plt.rcParams['text.usetex'] = True
    numbers = [str(i) if i%10==0 else  '.' for i in range(1,len(seq)+1)]
    numbers = [numbers[i] if i%10!=8 else '' for i in range(len(numbers))]
    numbers[0] = str(1)
    numbers[-2] = str(len(seq))[-2]
    numbers[-1] = str(len(seq))[-1]
    numbers = ''.join(numbers)
    
    # Nucleotide sequences (for three lines)
    texts = [numbers, ss, seq]
    
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(20, 1))  # Adjust the height to accommodate multiple lines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    line_spacing = 0.35
    letter_spacing = 0.025
    letter_spacing = [0.01, 0.0100, 0.01]
    lettering_minus = [0.1, 0.097,0.1]
    # Plot each nucleotide in each line with a color from the rainbow gradient
    for line_idx, text in enumerate(texts):
        for j, nt in enumerate(text):
            if line_idx == 0:
                color = 'black'
            else:
                color = find_color(j+1)
            ax.text(j*letter_spacing[line_idx]-lettering_minus[line_idx], line_spacing * line_idx, f"\\texttt{{{nt}}}", fontsize=18, color=color) #f"\\texttt{{{nt}}}"

    # Remove ticks and customize the plot
    ax.tick_params(bottom=False, left=False)
    ax.set_facecolor("white")
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])

    if save_plot:
        plt.savefig(save_plot)
    # Display the plot
    plt.show()

plt.rcParams['text.usetex'] = True
ss =        '..((((((..(((((((...<<<(((......))).(((((((..[[[[))))))).))))))).))))))........>>>..]]]]'
design_01 = 'CCGGCGGGUUGCAGUGGGCAGCACGCUAACACGCGACGGGAGUUUGGUCGCUCCCGACUACUGCCCCCGCCAAAAAAUUUGUGAGACC'
plot_colored_seq(design_01, ss, FV='TBFV')

ss =        '..((((((.((((((...<<<(((........))).(((((..[[[[[...))))).))))))..))))))...........>>>..]]]]]'
design_09 = 'CCGCCCAGGGCCUAGGCAGCACGCGAAAUUAGGCGACGGGAGAGGGUCGGAUCCCGACUGGGCAGCUGGGCAAUUUGUAAUUUGUGAGACCC'
plot_colored_seq(design_09, ss, FV='TBFV')


ss =        '<<.((((((((((....))))).(((((([[[[[...))))))>>))))).]]]]]..'
design_04 = 'AGUCAGGCCCAGUUAAUGCUGGCCACGGUAUUCCCCAACCGUGCUGCCUGUGGAAUUU'
#plot_colored_seq(design_04, ss, FV='MBFV')


ss =        '<<.(((((((((.......)))).((((((([[[[[[.)))))))>>)))))...]]]]]]..'
design_02 = 'UGUCAGGCCUAGGAAAAGACUAGCCACGGAUGCGUUCAAUUCGUGCAGCCUGUUUGAGUGUUU'
#plot_colored_seq(design_02, ss, FV='MBFV')

