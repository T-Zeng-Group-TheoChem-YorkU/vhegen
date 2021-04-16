trigonal_groups = ('CN','CNv','CNh','DN','DNh','DNd')

tetragonal_groups = ('SN','Cinfv','Dinfv')

cubic_groups = ()

all_groups = trigonal_groups+tetragonal_groups+cubic_groups

irrep_dct = {'C3':("E","A"),
             'C3V':("A1","A2","E"),
             'C3H':("A'","E'","A''","E''"),
             'D3':("A1","A2","E"),
             'D3H':("A1'","A2'","E'","A1''","A2''","E''"),
             'D3D':("A1G","A2G","EG","A1U","A2U","EU"),
             'C4':("A","B","E"),
             'S4':("A","B","E"),
             'C4V':("A1","A2","B1","B2","E"),
             'C4H':("AG","BG","EG","AU","BU","EU"),
             'D4':("A1","A2","B1","B2","E"),
             'D4H':("A1G","A2G","B1G","B2G","EG","A1U","A2U","B1U","B2U","EU"),
             'D2D':("A1","A2","B1","B2","E")}

irrep_priority = ['E','A','B']

state_components_dct = {'A':['A'],
                        'B':['B'],
                        'E':['+','-'],
                        'T':['X','Y','Z']}

escapes = ''.join([chr(i) for i in range(1,32)])

def replace_all(string, dct):
    for k in dct:
        string = string.replace(k,dct[k])
    return string
#EOF
