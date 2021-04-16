from argparse import ArgumentParser
import vhegen.modules.glbls as glo
import os.path

try:
    input = raw_input
except NameError:
    pass

def configure_parser():
    parser = ArgumentParser()
    parser.add_argument("--c",dest="c",
                        help="config file path, default = config.cfg",
                        metavar="CONFIG", default="config.cfg")
    parser.add_argument("--sym",dest="sym",
                        help="point group static input",
                        metavar="SYMMETRY")
    parser.add_argument("--states",dest="states",
                        help="electronic states static input",
                        metavar="STATES")
    parser.add_argument("--modes",dest="modes",
                        help="vibrational modes static input.",
                        metavar="MODES")
    parser.add_argument("--o",dest="o",
                        help="order(s) of expansion static input.",
                        metavar="ORDERS")
    parser.add_argument("--f",dest="f",
                        help="output filename. default = 'output'",
                        metavar="FILENAME", default="output")
    return parser

def read_config(path):
    if (os.path.exists(path)):
        configfile = open(path,"r").readlines()
        config = {}
        for line in configfile:
            line = line.replace(u'\n','')
            if line[0] == '#':
                pass
            else:
                arg = line.split('=')
                config[arg[0]] = arg[1]
        return config
    else:
        return {'input':'dynamic',
            'pdf_out':'true',
            'log_out':'true',
            'e_coords':'both',
            'basis':'both',
            'mctdh_out':'true'}

def symmetry_arg(arg):
    
    arg= arg.replace(" ","")
    arg = arg.upper()

    if arg == 'LIST':
        print('\nAvailable point groups:')
        for s in glo.all_groups:
            print(s)
        raise Exception()
    elif arg == 'EXIT':
        exit()
    #elif symmetry in glo.all_groups:
    #    return symmetry
    else:
        l=len(arg)
        if arg[0].isalpha():
            valid=('C','D','S')
            lettornum=[]
            lettornum.append('l')
            if arg[0] in valid:
                letter=arg[0]
            else:
                raise Exception('InputError: First letter must be in:'+'C'+'D'+'S')
        else:
                raise Exception('InputError: First input must be a letter')
        if (arg[1:4]=='INF'): 
            rot=0
            if (arg[0]=='C' and arg[4]=='V'):
                extra='V'
                refl=True
                return {'letter':letter,
                    'rot':rot,
                    'extra':extra,
                    'refl':refl,
                    'print':letter+'inf'+extra.lower()}
            elif (arg[0]=='D' and arg[4]=='H'):
                extra='H'
                refl=True
                return {'letter':letter,
                    'rot':rot,
                    'extra':extra,
                    'refl':refl,
                    'print':letter+'inf'+extra.lower()}
            else:
                raise Exception('Not a valid symmetry with inf rotation')
        for i in range(1,l):
            if arg[i].isdigit():
                lettornum.append('n')
            elif arg[i].isalpha():
                lettornum.append('l')
            else:
                raise Exception('InputError: Bad character')

        if lettornum[1]!='n':
            raise Exception('InputError: second character must be a number')
        if l==2:
            rot=int(arg[1])
            extra=''
        elif l==3:
            if lettornum[2]=='n':
                rot=int(arg[1:3])
                extra=''
            elif lettornum[2]=='l':
                rot=int(arg[1])
                extra=arg[2]
        elif l==4:
            if lettornum[2]=='n':
                rot=int(arg[1:3])
                if (lettornum[3]=='l'):
                    extra=arg[3]
                else:
                    rot=int(arg[1:4])
                    extra=''
            else:
                raise Exception('InputError: not a symmetry group')
        elif l==5:
            if lettornum[2]=='n' and lettornum[3]=='n':
                rot=int(arg[1:4])
                if (lettornum[4]=='l'):
                    extra=arg[4]
                else:
                    raise Exception('InputError: you want rotation symmetry over 999?')
            else:
                raise Exception('InputError: not a symmetry group') 
        else:
            raise Exception('InputError: you want rotation symmetry over 999?')       
        
        valid=('','H','V','D')
        if extra in valid:
            if extra=='D' and letter!='D':
                raise Exception('InputError: not a valid symmetry group')
            if extra=='V' and letter!='C':
                raise Exception('InputError: not a valid symmetry group')
            if letter=='S' and extra!='':
                raise Exception('InputError: not a valid symmetry group') 
            if letter=='S' and rot%2!=0:
                raise Exception('InputError: not a valid symmetry group') 

        else:
            raise Exception('InputError: Last character must be in'+' h,v,d, ') 
        
        if extra =='D' and rot%2==0:
            rot=rot*2
        
        refl=False
        return {'letter':letter,
            'rot':rot,
            'extra':extra,
            'refl':refl,
            'print':letter+str(rot)+extra.lower()}

class state:
    def __init__(self):
        self.full=[]
        self.labels=[]
        self.size=0

    def getsize(self):
        self.size=len(self.print)

def return_acceptable_states(sym):
    states=state()
    
    if (sym['letter']=='C'):
        if (sym['rot']==0): #Cinfv
            states.full.append('A1')
            states.labels.append(('A',1,0,''))
            states.full.append('A2')
            states.labels.append(('A',2,0,''))
            maxemode=999
            for i in range(maxemode):
                states.full.append('E'+str(i+1))
                states.labels.append(('E',i+1,0,''))
            sym['refl']=True
        elif (sym['rot']%2==0):
            if (sym['extra']==''):
                states.full.append('A')
                states.labels.append(('A',0,0,''))
                states.full.append('B')
                states.labels.append(('B',0,0,''))
                emax=sym['rot']//2-1
                for i in range(emax):
                    states.full.append('E'+str(i+1))
                    states.labels.append(('B',i+1,0,''))
            elif (sym['extra']=='V'):
                states.full.append('A1')
                states.labels.append(('A',1,0,''))
                states.full.append('B1')
                states.labels.append(('B',1,0,''))
                states.full.append('A2')
                states.labels.append(('A',2,0,''))
                states.full.append('B2')
                states.labels.append(('B',2,0,''))
                emax=sym['rot']//2-1
                for i in range(emax):
                    states.full.append('E'+str(i+1))
                    states.labels.append(('E',i+1,0,''))
                sym['refl']=True
                
            elif (sym['extra']=='H'):
                states.full.append('Ag')
                states.labels.append(('A',0,1,'g'))
                states.full.append('Bg')
                states.labels.append(('B',0,1,'g'))
                emax=sym['rot']//2-1
                for i in range(emax):
                    states.full.append('E'+str(i+1)+'g')
                    states.labels.append(('E',i+1,1,'g'))
                states.full.append('Au')
                states.labels.append(('A',0,2,'u'))
                states.full.append('Bu')
                states.labels.append(('B',0,2,'u'))
                for i in range(emax):
                    states.full.append('E'+str(i+1)+'u')
                    states.labels.append(('E',i+1,2,'u'))
        elif (sym['rot']%2==1):
            
            if (sym['extra']==''):
                #print(sym['print'])
                states.full.append('A')
                states.labels.append(('A',0,0,''))
                emax=(sym['rot']+1)//2-1
                for i in range(emax):
                    states.full.append('E'+str(i+1))
                    states.labels.append(('E',i+1,0,''))
            elif (sym['extra']=='V'):
                states.full.append('A1')
                states.labels.append(('A',1,0,''))
                states.full.append('A2')
                states.labels.append(('A',2,0,''))
                emax=(sym['rot']+1)//2-1
                for i in range(emax):
                    states.full.append('E'+str(i+1))
                    states.labels.append(('E',i+1,0,''))
                sym['refl']=True
            elif (sym['extra']=='H'):
                states.full.append("A'")
                states.labels.append(('A',0,1,"'"))
                emax=(sym['rot']+1)//2-1
                for i in range(emax):
                    states.full.append('E'+str(i+1)+"'")
                    states.labels.append(('E',i+1,1,"'"))
                states.full.append("A''")
                states.labels.append(('A',0,2,"''"))
                for i in range(emax):
                    states.full.append('E'+str(i+1)+"''")
                    states.labels.append(('E',i+1,2,"''"))
    elif (sym['letter']=='D'):
        if (sym['rot']==0):
            if (sym['extra']=='H'):
                states.full.append('A1g')
                states.labels.append(('A',1,1,'g'))
                states.full.append('A2g')
                states.labels.append(('A',2,1,'g'))
                emax=999
                for i in range(emax):
                    states.full.append('E'+str(i+1)+'g')
                    states.labels.append(('E',i+1,1,'g'))
                states.full.append('A1u')
                states.labels.append(('A',1,2,'u'))
                states.full.append('A2u')
                states.labels.append(('A',2,2,'u'))
                emax=999
                for i in range(emax):
                    states.full.append('E'+str(i+1)+'u')
                    states.labels.append(('E',i+1,2,'u'))
                sym['refl']=True

        elif (sym['rot']%2==0):
            if (sym['extra']=='H'):
                states.full.append('A1g')
                states.labels.append(('A',1,1,'g'))
                states.full.append('A2g')
                states.labels.append(('A',2,1,'g'))
                states.full.append('B1g')
                states.labels.append(('B',1,1,'g'))
                states.full.append('B2g')
                states.labels.append(('B',2,1,'g'))
                emax=sym['rot']//2-1
                for i in range(emax):
                    states.full.append('E'+str(i+1)+'g')
                    states.labels.append(('E',i+1,1,'g'))
                states.full.append('A1u')
                states.labels.append(('A',1,2,'u'))
                states.full.append('A2u')
                states.labels.append(('A',2,2,'u'))
                states.full.append('B1u')
                states.labels.append(('B',1,2,'u'))
                states.full.append('B2u')
                states.labels.append(('B',2,2,'u'))
                emax=sym['rot']//2-1
                for i in range(emax):
                    states.full.append('E'+str(i+1)+'u')
                    states.labels.append(('E',i+1,2,'u'))
                sym['refl']=True
                
            elif (sym['extra']=='D'):
                states.full.append('A1')
                states.labels.append(('A',1,0,''))
                states.full.append('A2')
                states.labels.append(('A',2,0,''))
                states.full.append('B1')
                states.labels.append(('B',1,0,''))
                states.full.append('B2')
                states.labels.append(('B',2,0,''))
                emax=sym['rot']-1
                for i in range(emax):
                    states.full.append('E'+str(i+1))
                    states.labels.append(('E',i+1,0,''))
                sym['refl']=True

        elif (sym['rot']%2==1):
            if (sym['extra']=='H'):
                states.full.append("A1'")
                states.labels.append(('A',1,1,"'"))
                states.full.append("A2'")
                states.labels.append(('A',2,1,"'"))
                emax=sym['rot']//2
                for i in range(emax):
                    states.full.append('E'+str(i+1)+"'")
                    states.labels.append(('E',i+1,1,"'"))
                states.full.append("A1''")
                states.labels.append(('A',1,2,"''"))
                states.full.append("A2''")
                states.labels.append(('A',2,2,"''"))
                emax=sym['rot']//2
                for i in range(emax):
                    states.full.append('E'+str(i+1)+"''")
                    states.labels.append(('E',i+1,1,"''"))
                sym['refl']=True
            elif (sym['extra']=='D'):
                states.full.append("A1g")
                states.labels.append(('A',1,1,'g'))
                states.full.append("A2g")
                states.labels.append(('A',2,1,'g'))
                emax=sym['rot']//2
                for i in range(emax):
                    states.full.append('E'+str(i+1)+'g')
                    states.labels.append(('E',i+1,1,'g'))
                states.full.append("A1u")
                states.labels.append(('A',1,2,'u'))
                states.full.append("A2u")
                states.labels.append(('A',2,2,'u'))
                emax=sym['rot']//2
                for i in range(emax):
                    states.full.append('E'+str(i+1)+'u')
                    states.labels.append(('E',i+1,2,'u'))
                sym['refl']=True
    elif (sym['letter']=='S'):
        if (sym['rot']%4==0):
            if (sym['extra']==''):
                states.full.append('A')
                states.labels.append(('A',0,0,''))
                states.full.append('B')
                states.labels.append(('B',0,0,''))
                emax=sym['rot']//2-1
                for i in range(emax):
                    states.full.append('E'+str(i+1))
                    states.labels.append(('E',i+1,0,''))
        elif (sym['rot']%4==2):
            if (sym['extra']==''):
                states.full.append('Ag')
                states.labels.append(('A',0,1,'g'))
                emax=sym['rot']//4
                for i in range(emax):
                    states.full.append('E'+str(i+1)+'g')
                    states.labels.append(('E',i+1,1,'g'))
                states.full.append('Au')
                states.labels.append(('A',0,2,'u'))
                for i in range(emax):
                    states.full.append('E'+str(i+1)+'u')
                    states.labels.append(('E',i+1,2,'u'))              
    return states

def printfancy(acceptable):
    l=len(acceptable)
    #print(l)
    j=0
    for i in range((l-1)//5):
        #print(acceptable[j:j+5])
        print(f'{j*5:3d}:{acceptable[j*5]:6} {j*5+1:3d}:{acceptable[j*5+1]:6} {j*5+2:3d}:{acceptable[j*5+2]:6} {j*5+3:3d}:{acceptable[j*5+3]:6} {j*5+4:3d}:{acceptable[j*5+4]:6}')
        j=j+1
    if (l-j*5==1):
        #print(acceptable[j*5::])
        print(f'{j*5:3d}:{acceptable[j*5]:6}')
    elif (l-j*5==2):
        print(f'{j*5:3d}:{acceptable[j*5]:6} {j*5+1:3d}:{acceptable[j*5+1]:6}')
    elif (l-j*5==3):
        print(f'{j*5:3d}:{acceptable[j*5]:6} {j*5+1:3d}:{acceptable[j*5+1]:6} {j*5+2:3d}:{acceptable[j*5+2]:6}')
    elif (l-j*5==4):
        print(f'{j*5:3d}:{acceptable[j*5]:6} {j*5+1:3d}:{acceptable[j*5+1]:6} {j*5+2:3d}:{acceptable[j*5+2]:6} {j*5+3:3d}:{acceptable[j*5+3]:6}')
    elif (l-j*5==5):
        print(f'{j*5:3d}:{acceptable[j*5]:6} {j*5+1:3d}:{acceptable[j*5+1]:6} {j*5+2:3d}:{acceptable[j*5+2]:6} {j*5+3:3d}:{acceptable[j*5+3]:6} {j*5+4:3d}:{acceptable[j*5+4]:6}')

def addones(arg,sym):
    #appends 1 to e states if input is 
    #Cnx, n<=4
    #Dnd, n<=3
    #Dnh, n<=4
    #Dn,  n<=4
    #Sn,  n<=4
    numarg=len(arg)
    if (sym['letter']=='C' and sym['rot']<=4):
        for i in range(numarg):
            if arg[i][0].upper()=='E':
                if (len(arg[i])>1):
                    if arg[i][1]!=1:
                        arg[i]=arg[i][0]+'1'+arg[i][1:]
                else:
                    arg[i]=arg[i]+'1'
    if (sym['letter']=='D' and sym['extra']=='' and sym['rot']<=4):
        for i in range(numarg):
            if arg[i][0].upper()=='E':
                if (len(arg[i])>1):
                    if arg[i][1]!=1:
                        arg[i]=arg[i][0]+'1'+arg[i][1:]
                else:
                    arg[i]=arg[i]+'1'
    if (sym['letter']=='D' and sym['extra']=='D' and sym['rot']<=3):
        for i in range(numarg):
            if arg[i][0].upper()=='E':
                if (len(arg[i])>1):
                    if arg[i][1]!=1:
                        arg[i]=arg[i][0]+'1'+arg[i][1:]
                else:
                    arg[i]=arg[i]+'1'
    if (sym['letter']=='H' and sym['extra']=='D' and sym['rot']<=4):
        for i in range(numarg):
            if arg[i][0].upper()=='E':
                if (len(arg[i])>1):
                    if arg[i][1]!=1:
                        arg[i]=arg[i][0]+'1'+arg[i][1:]
                else:
                    arg[i]=arg[i]+'1'
    if (sym['letter']=='S' and sym['rot']<=4):
        for i in range(numarg):
            if arg[i][0].upper()=='E':
                if (len(arg[i])>1):
                    if arg[i][1]!=1:
                        arg[i]=arg[i][0]+'1'+arg[i][1:]
                else:
                    arg[i]=arg[i]+'1'

    

def states_arg(arg,sym):
    arg = arg.upper()
    arg = arg.replace('"',"''")
    arg =  arg.replace(" ","")
    
    acceptable=return_acceptable_states(sym).full
    if arg == 'LIST':
        print('\nIrreps of '+sym['print']+':')
        #acceptable=return_acceptable_states(sym).full
        printfancy(acceptable)
        raise Exception()
    elif arg == 'EXIT':
        exit()
    if ',' in arg:
        states = arg.replace(" ","").split(',')
    elif '+' in arg:
        states = arg.replace(" ","").split('+')
    else:
        states = [arg.replace(" ","")]

    #print(states)
    addones(states,sym) 
    #print(states)
    #exit()
    max_cond = len(states)
    met_cond = 0
    for i in range(max_cond):
        if states[i].isdigit():
            if int(states[i])<len(acceptable):
                states[i]=acceptable[int(states[i])]
            else:
                raise Exception(f'InputError: Number for state {i+1} too large')
        
           
    #print(states)
    for i in range(len(states)):
        
        for a in acceptable:
            if (states[i].upper()==a.upper()):
                states[i]=a
                met_cond += 1
    if met_cond==max_cond:
        #print(states)
        return states
    if met_cond != max_cond:
        raise Exception('InputError: State(s) '+arg+' not valid irreps in '+str(sym["print"])+'.')

def modes_arg(arg,sym):
    arg = arg.upper()
    arg = arg.replace('"',"''")
    if arg == 'LIST':
        print('\nIrreps of '+sym+':')
        for e in glo.irrep_dct[sym]:
            print(e)
        raise Exception()
    elif arg == 'EXIT':
        exit()
    if ',' in arg:
        modes = arg.split(',')
    elif '+' in arg:
        modes = arg.split('+')
    else:
        modes = [arg]
    if len(modes) not in (1,2):
        raise Exception('InputError: Number of modes must be 1 or 2.')
    max_cond = len(modes)
    met_cond = 0
    for e in modes:
        if e in glo.irrep_dct[sym]:
            met_cond += 1
    if met_cond == max_cond:
        return reorder(modes)
    else:
        raise Exception('InputError: Modes(s) '+modes.lower()+' not valid irreps in '+symmetry+'.')

def orders_arg(arg):
    if arg.upper() == 'EXIT':
        exit()
    orders = arg.split(',')
    for i,order in enumerate(orders):
        try:
            orders[i] = int(order)
        except ValueError:
            raise Exception('InputError: Orders must be non-negative integers.')
    if len(orders) > 2:
        raise Exception('InputError: Invalid range of orders.')
    elif len(orders) == 2:
        if orders[1] <= orders[0]:
            raise Exception('InputError: Invalid range of orders.')
        if (orders[0]) < 0 or (orders[1] < 0):
            raise Exception('DynamicInputError: Orders must be non-negative integers.')
    elif len(orders) == 1:
        if orders[0] < 0:
            raise Exception('InputError: Orders must be non-negative integers.')
    if len(orders) == 2:
        orders = range(int(orders[0]),int(orders[1])+1)
    return orders
    
def filename_arg(arg):
    for char in (r"!@#$%^&*(){}[],."):
        if char in arg:
            print(char)
            raise Exception('InputError: Filename may only contain letters, numbers, "_", and "-".')
    if len(arg) == 0:
        return 'output'
    else:
        return arg

def format(inp_list):
    formatted_str = '('
    for i,n in enumerate(inp_list):
        if i>0 and i<len(inp_list):
            formatted_str+='+'
        formatted_str+=n
    formatted_str+=')'
    return formatted_str

def format_problem(states,modes):
    states_str = format(states)
    modes_str = format(modes)
    return states_str+'x'+modes_str.lower()

def reorder(inp_list):
    #re-orders states/modes by priority eg) (A+E)x(e+a) -> (E+A)x(e+a)
    if len(inp_list) > 1:
        if glo.irrep_priority.index(inp_list[0][0]) > glo.irrep_priority.index(inp_list[1][0]):
            inp_list = list(reversed(inp_list))
    return inp_list
    
def return_problem(symmetry,states,modes,orders,filename):
    return  ('\n---------------------------------------'
             '\n  VHEGEN instance parameters:'
             '\n  Symmetry: ' +symmetry['letter']+str(symmetry['rot'])+symmetry['extra']+
             '\n  Problem: ' + format_problem(states,modes)+
             '\n  Expansion orders: ' + str(orders)+
             '\n  Output filename: ' + filename+
             '\n---------------------------------------')

def dynamic_input():
    #symmetry
    dyn_inp = input('Enter symmetry: ')
    while True:
        try:
            sym_inp = symmetry_arg(dyn_inp)
        except Exception as e:
            print(e)
            dyn_inp = input('Re-enter symmetry: ')
            continue
        else:
            symmetry = sym_inp
            print(symmetry)
            print(symmetry['print']+' symmetry accepted.')
            break
    #states
    dyn_inp = input('Enter electronic state(s): ')
    while True:
        try:
            states_inp = states_arg(dyn_inp,symmetry)
            if len(states_inp) not in (1,2):
                raise Exception('InputError: Number of states must be 1 or 2.')
            #print(states_inp)
            states_inp = reorder(states_inp)
        except Exception as e:
            print(e)
            dyn_inp = input('Retry state(s): ')
            continue
        else:
            states = states_inp
            print('Electronic states '+ format(states)+' accepted.')
            break
    #modes
    dyn_inp = input('Enter vibrational mode(s): ')
    while True:
        try:
            modes_inp = states_arg(dyn_inp,symmetry)
            #modes_inp = reorder(modes_inp)
            #if symmetry in glo.trigonal_groups:
            #    modes_inp = modes_inp #limit_trigonal(modes_inp)
        except Exception as e:
            print(e)
            dyn_inp = input('Retry mode(s): ')
            continue
        else:
            modes = modes_inp
            print('Vibrational modes '+format(modes).lower()+' accepted.')
            break
    #orders
    dyn_inp = input('Order(s) of expansion:')
    while True:
        try:
            orders_inp = orders_arg(dyn_inp)
        except Exception as e:
            print(e)
            dyn_inp = input('Retry orders: ')
            continue
        else:
            orders = orders_inp
            print('Orders of expansion '+(str(list(orders)))+' accepted.')
            break
    #output filename
    dyn_inp = input('Enter filename: ')
    while True:
        try:
            filename_inp = filename_arg(dyn_inp)
        except Exception as e:
            print(e)
            dyn_inp = input('Retry filename: ')
            continue
        else:
            filename = filename_inp
            print('Filename: "'+filename+'".')
            break

    print(return_problem(symmetry,states,modes,orders,filename))

    while True:
        yn_prompt = input('Continue? [Y/N]: ')
        if yn_prompt.upper() == 'Y':
            print('')
            return symmetry, states, modes, orders, filename
        elif yn_prompt.upper() == 'N' or yn_prompt.upper() == 'EXIT':
            exit()

def read_input():
    argparser = configure_parser()
    options = argparser.parse_args()
    config = read_config('config.cfg')
    if 'static' in str(config[u'input']):
        print("Reading static inputs.")
        #static input
        try:
            symmetry = symmetry_arg(options.sym)
        except Exception as e:
            print('StaticInputError: Invalid symmetry argument.')
            exit()
        try:
            if len(options.states) not in (1,2):
                raise Exception('InputError: Number of states must be 1 or 2.')
            #states = reorder(states_arg(options.states,symmetry))
            #print(options.states)
            states = reorder(states_arg(options.states,symmetry))
            #states=reorder(states)
            #if len(states_inp) not in (1,2):
            #    raise Exception('InputError: Number of states must be 1 or 2.')
            #print(states_inp)
            #states_inp = reorder(states_inp)
            #if symmetry in glo.trigonal_groups:
            #    states = states #limit_trigonal(states)
        except Exception as e:
            print('StaticInputError: Invalid states argument.')
            exit()
        try:
            modes = states_arg(options.modes,symmetry)
            #modes = reorder(modes_arg(options.modes,symmetry))
            #if symmetry in glo.trigonal_groups:
            #    modes = modes #limit_trigonal(modes)
        except Exception as e:
            print(e)
            print('StaticInputError: Invalid modes argument.')
            exit()
        try:
            orders = orders_arg(options.o)
        except Exception as e:
            print('StaticInputError: Invalid orders argument.')
            exit()
        try:
            filename = filename_arg(options.f)
        except Exception as e:
            print('StaticInputError: Invalid filename argument.')
            exit()
    elif 'dynamic' in str(config[u'input']):
        print("Entering dynamic input.")
        #enter dynamic input
        symmetry,states,modes,orders,filename = dynamic_input()
    else:
        print("Error: 'input' value in config.cfg not found or recognized.")
        exit()
    problem_dct = {'sym':symmetry,
                   'states':states,
                   'modes':modes,
                   'o': orders,
                   'f': filename}
    return config,problem_dct

def prepare_input(sym,states,modes,orders,filename='output'):
    sym = symmetry_arg(sym)
    states = reorder(states_arg(states,sym))
    modes = states_arg(modes,sym)
    orders = orders_arg(orders)
    filename = filename_arg(filename)
    return {'sym':sym,
            'states':states,
            'modes':modes,
            'o': orders,
            'f':filename}
#EOF
