import sympy as sym
from vhegen.modules.vibrational import Term, Coeff, Mono, Sin, Cos
from vhegen.modules.input import return_acceptable_states

requirements_dct = {0:[1,['A','A']],
                    1:[1,['A','B']],
                    2:[1,['B','B']],
                    3:[1,['E','A']],
                    4:[1,['E','B']],
                    5:[1,['E','E']],
                    6:[-1,['A','A']],
                    7:[-1,['A','B']],
                    8:[-1,['B','B']],
                    9:[-1,['E','A']],
                    10:[-1,['E','B']],
                    11:[-1,['E','E']],
                    12:[1j,['A','A']],
                    13:[1j,['A','B']],
                    14:[1j,['B','B']],
                    15:[1j,['E','A']],
                    16:[1j,['E','B']],
                    17:[1j,['E','E']],
                    18:[0.5*(-1+1j*sym.sqrt(3)),['A','A']],
                    19:[0.5*(-1+1j*sym.sqrt(3)),['E','A']],
                    20:[0.5*(-1+1j*sym.sqrt(3)),['E','E']]}


def return_formula(n,modes,symmetry):
    #n: rotation number
    #modes: list of first letter of modes eg. a,b,e
    #formula_dct = {0:[Term(1,Coeff('ar',['I1','bI2']),[Mono('z','_alpha','I1'),Mono('z','_beta','bI2')]),
    #                  Term(1j,Coeff('ai',['I3','bI4']),[Mono('z','_alpha','I3'),Mono('z','_beta','bI4')])],
    #
    #               1:[Term(1,Coeff('ar',['I1','2J']),[Mono('z','','I1'),Mono('w','','2J')]),
    #                  Term(1j,Coeff('ai',['I2','2J']),[Mono('z','','I2'),Mono('w','','2J')])],
    #
    #               2:[Term(1,Coeff('ar',['2J1+1','2J2+1']),[Mono('w','_alpha','2J1+1'),Mono('w','_beta','2J2+1')]),
    #                  Term(1j,Coeff('ai',['2J1+1','2J2+1']),[Mono('w','_alpha','2J1+1'),Mono('w','_beta','2J2+1')]),
    #                  Term(1,Coeff('ar',['2J1','2J2']),[Mono('w','_alpha','2J1'),Mono('w','_beta','2J2')]),
    #                  Term(1j,Coeff('ai',['2J1','2J2']),[Mono('w','_alpha','2J1'),Mono('w','_beta','2J2')])],
    #
    #                 3:[Term(1,Coeff('ar',['I1','2K',n+'m']),[Mono('z','','I1'),Mono('rho','','Abs('+n+'m)+2K'),Cos([n+'m'])]),
    #                  Term(-1,Coeff('ai',['I2','2K',n+'m']),[Mono('z','','I2'),Mono('rho','','Abs('+n+'m)+2K'),Sin([n+'m'])]),
    #                  Term(1j,Coeff('ar',['I1','2K',n+'m']),[Mono('z','','I1'),Mono('rho','','Abs('+n+'m)+2K'),Sin([n+'m'])]),
    #                  Term(1j,Coeff('ai',['I2','2K',n+'m']),[Mono('z','','I2'),Mono('rho','','Abs('+n+'m)+2K'),Cos([n+'m'])])],
    #
    #               4:[Term(1,Coeff('ar',['2I','2K','2m']),[Mono('w','','Mod(2m/2,2)+2I'),Mono('rho','','Abs(2m)+2K'),Cos(['2m'])]),
    #                  Term(-1,Coeff('ai',['2I','2K','2m']),[Mono('w','','Mod(2m/2,2)+2I'),Mono('rho','','Abs(2m)+2K'),Sin(['2m'])]),
    #                  Term(1j,Coeff('ar',['2I','2K','2m']),[Mono('w','','Mod(2m/2,2)+2I'),Mono('rho','','Abs(2m)+2K'),Sin(['2m'])]),
    #                  Term(1j,Coeff('ai',['2I','2K','2m']),[Mono('w','','Mod(2m/2,2)+2I'),Mono('rho','','Abs(2m)+2K'),Cos(['2m'])])],
    #
    #               5:[Term(1,Coeff('ar',['2K1','2K2','m1',n+'n']),[Mono('rho','_alpha','Abs(m1)+2K1'),Mono('rho','_beta','Abs('+n+'n-m1)+2K2'),Cos(['m1',n+'n-m1'])]),
    #                  Term(-1,Coeff('ai',['2K1','2K2','m1',n+'n']),[Mono('rho','_alpha','Abs(m1)+2K1'),Mono('rho','_beta','Abs('+n+'n-m1)+2K2'),Sin(['m1',n+'n-m1'])]),
    #                  Term(1j,Coeff('ar',['2K1','2K2','m1',n+'n']),[Mono('rho','_alpha','Abs(m1)+2K1'),Mono('rho','_beta','Abs('+n+'n-m1)+2K2'),Sin(['m1',n+'n-m1'])]),
    #                  Term(1j,Coeff('ai',['2K1','2K2','m1',n+'n']),[Mono('rho','_alpha','Abs(m1)+2K1'),Mono('rho','_beta','Abs('+n+'n-m1)+2K2'),Cos(['m1',n+'n-m1'])])],
    #          
    #               6:[0],
    #
    #               7:[Term(1,Coeff('br',['I1','2J+1']),[Mono('z','','I1'),Mono('w','','2J+1')]),
    #                  Term(1j,Coeff('bi',['I2','2J+1']),[Mono('z','','I2'),Mono('w','','2J+1')])],
    #
    #               8:[Term(1,Coeff('br',['2J1+1','2J2']),[Mono('w','_alpha','2J1+1'),Mono('w','_beta','2J2')]),
    #                  Term(1j,Coeff('bi',['2J1+1','2J2']),[Mono('w','_alpha','2J1+1'),Mono('w','_beta','2J2')]),
    #                  Term(1,Coeff('br',['2J1','2J2+1']),[Mono('w','_alpha','2J1'),Mono('w','_beta','2J2+1')]),
    #                  Term(1j,Coeff('bi',['2J1','2J2+1']),[Mono('w','_alpha','2J1'),Mono('w','_beta','2J2+1')])],
    #
    #               9:[Term(1,Coeff('br',['I1','2K',n+'n+2']),[Mono('z','','I1'),Mono('rho','','Abs('+n+'n+2)+2K'),Cos([n+'n+2'])]),
    #                  Term(-1,Coeff('bi',['I2','2K',n+'n+2']),[Mono('z','','I2'),Mono('rho','','Abs('+n+'n+2)+2K'),Sin([n+'n+2'])]),
    #                  Term(1j,Coeff('br',['I1','2K',n+'n+2']),[Mono('z','','I1'),Mono('rho','','Abs('+n+'n+2)+2K'),Sin([n+'n+2'])]),
    #                  Term(1j,Coeff('bi',['I2','2K',n+'n+2']),[Mono('z','','I2'),Mono('rho','','Abs('+n+'n+2)+2K'),Cos([n+'n+2'])])],
    #
    #               10:[Term(1,Coeff('br',['2I+1','2K','2m']),[Mono('w','','Mod(2m/2,2)+2I+1'),Mono('rho','','Abs(2m)+2K'),Cos(['2m'])]),
    #                   Term(-1,Coeff('bi',['2I+1','2K','2m']),[Mono('w','','Mod(2m/2,2)+2I+1'),Mono('rho','','Abs(2m)+2K'),Sin(['2m'])]),
    #                   Term(1j,Coeff('br',['2I+1','2K','2m']),[Mono('w','','Mod(2m/2,2)+2I+1'),Mono('rho','','Abs(2m)+2K'),Sin(['2m'])]),
    #                   Term(1j,Coeff('bi',['2I+1','2K','2m']),[Mono('w','','Mod(2m/2,2)+2I+1'),Mono('rho','','Abs(2m)+2K'),Cos(['2m'])])],
    #
    #               11:[Term(1,Coeff('br',['2K1','2K2','m1',n+'n+2']),[Mono('rho','_alpha','Abs(m1)+2K1'),Mono('rho','_beta','Abs('+n+'n+2-m1)+2K2'),Cos(['m1',n+'n+2-m1'])]),
    #                   Term(-1,Coeff('bi',['2K1','2K2','m1',n+'n+2']),[Mono('rho','_alpha','Abs(m1)+2K1'),Mono('rho','_beta','Abs('+n+'n+2-m1)+2K2'),Sin(['m1',n+'n+2-m1'])]),
    #                   Term(1j,Coeff('br',['2K1','2K2','m1',n+'n+2']),[Mono('rho','_alpha','Abs(m1)+2K1'),Mono('rho','_beta','Abs('+n+'n+2-m1)+2K2'),Sin(['m1',n+'n+2-m1'])]),
    #                   Term(1j,Coeff('bi',['2K1','2K2','m1',n+'n+2']),[Mono('rho','_alpha','Abs(m1)+2K1'),Mono('rho','_beta','Abs('+n+'n+2-m1)+2K2'),Cos(['m1',n+'n+2-m1'])])],
    #               
    #               12:[0], 
    #
    #               13:[0],
    #
    #               14:[0],
    #
    #               15:[Term(1,Coeff('cr',['I1','2K',n+'n-1']),[Mono('z','','I1'),Mono('rho','','Abs('+n+'n-1)+2K'),Cos([n+'n-1'])]),
    #                   Term(-1,Coeff('ci',['I2','2K',n+'n-1']),[Mono('z','','I2'),Mono('rho','','Abs('+n+'n-1)+2K'),Sin([n+'n-1'])]),
    #                   Term(1j,Coeff('cr',['I1','2K',n+'n-1']),[Mono('z','','I1'),Mono('rho','','Abs('+n+'n-1)+2K'),Sin([n+'n-1'])]),
    #                   Term(1j,Coeff('ci',['I2','2K',n+'n-1']),[Mono('z','','I2'),Mono('rho','','Abs('+n+'n-1)+2K'),Cos([n+'n-1'])])],
    #               
    #               16:[Term(1,Coeff('cr',['2I','2K','2n-1']),[Mono('w','','Mod((2n-1+1)/2,2)+2I'),Mono('rho','','Abs(2n-1)+2K'),Cos(['2n-1'])]),
    #                   Term(-1,Coeff('ci',['2I','2K','2n-1']),[Mono('w','','Mod((2n-1+1)/2,2)+2I'),Mono('rho','','Abs(2n-1)+2K'),Sin(['2n-1'])]),
    #                   Term(1j,Coeff('cr',['2I','2K','2n-1']),[Mono('w','','Mod((2n-1+1)/2,2)+2I'),Mono('rho','','Abs(2n-1)+2K'),Sin(['2n-1'])]),
    #                   Term(1j,Coeff('ci',['2I','2K','2n-1']),[Mono('w','','Mod((2n-1+1)/2,2)+2I'),Mono('rho','','Abs(2n-1)+2K'),Cos(['2n-1'])])],
    #               
    #               17:[Term(1,Coeff('cr',['2K1','2K2','m1',n+'n-1']),[Mono('rho','_alpha','Abs(m1)+2K1'),Mono('rho','_beta','Abs('+n+'n-1-m1)+2K2'),Cos(['m1',n+'n-1-m1'])]),
    #                   Term(-1,Coeff('ci',['2K1','2K2','m1',n+'n-1']),[Mono('rho','_alpha','Abs(m1)+2K1'),Mono('rho','_beta','Abs('+n+'n-1-m1)+2K2'),Sin(['m1',n+'n-1-m1'])]),
    #                   Term(1j,Coeff('cr',['2K1','2K2','m1',n+'n-1']),[Mono('rho','_alpha','Abs(m1)+2K1'),Mono('rho','_beta','Abs('+n+'n-1-m1)+2K2'),Sin(['m1',n+'n-1-m1'])]),
    #                   Term(1j,Coeff('ci',['2K1','2K2','m1',n+'n-1']),[Mono('rho','_alpha','Abs(m1)+2K1'),Mono('rho','_beta','Abs('+n+'n-1-m1)+2K2'),Cos(['m1',n+'n-1-m1'])])],
    #               
    #               18:[0],
    #               
    #               19:[Term(1,Coeff('br',['I1','2K',n+'n-1']),[Mono('z','','I1'),Mono('rho','','Abs('+n+'n-1)+2K'),Cos([n+'n-1'])]),
    #                   Term(-1,Coeff('bi',['I2','2K',n+'n-1']),[Mono('z','','I2'),Mono('rho','','Abs('+n+'n-1)+2K'),Sin([n+'n-1'])]),
    #                   Term(1j,Coeff('br',['I1','2K',n+'n-1']),[Mono('z','','I1'),Mono('rho','','Abs('+n+'n-1)+2K'),Sin([n+'n-1'])]),
    #                   Term(1j,Coeff('bi',['I2','2K',n+'n-1']),[Mono('z','','I2'),Mono('rho','','Abs('+n+'n-1)+2K'),Cos([n+'n-1'])])],
    #               
    #               20:[Term(1,Coeff('br',['2K1','2K2','m1',n+'n-1']),[Mono('rho','_alpha','Abs(m1)+2K1'),Mono('rho','_beta','Abs('+n+'n-1-m1)+2K2'),Cos(['m1',n+'n-1-m1'])]),
    #                   Term(-1,Coeff('bi',['2K1','2K2','m1',n+'n-1']),[Mono('rho','_alpha','Abs(m1)+2K1'),Mono('rho','_beta','Abs('+n+'n-1-m1)+2K2'),Sin(['m1',n+'n-1-m1'])]),
    #                   Term(1j,Coeff('br',['2K1','2K2','m1',n+'n-1']),[Mono('rho','_alpha','Abs(m1)+2K1'),Mono('rho','_beta','Abs('+n+'n-1-m1)+2K2'),Sin(['m1',n+'n-1-m1'])]),
    #                   Term(1j,Coeff('bi',['2K1','2K2','m1',n+'n-1']),[Mono('rho','_alpha','Abs(m1)+2K1'),Mono('rho','_beta','Abs('+n+'n-1-m1)+2K2'),Cos(['m1',n+'n-1-m1'])])]}
    amodes=return_acceptable_states(symmetry)
    modelabels=[]
    for i in range(len(modes)):
       for j in range(len(amodes.full)):
          if modes[i]==amodes.full[j]:
             modelabels.append(amodes.labels[j])
    #print(modes)
    #print(modelabels)
    #exit()
    #label and count each modes
    posa=[]
    numa=0
    posb=[]
    numb=0
    pose=[]
    nume=0
    modes0=[mode[0] for mode in modes]
    for i in range(len(modes)):
       if modes0[i]=='E':
          pose.append(i)
          nume=nume+1
       if modes0[i]=='A':
          posa.append(i)
          numa=numa+1
       if modes0[i]=='B':
          posb.append(i)
          numb=numb+1

    indices=[]#list of indices
    
    if nume>0:
       #full expansion with cos and sin required
       funcs_rcos=[]
       funcs_rsin=[]
       funcs_isin=[]
       funcs_icos=[]
       for i in range(numa):
          func=Mono('z',f'_{i+1}',f'I_{i+1}')
          indices.append(f'I_{i+1}')
          funcs_rcos.append(func)
          funcs_rsin.append(func)
          funcs_isin.append(func)
          funcs_icos.append(func)
       for i in range(numb):
          func=Mono('w',f'_{i+1}',f'J_{i+1}')
          indices.append(f'J_{i+1}')
          funcs_rcos.append(func)
          funcs_rsin.append(func)
          funcs_isin.append(func)
          funcs_icos.append(func)
       for i in range(nume):
          func=Mono('rho',f'_{i+1}',f'Abs(M_{i+1})+2*K_{i+1}')
          indices.append(f'M_{i+1}')
          indices.append(f'K_{i+1}')
          subind=createind(modelabels[pose[i]],i+1)
          #print(subind)
          #exit()
          #func=Mono('\\rho','_'+subind,f'Abs(M_{i+1})+2*K_{i+1}')
          #indices.append(f'M_{i+1}')
          #indices.append(f'K_{i+1}')
          funcs_rcos.append(func)
          funcs_rsin.append(func)
          funcs_isin.append(func)
          funcs_icos.append(func)
       #create list of phi for cos and sin terms
       philist=[]
       for i in range(nume):
          philist.append(f'M_{i+1}')
          
       cfunc=Cos(philist)
       sfunc=Sin(philist)
       #cfunc.compute()
       #cfunc.compile()

       #print(cfunc.sympy_form)
       #exit()
       #cfunc.symcomp()
       #print(cfunc.symb)
       funcs_rcos.append(cfunc)
       funcs_rsin.append(sfunc)
       funcs_isin.append(sfunc)
       funcs_icos.append(cfunc)

       return [Term(1,Coeff('cr',indices),funcs_rcos),
                       Term(-1,Coeff('ci',indices),funcs_rsin),
                       Term(1j,Coeff('cr',indices),funcs_isin),
                       Term(1j,Coeff('ci',indices),funcs_icos)]
    else:
       #no cos and sin required
       funcs_r=[]
       funcs_i=[]
       for i in range(numa):
          func=Mono('z',f'_{i+1}',f'I_{i+1}')
          indices.append(f'I_{i+1}')
          funcs_r.append(func)
          funcs_i.append(func)
       for i in range(numb):
          func=Mono('w',f'_{i+1}',f'J_{i+1}')
          indices.append(f'J_{i+1}')
          funcs_r.append(func)
          funcs_i.append(func)


       return [Term(1,Coeff('cr',indices),funcs_r),
                       Term(1,Coeff('ci',indices),funcs_i)]
def createind(modelabel,i):
   if modelabel[1]==0 and modelabel[2]==0:
      return f'_{i}'

   if modelabel[1]>0 and modelabel[2]==0:
      return f'{{{modelabel[1]}_{i}}}'
#EOF
