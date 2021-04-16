import sympy as sym
import itertools
import os
from vhegen.modules.glbls import replace_all, escapes
from vhegen.modules.input import reorder
from vhegen.modules.input import return_acceptable_states
from copy import copy,deepcopy
import numpy as np

class Term: #ExpansionTerm
    def __init__(self,prefactor,coeff,funcs):
        self.prefactor = prefactor      #prefactor should be any sympy-digestible argument
        self.coeff = coeff              #coeff should be an instance of Coeff
        self.indices = coeff.indices
        self.funcs = funcs              #funcs should be a list of instances of Mono and/or Trig 
        self.form = [self.prefactor,self.coeff] + self.funcs
        l=len(funcs)
        self.sinorcos='neither'
        if isinstance(funcs[l-1],Cos):
            self.sinorcos='cos'
        if isinstance(funcs[l-1],Sin):
            self.sinorcos='sin'
    
    def change_to_exp(self):
        if self.sinorcos=='sin':
            self.prefactor=sym.I
            self.update_form()
            self.funcs[len(self.funcs)-1].change_to_exp()
        if (self.sinorcos=='cos'):
            self.funcs[len(self.funcs)-1].change_to_exp()
    #1  
    def update_indices(self,p):
        #need to update them in coeff and funcs
        for c,i in enumerate(self.indices): #count, index
            #print(self.indices[c])
            if (i[0]=='K'):
                self.coeff.indices[c] = self.coeff.indices[c].replace(i,str(p[c]//2))
            else:
                self.coeff.indices[c] = self.coeff.indices[c].replace(i,str(p[c]))
            #print(self.indices[c])
            for c2,f in enumerate(self.funcs):
                if isinstance(f,Mono):
                    if (i[0]=='K'):
                        f.arg = f.arg.replace(i,str(p[c]//2))
                    else:
                        f.arg = f.arg.replace(i,str(p[c]))
                elif isinstance(f,Trig):
                    for c3,a in enumerate(f.args):
                        f.args[c3] = a.replace(i,str(p[c]))

    def refresh_indices(self):
        self.indices = self.coeff.indices

    def update_form(self):
        self.form = [self.prefactor,self.coeff] + self.funcs

    def get_order(self):
        order = 0
        for f in self.funcs:
            if isinstance(f,Mono):
                order += sym.sympify(f.arg)
        self.order = order

    def zero(self):
        self = 0

    def grab_constraints(self,constraint_dct): #checks global constraints, acquires relevant("local") constraints, returns them in index-position basis.
        applicable_constraints = {}
        for c in constraint_dct:
            if 'nz' in constraint_dct[c]:
                applicable_constraints[c] = constraint_dct[c]
            if c == 'all':
                applicable_constraints[c] = constraint_dct[c]
            else: #summing-index based constraints
                if '&' in c:
                    c_parse = c.split('&')
                    c_index = []
                    for i in self.indices:
                        for t in c_parse:
                            if t in i:
                                c_index.append(self.indices.index(i))
                    if len(c_index) == 2:
                        applicable_constraints['&'.join([str(i) for i in c_index])] = constraint_dct[c]
                else:
                    for i in self.indices:
                        if c in i:
                            if self.indices.index(i) in applicable_constraints:
                                applicable_constraints[self.indices.index(i)].add(constraint_dct[c])
                            else:
                                applicable_constraints[self.indices.index(i)] = set([constraint_dct[c]])
        self.local_constraints = applicable_constraints
        return self.local_constraints

    #def constrained_product(val_list):


    #2  
    def expand(self,order,eigenvals,symmetry,modes):
        from modules.input import return_acceptable_states
        #order: order of expansion
        #eigenvals: list of eigenvalue numbers
        
        acceptable=return_acceptable_states(symmetry)
        modelabels=[]
        for k in range(len(modes)):
            for m in range(len(acceptable.full)):
                if acceptable.full[m]==modes[k]:
                    modelabels.append(acceptable.labels[m])
        
        indpos, indlabels = assign_indices(modes,modelabels)

        #print(self.prefactor)
        #print(self.coeff.stem)
        #print(indpos)
        #print(indlabels)
        #exit()
        #print(modes)
        #print(modelabels)
        #exit()

        val_list = gen_index_lists(self.indices,order,eigenvals,symmetry['rot'])
        partitions = []
        expansions = [] #index of the partitions matches index of the summing indices
        nm_constraint = False
        needrefl=False
        #print(symmetry)
        if symmetry['refl']:
            needrefl=True
        needexpand=True
        if needrefl and self.coeff.stem=='r' and self.sinorcos=='sin' and eigenvals[2]==0:
            needexpand=False
        if needrefl and self.coeff.stem=='i' and self.sinorcos=='cos' and eigenvals[2]==0:
            needexpand=False
        if needrefl and self.coeff.stem=='i' and self.sinorcos=='neither' and eigenvals[2]==0:
            needexpand=False
        
        #vals_iterprod = list(itertools.product(*val_list))
        if needexpand:
            for c,n in enumerate(itertools.product(*val_list)):
                meets_real_req = True
                nabs=[abs(n[i]) for i in range(len(n))]
                if (sum(nabs) == order) and (n not in partitions):
                    
                    meets_req = apply_constraints_new(n,self.prefactor,symmetry,indpos,indlabels,eigenvals,self.sinorcos)
                    partitions.append(n)              
                    #if meets_req == True:
                    #    #print(n)
                    #    #exit()
                    if meets_req == True and meets_real_req == True:
                        term = deepcopy(self)
                        term.update_indices(n)
                        term.get_order()
                        expansions.append(term)
        self.expansions = expansions

    def expand_new(self,order,eigenvals,symmetry,modes):
        
        #order: order of expansion
        #eigenvals: list of eigenvalue numbers
        
        acceptable=return_acceptable_states(symmetry)
        modelabels=[]
        for k in range(len(modes)):
            for m in range(len(acceptable.full)):
                if acceptable.full[m]==modes[k]:
                    modelabels.append(acceptable.labels[m])
        
        indpos, indlabels = assign_indices(modes,modelabels)
        mult, omult = assign_multipliers(modes,modelabels)
        #omult=gen_real_restraint(omult,self.indices,order,eigenvals,symmetry['rot'])
        #print(mult)
        #print(omult)

        #val_list = gen_index_lists(self.indices,order,eigenvals,symmetry['rot'])
        partitions = []
        expansions = [] #index of the partitions matches index of the summing indices
        nm_constraint = False
        needrefl=False
        #print(symmetry)
        realexpand=False
        if symmetry['refl']:
            needrefl=True
        needexpand=True
        if needrefl and self.coeff.stem=='r' and self.sinorcos=='sin' and eigenvals[2]==0:
            needexpand=False
        if needrefl and self.coeff.stem=='i' and self.sinorcos=='cos' and eigenvals[2]==0:
            needexpand=False
        if needrefl and self.coeff.stem=='i' and self.sinorcos=='neither' and eigenvals[2]==0:
            needexpand=False
        if needrefl and eigenvals[2]==0:
            realexpand=True
        
        expansions=[]
        indices=[]
        #vals_iterprod = list(itertools.product(*val_list))
        if needexpand:
            nmodes=len(self.indices)
            str_builder=-1*np.ones(nmodes,dtype=int)
            prevd=-1
            orderin=copy(order)
            previ=-1*np.ones(nmodes,dtype=int)
            list_of_terms=[]
            self.expand_for_order(nmodes,str_builder,previ,prevd,orderin,mult,omult,symmetry,indpos,indlabels,eigenvals,list_of_terms)
            #exit()
            for i in range(len(list_of_terms)):
                goodterm=True
                if realexpand:
                    goodterm=gen_real_restraint(list_of_terms[i],self.indices,order,eigenvals,symmetry['rot'])
                if (goodterm):
                    term = deepcopy(self)
                    term.update_indices(list_of_terms[i])
                    term.get_order()
                    expansions.append(term) 
                    indices.append(term.coeff.indices)
                    #if meets_req == True and meets_real_req == True:
                    #    term = deepcopy(self)
                    #    term.update_indices(n)
                    #    term.get_order()
                    #    expansions.append(term)
        self.expansions = expansions
        self.indiceslist=indices

    def expand_for_order(self,Nd,str_builder,previ,prevd,order,mult,omult,symm,ip,il,evals,lot):
        if (prevd==Nd-1):
            if (order==0):
                #print(str_builder)
                #print(self.coeff.stem)
                meets_req = apply_constraints_new(str_builder,self.prefactor,symm,ip,il,evals,self.sinorcos)
                
                if (meets_req):
                    #print([self.prefactor,self.coeff.stem,self.sinorcos,str_builder*mult]) 
                    lot.append(np.ndarray.tolist(str_builder*mult))  
                return
            else:
                return
        else:
            newi=np.ones(Nd,dtype=int)
            d=prevd+1 #go up a dimension
            newi=-1*np.ones(Nd,dtype=int)
            newi[0:d+1]=previ[0:d+1]
            for i in range(previ[d]+1,order*omult[d]-order//2*(mult[d]-1)+1):
                if (omult[d]==2):
                    str_builder[d]=(i+1)//2*(-1)**i
                else:
                    str_builder[d]=i
                newi[d]=i
                self.expand_for_order(Nd,str_builder,newi,d,order-abs(str_builder[d]*mult[d]),mult,omult,symm,ip,il,evals,lot)
            d=d-1 #go back a dimension


    #3
    def compile_expansions(self):   #compile expansion term as sympy expr
        collect_cartesian_terms = None
        for func in self.funcs:
            if isinstance(func,Trig):
                trigfunc = func
                if trigfunc.cartesian == True:
                    if len(trigfunc.args) > 1:
                        collect_cartesian_terms = 2
                    else:
                        collect_cartesian_terms = 1
        free_parameters = []
        compiled_expansions = []
        for term in self.expansions:
            term.coeff.compute()
            term.coeff.compile()
            
            for i,func in enumerate(term.funcs):
                func.compute()
                func.compile()
            compiled_term = 1
            for part in term.form:
                if hasattr(part,'sympy_form') == False:
                    compiled_term = compiled_term*part
                else:
                    compiled_term = compiled_term*part.sympy_form
            #if two e modes:
            if collect_cartesian_terms == 2:
                symlist=[]
                symlist.append('x_1')
                symlist.append('y_1')
                for i in range(1,len(trigfunc.args)):
                    symlist.append(f'x_{i+1}')
                    symlist.append(f'y_{i+1}')
                #compiled_term = sym.collect(compiled_term,['x_alpha','y_alpha','x_beta','y_beta'],exact=True)
                compiled_term = sym.collect(compiled_term,symlist,exact=True)
            #else if 1 e mode:
            elif collect_cartesian_terms == 1:
                compiled_term = sym.collect(compiled_term,['x_1','y_1'],exact=True) 
            compiled_expansions.append(compiled_term)

            if term.coeff.sympy_form not in free_parameters and compiled_term != 0:
                free_parameters.append(term.coeff.sympy_form)
        #self.expansions_old=deepcopy(self.expansions)
        self.expansions = compiled_expansions
        self.parameters = free_parameters

    def compile_formula(self):
        deepcopyself = deepcopy(self)
        sym_formula = ''
        if deepcopyself.prefactor != 1:
            prefactor = str(deepcopyself.prefactor)
            prefactor = prefactor.replace('1j','i')
            sym_formula += ' '+prefactor+'*'
        sym_formula += str(deepcopyself.coeff.compile())
        for func in deepcopyself.funcs:
            func.symcomp()
            sym_formula += func.symb
        return sym_formula
    
    def convert_to_cartesian(self):
        for part in self.funcs:
            part.convert_to_cartesian()

    def convert_to_polar(self):
        for part in self.funcs:
            part.convert_to_polar()

    def adapt_to_unimodal(self,modes):

        if modes[0][0] != 'A': #gamma+a to gamma
            del self.coeff.indices[0]
            self.refresh_indices()
            i = 0
            max_i = len(self.funcs) - 1

            while i <= max_i:
                if isinstance(self.funcs[i],Mono):
                    if self.funcs[i].coord == 'z':
                        del self.funcs[i]
                        i -= 1
                        max_i -= 1

                i += 1

        else: #a+a to a
            for c,i in enumerate(self.indices):
                if 'b' in i:
                    del self.coeff.indices[c]
            self.refresh_indices()

            i = 0
            max_i = len(self.funcs) - 1

            while i <= max_i:
                if isinstance(self.funcs[i],Mono):
                    if self.funcs[i].label == '_beta':
                        del self.funcs[i]
                        i -= 1
                        max_i -= 1
                    else:
                        self.funcs[i].label = ''
                i += 1

        self.update_form()

def assign_indices(modes,modelabels):
    #label and count each modes
    posa=[]
    posb=[]
    pose=[]
    modes0=[mode[0] for mode in modes]
    for i in range(len(modes)):
        if modes0[i]=='E':
            pose.append(i)
        if modes0[i]=='A':
            posa.append(i)
        if modes0[i]=='B':
            posb.append(i)
    
    indlabels=[]
    indpos=[]
    j=0
    for i in range(len(posa)):
        indlabels.append(modelabels[posa[i]])
        indpos.append(j)
        j=j+1
    for i in range(len(posb)):
        indlabels.append(modelabels[posb[i]])
        indpos.append(j)
        j=j+1
    for i in range(len(pose)):
        indlabels.append(modelabels[pose[i]])
        indpos.append(j)
        j=j+2 #skip K indices because they don't contribute to constraints outside of order
    return indpos, indlabels

def assign_multipliers(modes,modelabels):
    #label and count each modes
    posa=[]
    posb=[]
    pose=[]
    modes0=[mode[0] for mode in modes]
    tmodes=0
    for i in range(len(modes)):
        if modes0[i]=='E':
            pose.append(i)
            tmodes=tmodes+2
        if modes0[i]=='A':
            posa.append(i)
            tmodes=tmodes+1
        if modes0[i]=='B':
            posb.append(i)
            tmodes=tmodes+1
    

    mult=np.ones((tmodes),dtype=int)
    omult=np.ones((tmodes),dtype=int)
    j=0
    for i in range(len(posa)):
        mult[j]=1
        omult[j]=1
        j=j+1
    for i in range(len(posb)):
        mult[j]=1
        omult[j]=1
        j=j+1
    for i in range(len(pose)):
        mult[j]=1 #M index
        omult[j]=2 #
        mult[j+1]=2 #K index
        omult[j+1]=1
        j=j+2 #skip K indices because they don't contribute to constraints outside of order
    return mult, omult

def rdelta(p,q):
    if p==q:
        return 1
    else:
        return 0

def apply_constraints_new(n,prefactor,symmetry,indpos,indlabels,eigenvals,sinorcos):
    #n:  list of indices
    #prefactor: either +1,-1 or 1j
    #symmetry: dictionary of symmetry properties
    #indpos: list of index positions that matter, i.e. skipping K
    #indlabels: list of labels for each indices with symmetry properties
    #eigenvals: the list of four eigenvalues that define the system
    #sinorcos: states whether expansion function has sin or cos or neither
    kappa=eigenvals[0]
    rot=symmetry['rot']
    s=kappa
    #print(n)
    for i in range(len(indlabels)):
        if indlabels[i][0]=='B':
            s=s+n[indpos[i]]*(rot//2)
        if indlabels[i][0]=='E':
            s=s+n[indpos[i]]*indlabels[i][1]
        #print('rot')
        #print((s,s//rot))
    if (rot!=0):
        if (s//rot)*rot!=s:
            #print('rot not followed')
            return False
    else:
        if (s!=0):
            return False
        

    #print(isinstance(prefactor, complex))
    needrefl=False
    if symmetry['refl']:
        needrefl=True
        if isinstance(prefactor, complex):
            refleig=eigenvals[2]
        else:
            refleig=eigenvals[1]
        #print(refleig)
        if refleig==1:
            refl_start=0
        else: #refleig=-1
            refl_start=1

    
    if needrefl: #is reflection symmetry needed

        if sinorcos=='sin': #phi->-phi flips sign
            refl_start=refl_start+1
        else: #phi->-phi does nothing
            refl_start=refl_start+0
        s=refl_start
        for i in range(len(indlabels)):
            if indlabels[i][0]=='A' or indlabels[i][0]=='B':
                s=s+n[indpos[i]]*rdelta(2,indlabels[i][1])
        #print('refl')
        #print((s,s//2))
        if (s//2)*2!=s:
            #print('reflect not followed')
            #print(sinorcos)
            #print(prefactor)
            return False
    
    if indlabels[0][2]!=0: #inver symmetry needed
        if eigenvals[3]==1:
            inver_start=0
        else:
            inver_start=1
        s=inver_start
        for i in range(len(indlabels)):
            s=s+n[indpos[i]]*rdelta(2,indlabels[i][2])
        #print('inver')
        #print(n)
        #print(inver_start)
        #print((s,s//2))
        if (s//2)*2!=s:
            return False
    return True

        

class Coeff:
    def __init__(self,letterstem,indices):
        self.letter = letterstem[0]
        self.stem = letterstem[1]
        self.indices = indices
        self.indroot = [indi[0:2] for indi in indices]

    def compute(self):
        for i,index_arg in enumerate(self.indices):
            self.indices[i] = sym.sympify(index_arg)

    def compile(self):
        #compile into sympy symbol
        indexstr = '_{'
        for c,s in enumerate(self.indices):
            indexstr+=str(s)
            if c != (len(self.indices) -1):
                indexstr+=','
            else:
                indexstr += '}'
        self.sympy_form = sym.Symbol(self.letter+'^{'+self.stem+'}'+indexstr,real=True)
        return self.sympy_form

    #def sympy_coeff(self):
    #    indexstr = '_{'
    #   for c,s in enumerate(self.indices):
    #        indexstr+=str(s)
    #        if c != (len(self.indices) -1):
    #            indexstr+=','
    #        else:
    #            indexstr += '}'
    #    return sym.Symbol(self.letter+'^{'+self.stem+'}'+indexstr,real=True)

class Mono:
    def __init__(self,coord,label,arg):
        self.coord = coord
        self.label = label
        self.arg = arg
        self.cartesian = False

    def compute(self):
        self.arg = sym.simplify(self.arg)

    def convert_to_cartesian(self):
        if self.coord == 'rho':
            self.coord = 'x,y'
            self.cartesian = True

    def convert_to_polar(self):
        if self.coord == 'x,y':
            self.coord = 'rho'
            self.cartesian = False

    def symcomp(self):
        self.symb = '*'+self.coord+self.label + '**('+ self.arg +')'

    def compile(self):
        if self.cartesian != True:
            self.sympy_form = sym.Symbol(self.coord+self.label,real=True)**(self.arg)
        else:
            expr = (sym.sqrt(sym.Symbol('x'+self.label,real=True)**2 + sym.Symbol('y'+self.label,real=True)**2))**self.arg #(sqrt(x**2 + y**2))
            self.sympy_form = sym.factor(sym.expand(expr)) #sym.factor

class Trig:
    def __init__(self,args):
        self.args = args
        self.cartesian = False
        self.exp= False

    def change_to_exp(self):
        self.exp= True

    def compute(self):
        for i,arg in enumerate(self.args):
            self.args[i] = sym.simplify(arg)

    def convert_to_cartesian(self):
        self.cartesian = True

    def convert_to_polar(self):
        self.cartesian = False

class Sin(Trig):
    def __init__(self,args):
        Trig.__init__(self,args)

    def symcomp(self):
        #if len(self.args) > 1:
        #    self.symb = '*sin(('+self.args[0]+')'+'*phi_alpha +('+self.args[1]+')*phi_beta)'
        #else:
        #    self.symb = '*sin(('+self.args[0]+')*phi)'
        self.symb = '*sin('
        self.symb = self.symb+self.args[0]+f'*phi_{1}'
        for i in range(1,len(self.args)):
            self.symb=self.symb+'+'+self.args[i]+f'*phi_{i+1}'
        self.symb=self.symb+')'

    def compile(self):
        if len(self.args) > 1:
            if self.cartesian != True:
                arginsin=self.args[0]*sym.Symbol(f'phi_{1}',real=True)
                for i in range(1,len(self.args)):
                    arginsin=arginsin+self.args[i]*sym.Symbol(f'phi_{i+1}',real=True)
                if (self.exp==True):
                    self.sympy_form = sym.exp(sym.I*arginsin)
                else:
                    self.sympy_form = sym.sin(arginsin)
                #self.sympy_form = sym.sin(self.args[0]*sym.Symbol('phi_alpha',real=True)+self.args[1]*sym.Symbol('phi_beta',real=True))
            else:
                arginsin=self.args[0]*sym.atan2(sym.Symbol(f'y_{1}',real=True),sym.Symbol(f'x_{1}',real=True))
                for i in range(1,len(self.args)):
                    arginsin=arginsin+self.args[i]*sym.atan2(sym.Symbol(f'y_{i+1}',real=True),sym.Symbol(f'x_{i+1}',real=True))
                expr=sym.sin(arginsin)
                #expr = sym.sin(self.args[0]*sym.atan2(sym.Symbol('y_alpha',real=True),sym.Symbol('x_alpha',real=True))+self.args[1]*sym.atan2(sym.Symbol('y_beta',real=True),sym.Symbol('x_beta',real=True)))
                self.sympy_form = sym.factor(sym.expand(expr,trig=True)) #sym.factor
        else:
            if self.cartesian != True:
                if (self.exp==True):
                    self.sympy_form = sym.exp(-sym.I*self.args[0]*sym.Symbol(f'phi_{1}',real=True))
                else:
                    self.sympy_form = sym.sin(self.args[0]*sym.Symbol(f'phi_{1}',real=True))
            else:
                expr = sym.sin(self.args[0]*sym.atan2(sym.Symbol('y_1',real=True),sym.Symbol('x_1',real=True)))
                self.sympy_form = sym.factor(sym.expand(expr,trig=True)) #sym.factor

class Cos(Trig):
    def __init__(self,args):
        Trig.__init__(self,args)

    def symcomp(self):
        #if len(self.args) > 1:
        #    self.symb = '*cos(('+self.args[0]+')'+'*phi_alpha +('+self.args[1]+')*phi_beta)'
        #else:
        #    self.symb = '*cos(('+self.args[0]+')*phi)'
        self.symb = '*cos('
        self.symb = self.symb+self.args[0]+f'*phi_{1}'
        for i in range(1,len(self.args)):
            self.symb=self.symb+'+'+self.args[i]+f'*phi_{i+1}'
        self.symb=self.symb+')'
            
    def compile(self):
        if len(self.args) > 1:
            if self.cartesian != True:
                argincos=self.args[0]*sym.Symbol(f'phi_{1}',real=True)
                for i in range(1,len(self.args)):
                    argincos=argincos+self.args[i]*sym.Symbol(f'phi_{i+1}',real=True)
                if (self.exp==True):
                    self.sympy_form = sym.exp(sym.I*argincos)
                else:
                    self.sympy_form = sym.cos(argincos)
                #self.sympy_form = sym.cos(self.args[0]*sym.Symbol('phi_alpha',real=True)+self.args[1]*sym.Symbol('phi_beta',real=True))
            else:
                argincos=self.args[0]*sym.atan2(sym.Symbol(f'y_{1}',real=True),sym.Symbol(f'x_{1}',real=True))
                for i in range(1,len(self.args)):
                    argincos=argincos+self.args[i]*sym.atan2(sym.Symbol(f'y_{i+1}',real=True),sym.Symbol(f'x_{i+1}',real=True))
                expr=sym.cos(argincos)
                #expr = sym.cos(self.args[0]*sym.atan2(sym.Symbol('y_alpha',real=True),sym.Symbol('x_alpha',real=True))+self.args[1]*sym.atan2(sym.Symbol('y_beta',real=True),sym.Symbol('x_beta',real=True)))
                self.sympy_form = sym.factor(sym.expand(expr,trig=True)) #sym.factor
        else:
            if self.cartesian != True:
                if (self.exp==True):
                    self.sympy_form = sym.exp(sym.I*self.args[0]*sym.Symbol(f'phi_{1}',real=True))
                else:
                    self.sympy_form = sym.cos(self.args[0]*sym.Symbol(f'phi_{1}',real=True))
            else:
                expr = sym.cos(self.args[0]*sym.atan2(sym.Symbol('y_1',real=True),sym.Symbol('x_1',real=True)))
                self.sympy_form = sym.factor(sym.expand(expr,trig=True)) #sym.factor







def requires_signswap(rotational_eigenval):
    rot_eigenval_arg = copy(rotational_eigenval)
    if sym.im(rot_eigenval_arg).is_negative:
        return True
    else:
        return False    

from vhegen.tables.formulas import requirements_dct, return_formula

def get_root_formula(eigenvals,modes,symmetry):
    #symmetry symmetry dictionary
    #n_arg rotation number
    #modes list of vibrational modes
    #eigenvals, eigenvals from one or two electronic states
    #print(modes)
    #print(eigenvals)
    n_arg= symmetry['rot']
    n = str(n_arg)
    rotational_eigenval = eigenvals[0]
    refl_Re = eigenvals[1]
    refl_Im = eigenvals[2]
    inver=eigenvals[3]
    searchmodes = [mode[0] for mode in modes]

    #if len(modes) < 2:
    #    searchmodes.append('A')
    #    searchmodes = reorder(searchmodes)
    #if requires_signswap(rotational_eigenval) == True:
    #    reqs = [sym.re(rotational_eigenval) - sym.im(rotational_eigenval)*1j,searchmodes]
    #else:
    #    reqs = [rotational_eigenval,searchmodes]
    formula=return_formula(n,modes,symmetry)
    return(formula)
    #for k in requirements_dct:
    #    #print(reqs)
    #    #print(requirements_dct[k])
    #    if requirements_dct[k] == reqs:
    #        formula = return_formula(n,k)
    #        #print(formula)
    #        #exit()
    #        if requires_signswap == True:
    #            for term in formula:
    #                print(term)
    #                if isinstance(term,ET):
    #                    term.prefactor = sym.conjugate(term.prefactor)     
    #        #for term in formula:
    #        #    print(term.indices)     
    #        #exit()
    #        if refl_Im == 0:
    #            max_c = len(formula) - 1
    #            c = 0
    #            while c <= max_c:
    #                if isinstance(formula[c],Term):
    #                    if type(formula[c].prefactor) == complex:
    #                        del formula[c]
    #                        c -= 1
    #                        max_c -= 1
    #                c += 1
    #        return(formula)
    raise Exception('VHEGENError: Could not find root formula for rotational eigenvalue '+str(rotational_eigenval)+' and modes '+str(modes)+'.')

def read_index_arg(expr):
    #expr one type of index
    step = 1
    if expr[0] == 'K':
        step=2
    return step

def pos_ints(max,step):
    return list(range(0,max+1,step))

def any_ints(max):
    return list(range(-max,max+1))

def gen_index_attrs(expr,max_val):
    #expr: one index in list of indices
    #max_val: order of expansion (i.e. sum of indice values)
    step= read_index_arg(expr)
    if 'M' in expr:
        vals = any_ints(max_val)
    else:
        vals = pos_ints(max_val,step)

    return vals

def gen_index_lists(indices,order,eigenvals,nrot):
    val_list = []
    native_par_list = []
    for i in indices: #regular case
            vals = gen_index_attrs(i,order)
            val_list.append(vals)

    if eigenvals[2] == 0: #take real case
        m_found = False
        for c,i in enumerate(indices):
            if 'M' in i:
                m_found = True
                m_index = c #only need to keep first m positive
                break

        if eigenvals[0] == 0:
            if m_found == True:
                c = 0
                while c <= (len(val_list[m_index]) - 1):
                    if val_list[m_index][c] < 0:
                        del val_list[m_index][c]
                        c -= 1
                    c += 1
        elif abs(eigenvals[0]) == nrot//2:
            if m_found == True:
                c = 0
                while c <= (len(val_list[m_index]) - 1):
                    if val_list[m_index][c] < 0:
                        del val_list[m_index][c]
                        c -= 1
                    c += 1
        #print(eigenvals)
        #print(val_list)
        #exit()
    return val_list

def gen_real_restraint(list_of_terms,indices,order,eigenvals,nrot):

    if eigenvals[2] == 0: #take real case
        for c,i in enumerate(indices):
            #print(list_of_terms)
            #print(i,list_of_terms[c])
            if 'M' in i and list_of_terms[c]<0:
                return False
            if 'M' in i and list_of_terms[c]>0:
                return True
    return True
        
def adapt_to_unimodal(formula,modes):
    adapted_formula = []
    if formula != [0]:
        for term in formula:
            term.adapt_to_unimodal(modes)
            adapted_formula.append(term)
        
    return adapted_formula

def get_sym_props(pointgroup):
    sym_props = {'inver': False,
                 'refl': False,
                 'hrefl':False}
    if ('D' in pointgroup) or ('V' in pointgroup):
        sym_props['refl'] = True
    if 'H' in pointgroup:
        if '4' in pointgroup:
            sym_props['inver'] = True
        elif '3' in pointgroup:
            sym_props['hrefl'] = True
    return sym_props

def get_constraint(eigenvals, vib_modes, operation):
    princ_rot = eigenvals[0]
    refl_Re, refl_Im = eigenvals[1], eigenvals[2]
    inver_eigenval = eigenvals[3]
    try:
        os.chdir(os.getcwd() + '/constraints')
        constraintfiles = os.listdir(os.getcwd())
        constraintargs = [sym.sympify(i[0:i.index('_')].replace('X','*')) for i in constraintfiles]
        #Open refl.sym file
        for arg in constraintargs:
            if sym.simplify(princ_rot-arg) == 0:
                with open(replace_all(str(arg),{'*':'X',' ':''})+"_"+operation+".sym","r") as sym_file:
                    constraint_lines  = sym_file.readlines()
            elif sym.simplify(sym.re(princ_rot)-sym.im(princ_rot)*1j - arg) == 0:
                with open(replace_all(str(arg),{'*':'X',' ':''})+"_"+operation+".sym","r") as sym_file:
                    constraint_lines  = sym_file.readlines()
        os.chdir('..')
    except OSError as e:
        print(e)
        
    vib_modes = [i.replace("''",'"') for i in vib_modes]

    if operation == 'refl':
        eigen_req = '[' + str(refl_Re) + ',' + str(refl_Im) + ']'
        vib_modes = [replace_all(i,{'G':'','U':'',"'":'','"':''}) for i in vib_modes]

    elif operation == 'inver' or operation == 'hrefl':
        eigen_req = '[' + str(inver_eigenval) + ']'
        vib_modes = [replace_all(i,{'1':'','2':''}) for i in vib_modes]

    if len(vib_modes) == 1:        
        if operation == 'refl':
            vib_modes.append('A1')
        elif operation == 'inver':
            vib_modes.append('AG')
        elif operation == 'hrefl':
            vib_modes.append("A'")
    mode_req = vib_modes
    max_count = len(constraint_lines) - 1
    for count,line in enumerate(constraint_lines):
        modes_in_line = ''.join(line.split(',')[0])
        modes_in_line = modes_in_line[1:-1]
        modes_in_line = modes_in_line.split('.')
        line = line.replace('\n','')  
        if set([str(mode) for mode in mode_req]) == set(modes_in_line) and eigen_req in line:
            match_line = line
        if (count == max_count):
            try:
                match_line
            except NameError: #no matching constraints
                return {}
    constraints = match_line.split(': ',1)[1]

    constraints = constraints.split(',')
    constraints_dict = {}
    for constraint in constraints:
        index = constraint.split(' ')[0]
        restriction = ' '.join(constraint.split(' ')[1:])
        constraints_dict[index] = restriction
    if 'all' in constraints_dict:
        if constraints_dict['all'] == 'nr':
            del constraints_dict['all'] #check if any composite cases where 
    return constraints_dict
    
def load_constraints(sym, eigenvals, vib_modes):
    sym_props = get_sym_props(sym)
    matrix_element_constraints = {}
    for e in eigenvals:
        constraints_dict = {}
        if sym_props['refl'] == True:
            constraints_refl = get_constraint(eigenvals[e], vib_modes, 'refl')
            for key in constraints_refl:
                if key not in constraints_dict:
                    constraints_dict[key] = constraints_refl[key]
                else:
                    constraints_dict[key] += constraints_refl[key]
        if sym_props['inver'] == True:
            constraints_inver = get_constraint(eigenvals[e], vib_modes, 'inver')
            for key in constraints_inver:
                if key not in constraints_dict:
                    constraints_dict[key] = constraints_inver[key]
                else:
                    constraints_dict[key] += constraints_inver[key]
        if sym_props['hrefl'] == True:
            constraints_hrefl = get_constraint(eigenvals[e], vib_modes, 'hrefl')
            for key in constraints_hrefl:
                if key not in constraints_dict:
                    constraints_dict[key] = constraints_hrefl[key]
                else:
                    constraints_dict[key] += constraints_hrefl[key]
        matrix_element_constraints[e] = constraints_dict
    return matrix_element_constraints

def apply_constraints(local_constraints, p, fitted_term,unfitted_term): #returns whether constraints are sat (True) or not sat (False)
    sat_constraints = 0
    if local_constraints != {}:
        if 'all' in local_constraints.keys():
            if 'na' in local_constraints['all']: #all na case
                return False #return not sat
            else:
                if len(local_constraints) == 1: #all nr case. len 1 check ensures theres no other constraints to consider.
                    return True #return sat

        for c in local_constraints:
            #strip whitespace escape sequences
            if isinstance(local_constraints[c],set):
                local_constraints[c] = set([i.strip() for i in local_constraints[c]])
            else:
                local_constraints[c] = local_constraints[c].strip()
            
            #local_constraints[c] = local_constraints[c].strip() #remove whitespace 
            if type(c) == int: #handle summing-index based constraints
                if type(local_constraints[c]) == set:
                    len_cond = len(local_constraints[c])
                    met_cond = 0
                    for i in local_constraints[c]:
                        if constraint_funcs[i](c,p) == True:
                            met_cond += 1
                    if met_cond == len_cond:
                        sat_constraints += 1
                elif constraint_funcs[local_constraints[c]](c,p) == True:
                    sat_constraints += 1

            elif type(c.split('&')[0]) == int: #handle pairwise summing-index based constraints
                if constraint_funcs[local_constraints[c]]([int(i) for i in c.split('&')],p) == True:
                    sat_constraints += 1

            elif 'nz' in local_constraints[c]: #handle nz constraints
                if 'if' in local_constraints[c]: #has condition
                    cond = (local_constraints[c].split(' '))[-2:]
                else:
                    cond = []
                if '&' in c: #pairwise nz
                    c1 = c.split('&')[0]
                    c2 = c.split('&')[1]
                    if nz(c1,cond,fitted_term,unfitted_term,p) == True or nz(c2,cond,fitted_term,unfitted_term,p) == True:
                        sat_constraints +=1
                else:
                    if nz(c,cond,fitted_term,unfitted_term,p) == True:
                        sat_constraints += 1
    return sat_constraints == len(local_constraints)

#Constraint functions
def even(i,p):
    if p[i] % 2 != 0:
        return False
    else:
        return True

def odd(i,p):
    if p[i] % 2 == 0:
        return False
    else:
        return True

def ee_or_oo(i,p):
    if (p[i[0]] % 2 + p[i[1]] % 2) == 1:
        return False
    else:
        return True

def eo_or_oe(i,p):
    if (p[i[0]] % 2 + p[i[1]] % 2) != 1:
        return False
    else:
        return True

def non_neg(i,p):
    if p[i] >= 0:
        return True
    else:
        return False

def nm_postproc(i,p):
    if p[i][0] >  0:
        return True 
    if p[i][0] == 0 and p[i][1] >= 0:
        return True
    else:
        return False 

def na(i,p):
    return False

def nr(i,p):
    return True

def nz(arg,cond,fitted_term,unfitted_term,p): #coeff nz and trig nz constraints
    meets_req = False
    satisfy_arg = False
    if cond == []:
        meets_cond = True
    else: #if condition isn't true, constraint need not be applied.
        for c,i in enumerate(unfitted_term.indices):
            if cond[0] in i:
                cond_index = c
        if cond[1] == 'even':
            if p[cond_index] % 2 == 0:
                meets_cond = True
            else:
                meets_cond = False 
        else: #odd
            if p[cond_index] % 2 != 0:
                meets_cond = True
            else:
                meets_cond = False
    if meets_cond == False: #bypass
        meets_req = True
        return meets_req

    else: #nz constrain
        if arg in ['sin','cos']: #handle trig nz
            for part in fitted_term.funcs:
                if (arg == 'sin') and (isinstance(part,Sin)):
                    satisfy_arg = True
                elif (arg == 'cos') and (isinstance(part,Cos)):
                    satisfy_arg = True
        else: #handle coeff nz 
            satisfy_arg = coeff_nz(arg,fitted_term)
    if satisfy_arg == True:
        meets_req = True
    else:
        meets_req = False
    return meets_req

def coeff_nz(arg,term):
    validate = 0
    if arg[1] != '%':
        if term.coeff.stem == arg[1]:
            validate += 1
    else:
        validate += 1
    if not arg.endswith(R'%%'):
        index_parities = [arg[2],arg[3]]
        for c,i in enumerate(index_parities):
            if i == 'e':
                if int(term.indices[c]) % 2 == 0:
                    validate += 1
            elif i == 'o':
                if int(term.indices[c]) % 2 != 0:
                    validate += 1
    else:
        validate += 2
    if validate == 3:
        return True
    else:
        return False

def get_symbolic_formula(formula):
    symb_formula = ''
    if formula != [0]:
        for c,term_formula in enumerate(formula):
            symb_term_formula = term_formula.compile_formula()
            if symb_term_formula[1] != '-' and c != 0:
                symb_formula += '+'
            symb_formula += term_formula.compile_formula()
    return symb_formula

def get_matrix_element_expansion(formula,order,e_coord_system,eigenvals,symmetry,modes):
    
    expansions = []
    params = []
    
    if e_coord_system == 'both':
        iterates = 2
    else:
        iterates = 1
    for i in range(iterates):
        if (i==0):
            print('generating polar expansion')
        else:
            print('generating cartesian expansion')
        total_expansion = []
        if formula != [0]:
            for f,term_formula in enumerate(formula):
                #term_formula.grab_constraints(constraints)
                if e_coord_system == 'pol' or i == 0:
                    term_formula.convert_to_polar()
                if e_coord_system == 'cart' or i == 1:
                    term_formula.convert_to_cartesian()
                #if e_coord_system == 'expandedpol' or i == 2:
                #    term_formula.convert_to_expandedpol()
                #term_formula.expand(order,eigenvals,symmetry,modes)
                print(f'generating expansion for term {f}')
                term_formula.expand_new(order,eigenvals,symmetry,modes)
                
                #print(term_formula.expansions)   
                #print(term_formula.parameters)   
            for f,term_formula in enumerate(formula):
                if e_coord_system=='pol' or i==0:
                    if (len(formula)==4 and f<2):
                        for pos1 in range(len(term_formula.expansions)):
                            #print(pos1)
                            #print(term_formula.indiceslist[pos1])
                            if term_formula.indiceslist[pos1] in formula[f+2].indiceslist:
                                term_formula.expansions[pos1].change_to_exp()
                                pos2=formula[f+2].indiceslist.index(term_formula.indiceslist[pos1])
                            #    #print(term_formula.expansions[pos1])
                            #    #print(formula[f+2].expansions[pos2])
                            #    #term_formula.expansions[pos1]=(term_formula.expansions[pos1]+formula[f+2].expansions[pos2]).simplify()
                                del formula[f+2].expansions[pos2]#.prefactor=sym.S.Zero
                                del formula[f+2].indiceslist[pos2]
                            #    #print(term_formula.expansions[pos1])
                print(f'compiling term {f}')
                term_formula.compile_expansions()
                total_expansion += term_formula.expansions
                params += term_formula.parameters
        expansions.append(total_expansion)
        #print(' ')
    params = set(params)
    built_expansions = [0]*iterates
    for c,i in enumerate(expansions):
        for exp in i:
            built_expansions[c] += exp
    
    return built_expansions,params

constraint_funcs = {'even':even,
                    'odd':odd,
                    'ee_or_oo':ee_or_oo,
                    'eo_or_oe':eo_or_oe,
                    'non_neg': non_neg,
                    'nm_postproc': nm_postproc,
                    'na':na,
                    'nr':nr}

#EOF