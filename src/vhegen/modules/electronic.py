from sympy import Matrix, Symbol, Add, Mul, UnevaluatedExpr, conjugate, re, im, printing, simplify, sqrt, sympify, expand_complex
from sympy.physics.quantum import Ket, Bra
from vhegen.modules.glbls import state_components_dct, replace_all
from copy import copy,deepcopy

def kdelta(p,q):
    if p==q:
        return 1
    else:
        return 0

def refl_parity(state):
    if ('2' in state):
        return 2
    elif ('1' in state):
        return 1
    else:
        return 0

def inver_parity(state):
    if ('U' in state) or ("''" in state):
        return -1
    elif ('G' in state) or ("'" in state):
        return 1
    else:
        return 0

def get_state_components(states):

    label_dct = {0: R'_\alpha',
                 1: R'_\beta'}

    labelling = False
    state_types = [s[0] for s in states]
    if len(state_types) == 2:
        if state_types[0] == state_types[1]: #Gamma + Gamma 
            labelling = True

    state_components = []
    for c,s in enumerate(state_types):
        if labelling == True:
            state_components += [s_i+label_dct[c] for s_i in state_components_dct[s]]
        else:
            state_components += state_components_dct[s]
    
    return state_components
            
class VibronicMatrix:
    def __init__(self,symmetry,states):

        from vhegen.tables.matrices import requirements_dct, matrix_dct, return_matrices
        from vhegen.tables.eigenvalues import return_eigenvals

        self.original_mat, self.symmetrized_mat, self.transform_mat = return_matrices(symmetry,states)
        self.symmetry_eigenvals = return_eigenvals(symmetry,states) #Matrix element symmetry eigenvalues
        #print(self.symmetry_eigenvals)
        #exit()
        self.state_components = get_state_components(states) #Matrix ket-bra basis representation
        self.get_matrix_product_form()

    def get_matrix_product_form(self):
        self.ket_row = Matrix([Ket(s) for s in self.state_components]).transpose()
        self.bra_col = Matrix([Bra(s) for s in self.state_components])
        self.form = UnevaluatedExpr(self.ket_row)*UnevaluatedExpr(self.symmetrized_mat)*UnevaluatedExpr(self.bra_col)
        return self.form

    def __str__(self):
        return str(self.symmetrized_mat)

    def get_dependencies(self):

        independent_set = set()
        for m_e in self.symmetrized_mat:
            if not isinstance(m_e,conjugate) and m_e != 0:
                independent_set.add(m_e)

        self.dependencies = {}
        for indep_e in independent_set:
            self.dependencies[indep_e] = []
            for c,e in enumerate(self.symmetrized_mat):
                if e == indep_e:
                    self.dependencies[indep_e].append(self.original_mat[c])
                if conjugate(e) == indep_e:
                    self.dependencies[indep_e].append(conjugate(self.original_mat[c]))

        return self.dependencies

    def format_TeX_mat(self):
        tex_mat = deepcopy(self.original_mat)
        for c,i in enumerate(tex_mat):
            tex_mat[c] = format_matrix_element(i)

        tex_form = UnevaluatedExpr(self.ket_row)*UnevaluatedExpr(tex_mat)*UnevaluatedExpr(self.bra_col)
        tex_formatted_form = R'$\hat{H}='+printing.latex(tex_form,mat_delim='(',mat_str='matrix')+'$'
            
        return tex_formatted_form

    def returnpos(self):
        self.element_position={}
        scrublist=[]
        l=len(self.original_mat)
        for i in range(1,5):
            if (i*i==l):
                l=i
                break
        self.element_position['matsize']=l
        for c,e in enumerate(self.original_mat):
            if e != 0:
                print(c)
                print(e)
                cpos=(c)//l+1
                rpos=(c-(cpos-1)*l)+1
                pos=[cpos,rpos]
                if e not in scrublist:
                    scrublist.append(e)
                self.element_position[e]=pos
        #print(self.element_position)

    def change_basis(self): #Complex to real matrix transformation

        symmetrized_mat_expanded = deepcopy(self.symmetrized_mat)
        for c,e in enumerate(symmetrized_mat_expanded): #Expand all elements into real and imaginary parts
            if type(e) == conjugate:
                symmetrized_mat_expanded[c] = re(conjugate(e)) - 1j*im(conjugate(e))
            else:
                symmetrized_mat_expanded[c] = re(e) + 1j*im(e)

        evolved_mat = self.transform_mat*symmetrized_mat_expanded*(self.transform_mat.H)

        #Update matrix product form:
        for c in range(0,len(self.state_components)):
            self.state_components[c] = replace_all(self.state_components[c], {'+':'X','-':'Y'})
        
        self.element_position={}
        scrublist=[]
        l=len(self.original_mat)
        for i in range(1,5):
            if (i*i==l):
                l=i
                break
        self.element_position['matsize']=l
        #print('converting')
        for c,e in enumerate(self.original_mat):
            evolved_mat[c] = simplify(expand_complex(evolved_mat[c])) #Simplify evolved mat element-wise
            if e != 0:
                e = replace_all(str(e), {'+':'X','-':'Y'})
                #print(c)
                cpos=(c)//l+1
                rpos=(c-(cpos-1)*l)+1
                pos=[cpos,rpos]
                if e not in scrublist:
                    scrublist.append(e)
                self.element_position[Symbol(e)]=pos
                self.original_mat[c] = Symbol(e)

        #print(self.element_position)
        #print(self.original_mat)
        self.evolved_mat = evolved_mat
        self.get_matrix_product_form()

        #Element mapping
        mapping = {}
        for c,e in enumerate(self.original_mat):
            if e != 0:
                e_str = replace_all(str(e), {'+':'X','-':'Y'})
                mapping[Symbol(e_str)] = self.evolved_mat[c]

        return mapping
        
def get_dependent_elements(expansion_dct,count,obj):
    #expansion_dct: dictionary of expansions with label o for order
    #count is dictionary of number of terms for each order o
    #obj vhegen class
    print('\nApplying expansion dependencies...\n')
    all_expansions_dct = {} 
    dependent_count = {}
    for o in expansion_dct: #for order
        all_expansions_dct[o] = {}
        dependent_count[o] = {}
        for e in expansion_dct[o]:
            for dep_e in obj.matrix.dependencies[e]:
                if isinstance(dep_e,conjugate):
                    all_expansions_dct[o][conjugate(dep_e)] = [conjugate(expansion_dct[o][e][c]) for c,i in enumerate(expansion_dct[o][e])]
                    dependent_count[o][conjugate(dep_e)] = count[o][e]
                else:
                    all_expansions_dct[o][dep_e] = expansion_dct[o][e]
                    dependent_count[o][dep_e] = count[o][e]
    return all_expansions_dct

def format_matrix_element(m_e):
    if m_e != 0:
        if isinstance(m_e,conjugate):
            m_e_str = replace_all(str(conjugate(m_e)), {'__':'', 'alpha':R'{\alpha}', 'beta':R'{\beta}'})
            return (Symbol(R'{H_{'+str(m_e_str)+R'}}^*'))
        else:
            m_e_str = replace_all(str(m_e), {'__':'', 'alpha':R'{\alpha}', 'beta':R'{\beta}'})
            return Symbol(R'H_{'+m_e_str+R'}')
    else: return 0

def prune_dependent_elements(real_matrix_elements):
    indep_m_e = []
    m_e = [str(i) for i in real_matrix_elements]

    if 'X__A' in m_e or 'X__B':
        for i in m_e:
            if i[0] == 'X' or i[0] == 'Y':
                indep_m_e.append(Symbol(i))

    elif 'X__X' in m_e:
        for i in m_e:
            if i[0] == 'X':
                indep_m_e.append(Symbol(i))

    elif '__' in m_e[0]:
        m_e_lists = []
        for i in m_e:
            m_e_lists.append(i.split('__'))
        for i in m_e_lists:
            if 'alpha' in i[0] and 'X' in i[0]:
                indep_m_e.append(Symbol('__'.join(i)))

    
    else:
        indep_m_e = real_matrix_elements

    return indep_m_e

def map_elements(expansions, mapping):
    coords = 1
    if len(list(list(expansions.values())[0].values())[0]) > 1:
        coords = 2
    mapped_expansions = {}
    term_inheritance = {}
    independent_elements = prune_dependent_elements([i for i in mapping])
    for o in expansions:
        sub_tuples = [(e,expansions[o][e][0]) for e in expansions[o].keys()]
        if coords == 2:
            sub_tuples_2 = [(e,expansions[o][e][1]) for e in expansions[o].keys()]
        mapped_expansions[o] = {}
        term_inheritance[o] = {}
        for k in mapping:
            get_inheritance = False
            if k in independent_elements:
                get_inheritance = True
                inheritance_count = {}
            expanded_terms = []
            if type(mapping[k]) == Add:
                k_args = list(mapping[k].args)
            else:
                k_args = [mapping[k]]
            k_args = (k_args, [i.free_symbols for i in k_args])
            for c,element in enumerate(k_args[1]):
                for c2,item in enumerate(sub_tuples):
                    if item[0].free_symbols == element:
                        if coords == 1:
                            substituted_args = [k_args[0][c].subs(item[0],item[1])]
                        else:
                            substituted_args = [k_args[0][c].subs(item[0],item[1]),k_args[0][c].subs(sub_tuples_2[c2][0],sub_tuples_2[c2][1])]
                        expanded_terms.append(substituted_args)
                        term_arguments = substituted_args[0]

                        if type(term_arguments) == Symbol:
                            term_count = 1
                        elif type(term_arguments) == Mul: #handle numeric prefactor and single terms
                            sum_exists = False
                            for part in term_arguments.args:
                                if type(part) == Add:
                                    sum_exists = True
                                    term_count = len(part.args)
                            if sum_exists == False:
                                term_count = 1
                        else: #sym.Add
                            try:
                                term_count = len(substituted_args[0].args)
                            except AttributeError:
                                term_count = 0
                        if get_inheritance == True:
                            inheritance_count[k_args[1][c].pop()] = term_count

            sum_terms = 0
            for term in expanded_terms:
                sum_terms += term[0]
            if coords == 2:
                sum_terms2 = 0
                for term in expanded_terms:
                    sum_terms2 += term[1]
                sum_terms = [sum_terms,sum_terms2]
            else:
                sum_terms = [sum_terms]
            mapped_expansions[o][k] = sum_terms
            if get_inheritance == True:
                term_inheritance[o][k] = inheritance_count
    return mapped_expansions, term_inheritance
    
#EOF
