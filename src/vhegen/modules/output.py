import sympy as sym
from sympy.printing.latex import LatexPrinter, print_latex
import os
import subprocess
import vhegen.modules.electronic as el

log_str = ''

def log_append(string):
    #add contents to log file
    global log_str
    log_str += str(string) + '\n'

def format_mulliken_scripts(irrep,t):
    if t == 's':
        irrep_str = irrep[0].upper()
    else:
        irrep_str = irrep[0].lower()
    primes = 0
    for char in irrep:
        if char == "'":
            primes += 1
    if primes == 0:
        irrep_str += '_{' + irrep[1:].lower() + '}'
    else:
        irrep_str += primes*"'" + '_{' + irrep[1:].replace("'","") + '}'
    return irrep_str

def problem_format2TeX(states,modes):
    problem_str = '$'
    if len(states) > 1:
        problem_str += '('
    for c,s in enumerate(states):
        if c>0:
            problem_str +='+'
        problem_str += format_mulliken_scripts(s,'s')
    if len(states) > 1:
        problem_str += ')'
    problem_str += R' \otimes '
    if len(modes) > 1:
        problem_str += '('
    for c2,m in enumerate(modes):
        if c2>0:
            problem_str+='+'
        problem_str += format_mulliken_scripts(m,'m')
    if len(modes) > 1:
        problem_str += ')'
    problem_str += '$'
    return problem_str

def prune_dependent_terms(count):
  keys = list(count.keys())
  unique_keys = []
  unique_count = {}
  for i in keys:
    if sym.Symbol(str(i)[::-1]) not in unique_keys:
      unique_keys.append(i)
  for k in unique_keys:
    unique_count[k] = count[k]
  return(unique_count)

def count_format2TeX(count): 
  #{+A: 8, -A: 8, A+: 8, A-: 8}
  count = prune_dependent_terms(count)
  count_str = ''
  for c,m_e in enumerate(count):
    count_str += '$'+str(el.format_matrix_element(m_e))+'$: ' + '$'+str(count[m_e]) +'$'
    if c != len(count)-1:
      count_str += ', '
    else:
      count_str += '.'
  return count_str

def realcount_format2TeX(count): #includes inheritance from complex
  #{XA: {+A: 6, ++: 10}, YA: {+A: 6}, AX: {+A: 6}, AY: {+A: 6}}
  count = prune_dependent_terms(count)
  count_str = ''
  for c,m_e in enumerate(count):
    sumval = sum(count[m_e][i] for i in count[m_e])
    count_str += '$'+str(el.format_matrix_element(m_e))+'$: ' + '$'+str(sumval) +'$'
    if sumval != 0:
      count_str += ' ('
      for c2,inherited_term in enumerate(count[m_e]):
        added= False
        if sumval == count[m_e][inherited_term]:
          count_str += 'all from $'+str(el.format_matrix_element(inherited_term))+'$)'
          break
        else:
          if count[m_e][inherited_term] != 0:
            count_str += '$'+str(count[m_e][inherited_term])+'$ from $'+str(el.format_matrix_element(inherited_term))+'$'
            added=True

        if c2 != len(count[m_e])-1:
          if added==True:
            count_str += ', '
        else:
          count_str += ')'

    if c!= len(count)-1:
      count_str += ', '
    else:
      count_str += '.'
  return count_str 

def compose_TeX_both_bases(sym,modes,TeX_matrix_complex,TeX_matrix_real,TeX_problem,TeX_expansions_complex,TeX_expansions_real):
    #print(TeX_expansions_complex)
    from vhegen.modules.input import return_acceptable_states
    amodes=return_acceptable_states(sym)
    modelabels=[]
    for i in range(len(modes)):
       for j in range(len(amodes.full)):
          if modes[i]==amodes.full[j]:
             modelabels.append(amodes.labels[j])
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
    indices=[]
    for i in range(numa):
        indices.append(Rf'z_{i+1} &\rightarrow& '+format_mulliken_scripts(modes[posa[i]],'m')+'\\\\ \n')
    for i in range(numb):
        indices.append(Rf'w_{i+1} &\rightarrow& '+format_mulliken_scripts(modes[posb[i]],'m')+'\\\\ \n')
    for i in range(nume):
        indices.append(Rf'\rho_{i+1},\phi_{i+1},x_{i+1},y_{i+1} &\rightarrow& '+format_mulliken_scripts(modes[pose[i]],'m')+'\\\\ \n')
        #indices.append(Rf'\phi_{i+1} &\rightarrow& '+format_mulliken_scripts(modes[pose[i]],'m')+'\\\\ \n')
    latexscript=R'\begin{eqnarray*}'+'\n'
    for i in range(len(indices)):
        latexscript=latexscript+indices[i]
    latexscript=latexscript+'\end{eqnarray*}'
    #print(latexscript)
    #exit()
    TeX_str = (R'\batchmode'+'\n'
           R'\documentclass[fleqn]{article}'+'\n'
           R'\usepackage{amsmath}'+'\n'
           #R'\usepackage{breqn}'+'\n'
           R'\usepackage{hyperref}'+'\n'
           R'\usepackage[margin=0.7in]{geometry}' + '\n'
           R'\setlength\mathindent{0pt}'+'\n'
           R'\allowdisplaybreaks'+'\n'
           R'\title{VHEGEN: A vibronic Hamiltonian expansion generator for trigonal and tetragonal polyatomic systems}'+'\n'
           R'\author{James Brown \and Robert A. Lang \and Riley J. Hickman \and Tao Zeng}'+'\n'
           R'\date{}'+'\n'
           R'\begin{document}'+'\n'
           R'\maketitle'+'\n'
           R'Thank you for using \texttt{VHEGEN}, the \texttt{V}-ibronic \texttt{H}-amiltonian \texttt{E}-xpansion \texttt{GEN}-erator for trigonal and tetragonal polyatomic systems. '
           R'This is a \texttt{VHEGEN} output file compiled by \texttt{pdflatex}. '
           R'If the \texttt{VHEGEN} package was used in research resulting in a publication, please reference the article in \textit{Computer Physics Communications} which describes the program ([doi here]). '
           R'Additional information regarding the matrix element expansion process, including the independent matrix element eigenvalues, their root formulas and constraints, and their transformation to the real basis (if applicable), can be found in the \texttt{log} output file. '
           R'For questions, bugs, or comments, please contact jbrown88@yorku.ca.\\\\'+'\n'
           R'\tableofcontents'+'\n'
           R'\newpage'+'\n'
           R'\section{Vibronic interaction}'+'\n'
           R''+TeX_problem+' in $'+sym['letter']+'_{'+sym['print'][1:].lower()+R'}$'+'\n'
           R''+latexscript+'\n'
           R'\section{Vibronic Hamiltonian operator in the complex $E$ basis}'+'\n'
           R''+TeX_matrix_complex.replace('__','')+'\n'
           R'\section{Matrix element expansions in the complex $E$ basis}'+'\n'
           R''+TeX_expansions_complex.replace('__','')+'\n'
           R'\section{Vibronic Hamiltonian operator in the real $E$ basis}'+'\n'
           R''+TeX_matrix_real.replace('__','')+'\n'
           R'\section{Matrix element expansions in the real $E$ basis}'+'\n'
           R''+TeX_expansions_real.replace('__','')+'\n'
           R'\end{document}')
    return TeX_str

def compose_TeX(sym,modes,TeX_matrix,TeX_problem,TeX_expansions):
    


    TeX_str = (R'\batchmode'+'\n'
               R'\documentclass[fleqn]{article}'+'\n'
               R'\usepackage{amsmath}'+'\n'
               #R'\usepackage{breqn}'+'\n'
               R'\usepackage[margin=0.7in]{geometry}' + '\n'
               R'\setlength\mathindent{0pt}'+'\n'
               R'\allowdisplaybreaks'+'\n'
               R'\usepackage{hyperref}'+'\n'
               R'\title{VHEGEN: A vibronic Hamiltonian expansion generator for trigonal and tetragonal polyatomic systems}'+'\n'
               R'\author{James Brown \and Robert A. Lang \and Riley J. Hickman \and Tao Zeng}'+'\n'
               R'\date{}'+'\n'
               R'\begin{document}'+'\n'
               R'\maketitle'+'\n'
               R'Thank you for using \texttt{VHEGEN}, the \texttt{V}-ibronic \texttt{H}-amiltonian \texttt{E}-xpansion \texttt{GEN}-erator for trigonal and tetragonal polyatomic systems. '
               R'This is a \texttt{VHEGEN} output file compiled by \texttt{pdflatex}. '
               R'If the \texttt{VHEGEN} package was used in research resulting in a publication, please reference the article in \textit{Computer Physics Communications} which describes the program. '
               R'Additional information regarding the matrix element expansion process, including the independent matrix element eigenvalues, their root formulas and constraints, and their transformation to the real basis (if applicable), can be found in the \texttt{log} output file. '
               R'For questions, bugs, or comments, please contact jbrown88@yorku.ca.\\\\'+'\n'
               R'\tableofcontents'+'\n'
               R'\newpage'+'\n'
               R'\section{Vibronic interaction}'+'\n'
               R''+TeX_problem+' in $'+sym['letter']+'_{'+sym['print'][1:].lower()+R'}$'+'\n'
               #R''+latexscript+'\n'
               R'\section{Vibronic Hamiltonian operator}'+'\n'
               R''+TeX_matrix.replace('__','')+'\n'
               R'\section{Matrix element expansions}'+'\n'
               R''+TeX_expansions.replace('__','')+'\n'
               R'\end{document}')
    return TeX_str

def format_constraint(constraint): #{'a%ee': 'nz','broe&bieo': 'nz'}
    string = ''
    if constraint == {}:
        return ''
    for c,con in enumerate(constraint):
        cons = con.split('&')
        formatted_cons = []
        for const in cons:
            formatted_con = const
            if 'nz' in constraint[con] and con not in ['cos','sin']:
                add = 0
                if const[1] != '%':
                    formatted_con = const[:1] + '^' + const[1:]
                    add += 1
                if const[2] != '%':
                    formatted_con = formatted_con[:(2+add)] + '_' + formatted_con[(add+2):]
            formatted_cons.append(formatted_con)
        formatted_constraints = ' '.join(formatted_cons)
        formatted_constraints = formatted_constraints.replace('%','')
        string += formatted_constraints +' '+ constraint[con]
        if c != len(constraint)-1:
            string += ', '
    return string

def format_constraints(constraints):
    string = 'Constraints:\n\n'
    for e in constraints:
        if constraints[e] == {}:
            formatted_constraint = 'all nr'
        else:
            formatted_constraint = format_constraint(constraints[e])
        string += 'H_'+str(e) +' : ' + formatted_constraint +'\n'
    string +='\n'
    return string
    
def new_file(filename,extension,path):
	
    if path.endswith('/') != True:
        if (os.path.exists(path)!=True):
            os.mkdir(path)
        path += '/'
    file = open(path+filename+'.'+extension, "w")
    return file

def TeX_write(filename,path,TeX_str):
    TeX_file = new_file(filename,'tex',path)
    TeX_file.write(TeX_str)
    TeX_file.close()

def mctdh_write(filename,path,mctdh_str,mctdh_str2):
    mctdh_file = new_file(filename,'op',path)
    mctdh_file.write(mctdh_str)
    mctdh_file.close()
    mctdh_file = new_file(filename,'inp',path)
    mctdh_file.write(mctdh_str2)
    mctdh_file.close()

def log_write(filename,path):
    #write log file to outputs
    global log_str
    log_file = new_file(filename,'log',path)
    log_file.write(log_str)
    log_file.close()

def convert_syntax(sympy_expr,syntax):
    if syntax.upper() not in ['LATEX','MATHEMATICA']:
        raise Exception('ConvertSyntaxError: could not recognize syntax.')
    else:
        if syntax.upper() == 'LATEX':
            return convert_to_LaTeX(sympy_expr)
        elif syntax.upper() == 'MATHEMATICA':
            return convert_to_Mathematica(sympy_expr)

def convert_to_MCTDH(sympy_expr,e,element_position,vhegenclass,e_coords):
    #print()
    #print(f'creating mctdh operator for matrix element {e}')
    #print(element_position[e])
    sympy_str=str(sympy_expr).replace('.0','').split(' ')
    if (e_coords=='cart'):
        if hasattr(vhegenclass, 'mctdhout')==False:
            vhegenclass.mctdhout=mctdh(vhegenclass.modes,e_coords,vhegenclass.filename)
        if (element_position[e][0]<=element_position[e][1]):
            vhegenclass.mctdhout.appendmctdh(sympy_str,element_position[e],e_coords)
    else:
        if hasattr(vhegenclass, 'mctdhoutp')==False:
            vhegenclass.mctdhoutp=mctdh(vhegenclass.modes,e_coords,vhegenclass.filename)
        if (element_position[e][0]<=element_position[e][1]):
            vhegenclass.mctdhoutp.appendmctdh(sympy_str,element_position[e],e_coords)

class MyLatexPrinter(LatexPrinter):
    """Print exp as \exp instead of e^.
    """
    
    
    _default_settings = {
        "full_prec": False,
        "fold_frac_powers": False,
        "fold_func_brackets": False,
        "fold_short_frac": True,
        "inv_trig_style": "abbreviated",
        "itex": False,
        "ln_notation": False,
        "long_frac_ratio": None,
        "mat_delim": "[",
        "mat_str": None,
        "mode": "plain",
        "mul_symbol": None,
        "order": None,
        "symbol_names": {},
        "root_notation": True,
        "mat_symbol_style": "plain",
        "imaginary_unit": "i",
        "gothic_re_im": False,
        "decimal_separator": "period",
        "perm_cyclic": True,
        "parenthesize_super": True,
        "min": None,
        "max": None,
    }  # type: Dict[str, Any]

    def _print_ExpBase(self, expr, exp=None):
        # TODO should exp_polar be printed differently?
        #      what about exp_polar(0), exp_polar(1)?
        tex = r"\exp \left(%s\right)" % self._print(expr.args[0])
        return self._do_exponent(tex, exp)

    def _print_Function(self, expr, exp=None):
        r'''
        Render functions to LaTeX, handling functions that LaTeX knows about
        e.g., sin, cos, ... by using the proper LaTeX command (\sin, \cos, ...).
        For single-letter function names, render them as regular LaTeX math
        symbols. For multi-letter function names that LaTeX does not know
        about, (e.g., Li, sech) use \operatorname{} so that the function name
        is rendered in Roman font and LaTeX handles spacing properly.

        expr is the expression involving the function
        exp is an exponent
        '''
        func = expr.func.__name__
        if hasattr(self, '_print_' + func) and \
                not isinstance(expr, AppliedUndef):
            return getattr(self, '_print_' + func)(expr, exp)
        else:
            args = [str(self._print(arg)) for arg in expr.args]
            # How inverse trig functions should be displayed, formats are:
            # abbreviated: asin, full: arcsin, power: sin^-1
            inv_trig_style = self._settings['inv_trig_style']
            # If we are dealing with a power-style inverse trig function
            inv_trig_power_case = False
            # If it is applicable to fold the argument brackets
            can_fold_brackets = self._settings['fold_func_brackets'] and \
                len(args) == 1 and \
                not self._needs_function_brackets(expr.args[0])

            inv_trig_table = [
                "asin", "acos", "atan",
                "acsc", "asec", "acot",
                "asinh", "acosh", "atanh",
                "acsch", "asech", "acoth",
            ]

            # If the function is an inverse trig function, handle the style
            if func in inv_trig_table:
                if inv_trig_style == "abbreviated":
                    pass
                elif inv_trig_style == "full":
                    func = "arc" + func[1:]
                elif inv_trig_style == "power":
                    func = func[1:]
                    inv_trig_power_case = True

                    # Can never fold brackets if we're raised to a power
                    if exp is not None:
                        can_fold_brackets = False

            if inv_trig_power_case:
                if func in accepted_latex_functions:
                    name = r"\%s^{-1}" % func
                else:
                    name = r"\operatorname{%s}^{-1}" % func
            elif exp is not None:
                func_tex = self._hprint_Function(func)
                func_tex = self.parenthesize_super(func_tex)
                name = r'%s^{%s}' % (func_tex, exp)
            else:
                name = self._hprint_Function(func)

            if can_fold_brackets:
                if func in accepted_latex_functions:
                    # Wrap argument safely to avoid parse-time conflicts
                    # with the function name itself
                    name += r" {%s}"
                else:
                    name += r"%s"
            else:
                name += r" \left(%s\right)"

            if inv_trig_power_case and exp is not None:
                name += r"^{%s}" % exp

            return name % ",".join(args)

def print_my_latex(expr):
    """ Most of the printers define their own wrappers for print().
    These wrappers usually take printer settings. Our printer does not have
    any settings.
    """
    print(MyLatexPrinter().doprint(expr))


def convert_to_LaTeX(sympy_expr):
    #print(sympy_expr)
    latex_expr = MyLatexPrinter().doprint(sympy_expr)#,fold_short_frac=True)
    #latex_expr =  sym.printing.latex(sympy_expr,fold_short_frac=True)
    latex_expr = latex_expr.replace('1.0','')
    latex_expr = latex_expr.replace('.0','')
    latex_expr = latex_expr.replace(R'\left(',R'(')
    latex_expr = latex_expr.replace(R'\right)',R')')
    latex_expr = add_dashes_to_latex(latex_expr)
    return latex_expr

def add_dashes_to_latex(latex_expr):
    tcount=0
    total_length=len(latex_expr)
    while tcount+200<total_length:
        tcount=tcount+200
        #go backwards to find split
        
        for i in range(200):
            tcount=tcount-1
            #print(tcount)
            #print(total_length)
            if (latex_expr[tcount]=='+' or latex_expr[tcount]=='-'):
                if (latex_expr[tcount-1]!=',' and latex_expr[tcount-1]!=R'{'):
                    #print(len(latex_expr))
                    latex_expr=latex_expr[:tcount+1]+R'\\'+'\n'+'&'+latex_expr[tcount+1:]
                    #print(latex_expr[tcount+1:tcount+8])
                    tcount=tcount+4
                    #print(len(latex_expr))
                    total_length+=4
                    break
    return latex_expr


def exec_pdflatex(filename,path):
    if not path.endswith('/'):
        path += '/'
    print('Sending output to '+path)
    if os.path.isfile(path+filename + '.pdf') == True:
        print('Overwriting '+ path+filename + '.pdf')
        os.remove(path+filename + '.pdf')
    print('Generating output PDF:\n')
    try:
        proc = subprocess.Popen(['pdflatex --interaction=batchmode -output-directory='+path, path+filename+'.tex'])
        proc.communicate()
        proc.communicate()
    except OSError:
        try:
            os.system('pdflatex --interaction=batchmode -output-directory='+path+' '+path+filename +'.tex')
            os.system('pdflatex --interaction=batchmode -output-directory='+path+' '+path+filename +'.tex') #double compile for ToC
        except OSError:
            print('Error making '+filename+'.pdf.\n')   
    '''Remove unwanted LaTeX output files .aux and .log'''
    if os.path.isfile(path+filename + '.aux') == True:
        os.remove(path+filename + '.aux')
    if os.path.isfile(path+filename + '.log') == True:
        os.remove(path+filename + '.log')

class mctdh:
    def __init__(self,modes,e_coords,filename):
        modes0=[mode[0] for mode in modes]
        nume=0
        numb=0
        numa=0
        tot=0
        for i in range(len(modes)):
            if modes0[i]=='E':
                nume=nume+1
            if modes0[i]=='A':
                numa=numa+1
            if modes0[i]=='B':
                numb=numb+1
        self.modeargs=[]
        self.massargs=[]
        self.masspos=[]
        for i in range(numa):
            self.modeargs.append(f'z{i+1}')
            self.massargs.append(f'mass_z{i+1}')
            self.masspos.append((tot))
            tot=tot+1
        for i in range(numb):
            self.modeargs.append(f'w{i+1}')
            self.massargs.append(f'mass_w{i+1}')
            self.masspos.append((tot))
            tot=tot+1
        for i in range(nume):
            self.massargs.append(f'mass_e{i+1}')
            if (e_coords=='cart'):
                self.modeargs.append(f'x{i+1}')
                self.modeargs.append(f'y{i+1}')
            else:
                self.modeargs.append(f'rho{i+1}')
                self.modeargs.append(f'phi{i+1}')
            self.masspos.append((tot,tot+1))
            tot=tot+2
        
        self.modeargs.append('el')
        self.nmodes=len(self.modeargs)
        self.dashes=''
        
        self.printmodes=(R' modes').ljust(25)+R' |'
        for j in range(27):
            self.dashes=self.dashes+R'-'
        for i in range(self.nmodes-1):
            if ('phi' in self.modeargs[i]):
                self.printmodes=self.printmodes+' '+self.modeargs[i].ljust(11)+R' |'
                for j in range(14):
                    self.dashes=self.dashes+R'-'
            else:
                self.printmodes=self.printmodes+' '+self.modeargs[i].ljust(4)+R' |'
                for j in range(7):
                    self.dashes=self.dashes+R'-'

        self.printmodes=self.printmodes+' '+self.modeargs[self.nmodes-1]
        for j in range(5):
            self.dashes=self.dashes+R'-'
        if (e_coords=='cart'):
            self.inpfile=(R'RUN-SECTION'+'\n'
                            R'name='+filename+'-cart'+'\n'
                            R'title='+filename+'-cart'+'\n'
                            R''+'\n'
                            R'OPERATOR-SECTION'+'\n'
                            R'opname='+filename+'-cart'+'\n'
                            R'end-operator-section'+'\n'
                            R''+'\n')
            self.title_section=(R'OP_DEFINE-SECTION'+'\n'
                            R'title'+'\n'
                            R'   '+filename+'-cart'+'\n'
                            R'end-title'+'\n'
                            R'end-op_define-section'+'\n'
                            R''+'\n')
        else:
            self.inpfile=(R'RUN-SECTION'+'\n'
                            R'name='+filename+'-pol'+'\n'
                            R'title='+filename+'-pol'+'\n'
                            R''+'\n'
                            R'OPERATOR-SECTION'+'\n'
                            R'opname='+filename+'-pol'+'\n'
                            R'end-operator-section'+'\n'
                            R''+'\n')
            self.title_section=(R'OP_DEFINE-SECTION'+'\n'
                            R'title'+'\n'
                            R'   '+filename+'-pol'+'\n'
                            R'end-title'+'\n'
                            R'end-op_define-section'+'\n'
                            R''+'\n')
        self.addpbasis()
        self.parameter_section=(R'PARAMETER-SECTION'+'\n')
        self.parameters_list=[]
        for i in range(len(self.massargs)):
            self.parameters_list.append(self.massargs[i])
        self.hamiltonian_section=(self.dashes+'\n'
                            R'HAMILTONIAN-SECTION'+'\n'
                            +self.dashes+'\n'
                            +self.printmodes+'\n'
                            +self.dashes+'\n')
        self.add_keo(e_coords)
        
        #self.modecombolist=[]
    
    def add_keo(self,e_coords):
        colsize=[]
        for j in range(self.nmodes):
            if ('phi' in self.modeargs[j]):
                colsize.append(12)
            else:
                colsize.append(5)
        for i in range(len(self.massargs)):
            if ('e' in self.massargs[i]): #need two or three lines for this mass
                if (e_coords=='cart'):
                    for j in range(2): #need two lines for each e-mode
                        self.hamiltonian_section+=(R'-1/2*'+self.massargs[i]).ljust(25)+R' |'
                        for k in range(self.nmodes-1): #check all modes for correct position
                            if (k!=self.masspos[i][j]):
                                self.hamiltonian_section+=(' 1').ljust(colsize[k])+R' |'
                            else:
                                self.hamiltonian_section+=' dq^2'.ljust(colsize[k])+R' |'
                        self.hamiltonian_section+=' 1'.ljust(colsize[k])+'\n'
                else:
                    for j in range(3):
                        if (j!=1):
                            self.hamiltonian_section+=(R'-1/2*'+self.massargs[i]).ljust(25)+R' |'
                        else:
                            self.hamiltonian_section+=(R'-1/8*'+self.massargs[i]).ljust(25)+R' |'
                        l=0
                        for k in range(len(self.massargs)): #check mass args
                            #print([l,self.masspos[i]])
                            if (l!=self.masspos[i][0]):
                                if 'e' in self.massargs[k]:
                                    self.hamiltonian_section+=(' 1').ljust(colsize[l])+R' |'
                                    l=l+1
                                    self.hamiltonian_section+=(' 1').ljust(colsize[l])+R' |'
                                    l=l+1
                                else:
                                    self.hamiltonian_section+=(' 1').ljust(colsize[l])+R' |'
                                    l=l+1
                            else:
                                if (j==0):
                                    self.hamiltonian_section+=' dq^2'.ljust(colsize[l])+R' |'
                                    l=l+1
                                    self.hamiltonian_section+=' 1'.ljust(colsize[l])+R' |'
                                    l=l+1
                                elif (j==1):
                                    self.hamiltonian_section+=' q^-2'.ljust(colsize[l])+R' |'
                                    l=l+1
                                    self.hamiltonian_section+=' 1'.ljust(colsize[l])+R' |'
                                    l=l+1
                                else:
                                    self.hamiltonian_section+=' q^-2'.ljust(colsize[l])+R' |'
                                    l=l+1
                                    self.hamiltonian_section+=' dq^2'.ljust(colsize[l])+R' |'
                                    l=l+1
                                
                        self.hamiltonian_section+=' 1'.ljust(colsize[k])+'\n'
            else:
                self.hamiltonian_section+=(R'1').ljust(25)+' |'
                for k in range(self.nmodes-1):
                    if (k!=self.masspos[i]):
                        self.hamiltonian_section+=(' 1').ljust(colsize[k])+R' |'
                    else:
                        self.hamiltonian_section+=(' KE').ljust(colsize[k])+R' |'
                self.hamiltonian_section+=(' 1').ljust(colsize[k])+'\n'
    
    def addpbasis(self):
        self.inpfile=self.inpfile+'PBASIS-SECTION'+'\n'
        for i in range(self.nmodes):
            if ('rho' in self.modeargs[i]):
                self.inpfile+=self.modeargs[i].ljust(8)+'lagu1'.ljust(8)+'35'.ljust(8)+' 0.0 0.5'+'\n'
            elif ('phi' in self.modeargs[i]):
                self.inpfile+=self.modeargs[i].ljust(8)+'exp'.ljust(8)+'35'.ljust(8)+' 0.00 6.283185307179586 periodic'+'\n'
            elif ('el' in self.modeargs[i]):
                self.inpfile+='el'.ljust(8)+'el'.ljust(8)+'4'.ljust(8)+'\n'
            else:
                self.inpfile+=self.modeargs[i].ljust(8)+'HO'.ljust(8)+'21'.ljust(8)+' 0.0 1.0 1.0'+'\n'
        self.inpfile+='end-pbasis-section'+'\n'


    def appendmctdh(self,sympy_str,pos,e_coords):
        #j=0 is treated differently
        js=0
        
        if (sympy_str[0]!='0'):
            vtt=sympy_str[0].split('*')
            #print('in append')
            #print(sympy_str)
            
            #print(vtt)
            #compose part before coefficient
        
            if (vtt[0][0]=='-'):
                firstcoef='-'
                vtt[0]=vtt[0][1:]
            else:
                firstcoef=None
        
        
            #print('first term')
            #print(vtt[0])
            js=self.startline(vtt,firstcoef)
            #print(vtt[js:])
            if (e_coords=='cart'):
                self.add_hamiltonian_line(vtt[js:],pos)
            else:
                self.add_hamiltonian_line_p(vtt[js:],pos)
            #print(self.hamiltonian_section)
            #exit()
        

        for j in range(1,len(sympy_str)):
            if (sympy_str[j]=='-'):
                firstcoef='-'
            elif (sympy_str[j]=='+'):
                firstcoef=None
            else:
                vtt=sympy_str[j].split('*')
                js=self.startline(vtt,firstcoef)
                #print(vtt[js:])
                if (e_coords=='cart'):
                    self.add_hamiltonian_line(vtt[js:],pos)
                else:
                    self.add_hamiltonian_line_p(vtt[js:],pos)
            
            
        #print(self.hamiltonian_section)

    def startline(self,vtt,firstcoef):
        js=0
        if (vtt[0][0].isdigit()):
            if firstcoef is None:
                if (vtt[0]!='1'):
                    firstcoef=vtt[0]+'*'
            else:
                if (vtt[0]!='1'):
                    firstcoef=firstcoef+vtt[0]+'*'
            js=js+1
        else:
            firstcoef==None
        if (vtt[js]=='sqrt(2)'):
            if (firstcoef is not None):
                firstcoef=firstcoef+'sq2'+'*'
            else:
                firstcoef='sq2'+'*'
            if ('sq2' not in self.parameters_list):
                self.parameters_list.append('sq2')
            js=js+1
        newparam=self.convert_param(vtt[js])
        if (newparam not in self.parameters_list):
            self.parameters_list.append(newparam)
        if (firstcoef is not None):
            firstcoef=firstcoef+newparam
        else:
            firstcoef=newparam
        
        #if vtt[0][0].isdigit():
        #    firstc=vtt[0]+'*'+param
        #elif vtt[0][0]=='-':
        #    if len(vtt[0])>1:
        #        firstc=vtt[0]+'*'+param
        #    else:
        #        firstc=vtt[0]+param
        #else:
        #    firstc=param
        self.hamiltonian_section+=firstcoef.ljust(25)+R' |'
        return js+1

    def add_hamiltonian_line(self,vt,pos):
        lenvt1=len(vt)
        for i in range(lenvt1):
            vt[i]=vt[i].replace('_','')
        #print('----------------')
        #print(self.nmodes)
        for i in range(self.nmodes-1):
            #print(i,self.modeargs[i],vt)
            if self.modeargs[i] in vt:
                mind=vt.index(self.modeargs[i])
                if mind<lenvt1-2:
                    if vt[mind+1]=='':
                        expo=int(vt[mind+2])
                    else:
                        expo=1
                else:
                    expo=1
                if expo>1:
                    self.hamiltonian_section+=(' q'+R'^'+f'{expo}').ljust(5)
                else:
                    self.hamiltonian_section+=(' q').ljust(5)
            else:
                self.hamiltonian_section+=(' 1').ljust(5)
            if i<self.nmodes-1:
                self.hamiltonian_section+=R' |'
        
        self.hamiltonian_section+=f' S{pos[0]}'
        self.hamiltonian_section+=R'&'
        self.hamiltonian_section+=f'{pos[1]}'
        self.hamiltonian_section+='\n'                
                #for i in range(vs,len(vt)):
                #    if vt[i][0].isalpha():
                #        mctdh=mctdh+'*'+vt[i]
                #    if vt[i][0]
    def convert_param(self,param):
        sparam=list(str(param))
        paramout=sparam[0].capitalize()
        for i in range(1,len(sparam)):
            if sparam[i].isalpha() or sparam[i].isdigit() or sparam[i]=='-':
                if sparam[i]=='-':
                    paramout=paramout+'m'

                else:
                    paramout=paramout+sparam[i]
        return paramout

    def finalize(self):
        self.hamiltonian_section+='end-hamiltonian-section'+'\n'
        self.hamiltonian_section+=''+'\n'
        for i in range(len(self.parameters_list)):
            if (self.parameters_list[i]=='sq2'):
                self.parameter_section+=(self.parameters_list[i]+'='+'1.414213562373095'+'\n')
            else:
                self.parameter_section+=(self.parameters_list[i]+'='+'\n')
        self.parameter_section+='end-parameter-section'+'\n'
        self.parameter_section+=''+'\n'
        self.fileform=self.title_section+self.parameter_section+self.hamiltonian_section+'end-operator'+'\n'

    def write_to_screen(self):
        output=self.title_section+self.parameter_section+self.hamiltonian_section+'end-operator'+'\n'
        print(output)

    def add_hamiltonian_line_p(self,vt,pos):
        lenvt1=len(vt)
        for i in range(lenvt1):
            vt[i]=vt[i].replace('_','')
        #print('----------------')
        #print(self.nmodes)
        #print(vt)
        #exit()
        for i in range(self.nmodes-1):
            #print(i,self.modeargs[i],vt)
            phitrue=False
            if ('phi' in self.modeargs[i]):
                phitrue=True
                colsize=12
                sinphi=R'sin('+self.modeargs[i]+R')'
                cosphi=R'cos('+self.modeargs[i]+R')'
            else:
                colsize=5
            if phitrue==False:
                if self.modeargs[i] in vt: 
                    mind=vt.index(self.modeargs[i])
                    if mind<lenvt1-2:
                        if vt[mind+1]=='':
                            expo=int(vt[mind+2])
                        else:
                            expo=1
                    else:
                        expo=1
                    if expo>1:
                        self.hamiltonian_section+=(' q'+R'^'+f'{expo}').ljust(colsize)
                    else:
                        self.hamiltonian_section+=(' q').ljust(colsize)
                else:
                    self.hamiltonian_section+=(' 1').ljust(colsize)
            if phitrue==True:
                #both sin(phin) and cos(phin) can be present
                if sinphi in vt or cosphi in vt:
                    hamadd=None
                    if (sinphi in vt):
                        mind=vt.index(sinphi)
                        if mind<lenvt1-2:
                            if vt[mind+1]=='':
                                expo=int(vt[mind+2])
                            else:
                                expo=1
                        else:
                            expo=1
                        if expo>1:
                            hamadd=(' sin'+R'^'+f'{expo}')
                        else:
                            hamadd=(' sin')
                    if (cosphi in vt):
                        mind=vt.index(cosphi)
                        if mind<lenvt1-2:
                            if vt[mind+1]=='':
                                expo=int(vt[mind+2])
                            else:
                                expo=1
                        else:
                            expo=1
                        if expo>1:
                            if hamadd is not None:
                                hamadd+='*'+('cos'+R'^'+f'{expo}')
                            else:
                                hamadd=(' cos'+R'^'+f'{expo}')
                        else:
                            if hamadd is not None:
                                hamadd+='*'+('cos')
                            else:
                                hamadd=(' cos')
                    self.hamiltonian_section+=hamadd.ljust(colsize)
                else:
                    self.hamiltonian_section+=(' 1').ljust(colsize)
            
            if i<self.nmodes-1:
                self.hamiltonian_section+=R' |'
        
        self.hamiltonian_section+=f' S{pos[0]}'
        self.hamiltonian_section+=R'&'
        self.hamiltonian_section+=f'{pos[1]}'
        self.hamiltonian_section+='\n'                
                #for i in range(vs,len(vt)):
                #    if vt[i][0].isalpha():
                #        mctdh=mctdh+'*'+vt[i]
                #    if vt[i][0]
#EOF
