from sympy import Matrix, Symbol, conjugate, sqrt

requirements_dct = {0: ['A','A'],
                    1: ['A','B'],
                    2: ['E','A'],
                    3: ['E','B'],
                    4: ['E','E'],
                    5: ['B','B'],
                    6: ['A'],
                    7: ['B'],
                    8: ['E']}

#Structure of matrix_dct values: [original H ,symmetrized H, transformation U]
matrix_dct = {0: [Matrix([[0,'A_alpha__A_beta'],
                          ['A_beta__A_alpha',0]]),

                  Matrix([[0,'A_alpha__A_beta'],
                          ['A_alpha__A_beta',0]]),

                  Matrix([[1,0],
                          [0,1]])],

              1: [Matrix([[0,'A__B'],
                          ['B__A',0]]),

                  Matrix([[0,'A__B'],
                          ['A__B',0]]),

                  Matrix([[1,0],
                          [0,1]])],

              2: [Matrix([[0,0,Symbol('+__A')],
                          [0,0,Symbol('-__A')],
                          [Symbol('A__+'),Symbol('A__-'),0]]),

                  Matrix([[0,0,Symbol('+__A')],
                          [0,0,conjugate(Symbol('+__A'))],
                          [conjugate(Symbol('+__A')),Symbol('+__A'),0]]),

                  Matrix([[1,1,0],
                          [1j,-1j,0],
                          [0,0,sqrt(2)]])*1/sqrt(2)],

              3: [Matrix([[0,0,Symbol('+__B')],
                          [0,0,Symbol('-__B')],
                          [Symbol('B__+'),Symbol('B__-'),0]]),

                  Matrix([[0,0,Symbol('+__B')],
                          [0,0,conjugate(Symbol('+__B'))],
                          [conjugate(Symbol('+__B')),Symbol('+__B'),0]]),

                  Matrix([[1,1,0],
                          [1j,-1j,0],
                          [0,0,sqrt(2)]])*1/sqrt(2)],

              4: [Matrix([[0,0,Symbol('+_alpha__+_beta'),Symbol('+_alpha__-_beta')],
                          [0,0,Symbol('-_alpha__+_beta'),Symbol('-_alpha__-_beta')],
                          [Symbol('+_beta__+_alpha'),Symbol('+_beta__-_alpha'),0,0],
                          [Symbol('-_beta__+_alpha'),Symbol('-_beta__-_alpha'),0,0]]),

                  Matrix([[0,0,Symbol('+_alpha__+_beta'),Symbol('+_alpha__-_beta')],
                          [0,0,conjugate(Symbol('+_alpha__-_beta')),conjugate(Symbol('+_alpha__+_beta'))],
                          [conjugate(Symbol('+_alpha__+_beta')),Symbol('+_alpha__-_beta'),0,0],
                          [conjugate(Symbol('+_alpha__-_beta')),Symbol('+_alpha__+_beta'),0,0]]),

                  Matrix([[1,1,0,0],
                          [1j,-1j,0,0],
                          [0,0,1,1],
                          [0,0,1j,-1j]])*1/sqrt(2)],

              5: [Matrix([[0,'B_alpha__B_beta'],
                          ['B_beta__B_alpha',0]]),

                  Matrix([[0,'B_alpha__B_beta'],
                          ['B_alpha__B_beta',0]]),

                  Matrix([[1,0],
                          [0,1]])],

              6: [Matrix([['A__A']]),

                  Matrix([['A__A']]),

                  Matrix([[1]])],

              7: [Matrix([['B__B']]),

                  Matrix([['B__B']]),

                  Matrix([[1]])],

              8: [Matrix([[Symbol('+__+'), Symbol('+__-')],
                          [Symbol('-__+'), Symbol('-__-')]]),

                  Matrix([[Symbol('+__+'), Symbol('+__-')],
                          [conjugate(Symbol('+__-')), Symbol('+__+')]]),

                  Matrix([[1,1],
                          [1j,-1j]])*1/sqrt(2)]

              }

def return_matrices(symmetry,states):
    #from modules.input import return_acceptable_states
    requirements = [s[0] for s in states]
    print(requirements)
    for k in requirements_dct:
        if requirements_dct[k] == requirements:
            
            return matrix_dct[k][0], matrix_dct[k][1], matrix_dct[k][2]
    #matrix3=matrix_dct[k][2]
    #print(states)
    #if requirements[0]=='E' or requirements[1]=='E':
    #    ssub=[]
    #    acceptable=return_acceptable_states(symmetry)
    #    for k in range(len(acceptable.full)):
    #            if requirements[0]=='E' and requirements[1]!='E':
    #                    if states[0]==acceptable.full[k]:
    #                            ssub.append(acceptable.labels[k][1])
    #                            break
    #            if requirements[0]!='E' and requirements[1]=='E':
    #                    if states[1]==acceptable.full[k]:
    #                            ssub.append(acceptable.labels[k][1])
    #                            break
    #            if requirements[0]=='E' and requirements[1]=='E':
    #                    if states[0]==acceptable.full[k]:
    #                            ssub.append(acceptable.labels[k][1])
    #                    if states[1]==acceptable.full[k]:
    #                            ssub.append(acceptable.labels[k][1])
        
    #    if len(ssub)==1:


    
    print('Error: Could not find matrix form for states '+str(states)+'.')
    exit()

