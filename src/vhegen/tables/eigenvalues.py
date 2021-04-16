from vhegen.modules.electronic import kdelta 
from vhegen.modules.input import return_acceptable_states
#, refl_parity, inver_parity
from sympy import sqrt, Symbol
from math import exp,pi

def return_eigenvals(symmetry, states):
	
	#print(states)
	acceptable=return_acceptable_states(symmetry)
	requirements = [s[0] for s in states]
	if len(states)>1:
		if requirements[0]=='E' or requirements[1]=='E':
			ssub=[]
			
			for k in range(len(acceptable.full)):
				if requirements[0]=='E' and requirements[1]!='E':
					if states[0]==acceptable.full[k]:
						ssub.append(acceptable.labels[k][1])
						break
				if requirements[0]!='E' and requirements[1]=='E':
					if states[1]==acceptable.full[k]:
						ssub.append(acceptable.labels[k][1])
						break
				if requirements[0]=='E' and requirements[1]=='E':
					if states[0]==acceptable.full[k]:
						ssub.append(acceptable.labels[k][1])
					if states[1]==acceptable.full[k]:
						ssub.append(acceptable.labels[k][1])
	else:
		if requirements[0]=='E':
			ssub=[]
			for k in range(len(acceptable.full)):
				if states[0]==acceptable.full[k]:
					ssub.append(acceptable.labels[k][1])
		

	n = symmetry['rot']
	labels=[]
	for i in range(len(states)):
		for k in range(len(acceptable.full)):
			if states[i]==acceptable.full[k]:
				labels.append(acceptable.labels[k])

	refl_parity=[]
	inver_parity=[]
	for i in range(len(states)):
		refl_parity.append(labels[i][1])

		if (labels[i][2]!=0):
			inver_parity.append(pow((-1),labels[i][2]+1))
		else:
			inver_parity.append(0)
	#print(refl_parity)
	#exit()
	if (len(states)==1):
		refl_parity.append(0)
		inver_parity.append(0)
		#print(labels)
		#print(refl_parity)
		#print(inver_parity)


	added_dummy = False
	if len(states) == 1:
		added_dummy = True
		states.append('dummystate')

	if n % 2 == 1: #(odd n-axial) problems

		requirements_dct = {0: ['A'],
							1: ['E'],
							2: ['A','A'],
							3: ['E','A'],
							4: ['E','E']}

		eigenvals_dct = {0: {Symbol('A__A'):[0,1,0,1]},

        				 1: {Symbol('+__+'): [0,1,0,1],
                             Symbol('+__-'): [2*refl_parity[0],1,-1,1]},

						 2: {Symbol('A_alpha__A_beta'):[0,(-1)**(kdelta(refl_parity[0],refl_parity[1])+1),0,(-1)**(kdelta(inver_parity[0],inver_parity[1])+1)]},
						   
						 3: {Symbol('+__A'):[refl_parity[0],(-1)**kdelta(refl_parity[1],2),(-1)**kdelta(refl_parity[1],1),(-1)**(kdelta(inver_parity[0],inver_parity[1])+1)]},
						   
						 4: {Symbol('+_alpha__+_beta'):[refl_parity[0]-refl_parity[1],1,-1,(-1)**(kdelta(inver_parity[0],inver_parity[1])+1)],
                             Symbol('+_alpha__-_beta'):[refl_parity[0]+refl_parity[1],1,-1,(-1)**(kdelta(inver_parity[0],inver_parity[1])+1)]}
        				   
                        }

	else: #(even n-axial) problems

		requirements_dct = {0: ['A'],
							1: ['B'],
							2: ['E'],
							3: ['A','A'],
							4: ['A','B'],
							5: ['E','A'],
							6: ['E','B'],
							7: ['E','E'],
							8: ['B','B']}

		eigenvals_dct = {0: {Symbol('A__A'):[0,1,0,1]},

						 1: {Symbol('B__B'):[0,1,0,1]},

						 2: {Symbol('+__+'): [0,1,0,1],
                             Symbol('+__-'): [refl_parity[0]*2,1,-1,1]},

						 3: {Symbol('A_alpha__A_beta'):[0,(-1)**(kdelta(refl_parity[0],refl_parity[1])+1),0,(-1)**(kdelta(inver_parity[0],inver_parity[1])+1)]},
							
						 4: {Symbol('A__B'):[n//2,(-1)**(kdelta(refl_parity[0],refl_parity[1])+1),0,(-1)**(kdelta(inver_parity[0],inver_parity[1])+1)]},

						 5: {Symbol('+__A'):[refl_parity[0],(-1)**kdelta(refl_parity[1],2),(-1)**kdelta(refl_parity[1],1),(-1)**(kdelta(inver_parity[0],inver_parity[1])+1)]},
						
						 6: {Symbol('+__B'):[refl_parity[0]+n//2,(-1)**kdelta(refl_parity[1],2),(-1)**kdelta(refl_parity[1],1),(-1)**(kdelta(inver_parity[0],inver_parity[1])+1)]},

						 7: {Symbol('+_alpha__+_beta'):[refl_parity[0]-refl_parity[1],1,-1,(-1)**(kdelta(inver_parity[0],inver_parity[1])+1)],
                             Symbol('+_alpha__-_beta'):[refl_parity[0]+refl_parity[1],1,-1,(-1)**(kdelta(inver_parity[0],inver_parity[1])+1)]}, 

						 8: {Symbol('B_alpha__B_beta'):[0,(-1)**(kdelta(refl_parity[0],refl_parity[1])+1),0,(-1)**(kdelta(inver_parity[0],inver_parity[1])+1)]}
						}
						
	if added_dummy == True:
		del states[-1]
		
	problem_reqs = [s[0] for s in states]

	for k in requirements_dct:
		if requirements_dct[k] == problem_reqs:
			return eigenvals_dct[k]

	print('Error: Could not find eigenvalues for states '+str(states)+'.')
	exit()

#EOF
