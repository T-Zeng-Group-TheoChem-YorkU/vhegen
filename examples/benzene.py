#VHEGEN testrun script

import vhegen as VHE
from time import time

problem = VHE.modules.input.prepare_input('D6h',"E1g","e2g+e2g+e1u+e1g+e2u",'0,6','benzene')
t1 = time()

vhegen_instance = VHE.VHEGEN(problem)

print(vhegen_instance.return_init())

vhegen_instance.set_e_coordinates('both')

vhegen_instance.set_basis('both')

vhegen_instance.auto()

vhegen_instance.pdflatex()

vhegen_instance.mctdhoutput()

print('\nJob complete after '+str(round(time()-t1,3))+' seconds.')
print('\nTestrun complete without errors.')
