#VHEGEN testrun script

import vhegen as VHE
from time import time

problem = VHE.modules.input.prepare_input('C5h',"E1''","E1'+E1''",'0,3','C5h')
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
