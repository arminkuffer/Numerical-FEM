from factorial import facultaet

def tester():
    if((abs(facultaet(0))-1) > 10**(-12)): 
        print('facultaet function failed for input 0: tolerance: {}'.format((facultaet(0))/1))
    else:
        print('facultaet function passed for input 0')
    if((abs(facultaet(5))-120) > 10**(-12)): 
        print('facultaet function failed for input 5: tolerance: {}'.format((facultaet(5))/120))
    else: print('facultaet modul passed for input 5')
tester()