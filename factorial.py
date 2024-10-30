import numpy as np
#commented solutions not important for FEM programm
"""j=0
a = np.random.rand(1000)
if a[0] >= 0.5:
    print('a1>=0.5')
else:
    print('a1<0.5')
n=0
for i in a:
    if i >= 0.5:
        n+=1
print(n)
n=-1
while j in range(len(a)):
    if 0.499 <= a[j] and a[j]<=0.501:
        n = j
        break
    j+=1
if n == -1:
    print('Kein Element 0.499 <= a_i <= 0.501')
else: print('Index {} und Wert {}'.format(n,a[n]))"""

def facultaet(n):
    if(n<=1):
        return 1
    else:
        return n*facultaet(n-1)
