import numpy as np

P = np.array([[0.4, 0.3, 0.2],
              [0.5, 0.5, 0.2],
              [0.1, 0.2, 0.6]])

v = np.array([20000, 25000, 8000])

for i in range(20):
    print i
    v = np.dot(P, v)

    if i == 3 or i == 9 or i == 19:
        print v

#Dette printer:
"""
[ 16078.1  22179.7  14742.2]
[ 16000.4953553  22001.1406921  14998.3639526]
[ 16000.00010775  22000.00024812  14999.99964414]
"""

#Kan skrives uten bruk av np.dot()
"""
for i in range(20):
    new_q = np.zeros(3)
    for j in range(3):
        for k in range(3):
            new_q[j] += q[k] * P[j][k]
    q = new_q

print q
"""