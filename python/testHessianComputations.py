import numpy as np

G1 = np.zeros((3,3))
G1[1,2] = -1
G1[2,1] = 1
G2 = np.zeros((3,3))
G2[0,2] = 1
G2[2,0] = -1
G3 = np.zeros((3,3))
G3[0,1] = -1
G3[1,0] = 1

Gs = [G1, G2, G3]

for l,Gl in enumerate(Gs):
  print Gl + Gl.T

for l,Gl in enumerate(Gs):
  for m,Gm in enumerate(Gs):
    print m,l, "\n", Gl.dot(Gm)

a = np.array([1,2,3])
b = np.array([4,5,6])
H = np.zeros((3,3))
for l,Gl in enumerate(Gs):
  for m,Gm in enumerate(Gs):
    H[l,m] =  0.5*a.dot((Gl.dot(Gm)+Gm.dot(Gl)).dot(b))
print H
raw_input()

A = np.resize(np.arange(9),(3,3))
A = A.T+A
print A
for l,Gl in enumerate(Gs):
  for m,Gm in enumerate(Gs):
    print m,l
    print A.dot(Gl.dot(Gm))
    print Gl.dot(A.dot(Gm))
    print Gm.dot(A.dot(Gl)) - Gl.dot(A.dot(Gm))
    print (A.dot(Gl) - Gl.dot(A))
