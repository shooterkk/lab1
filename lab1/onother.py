import numpy as np
class illia():
    i = 24
    z = 35
    def tostring(self):
        rep = str(self.i)+ str(self.z)
        return rep
a = illia()
f = open('output.txt', 'a+')
f.write(a.tostring())
f.close()


a = np.matrix([[1, 2], [3, 4]])
b = np.matrix([[10], [20]])
dt = np.dot(a, b)
print(np.dot(np.transpose(dt), b))

print(np.linalg.norm(b, np.inf))