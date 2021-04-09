import pickle

class A:
    def __init__(self):
        self.a = 10
        self.b = [0, 1]
        self.c = {'a':10,'b':11}
        self.d = "test"
    def __repr__(self):
        return f"a={self.a}, b={self.b}, c={self.c}, d={self.d}"


aa = A()

print(aa)

with open('test.pickle', 'wb') as f:
    pickle.dump(aa, f)

