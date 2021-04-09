from multiprocessing import Pool

def f(x):
    fn = x[1]
    f = open(fn,"w")
    f.write(f"{x[0]*x[0]}\n")
    f.close()
    return


a = list(range(30))
b = list(map(lambda x: f"{x}.out", a))

c = zip(a,b)

with Pool(5) as p:
        p.map(f, c)
