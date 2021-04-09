import pickle
from OutlierDB import *


fn = 'OutlierDB.pickle'
my_lock(fn)

with open(fn, 'rb') as f:
    db = pickle.load(f)

print(dir(db))

md5 = db.sorted_index[0]

for i in range(5):
    pdbfile = db.next()
    print(pdbfile)

with open(fn, 'wb') as f:
    pickle.dump(db, f)


my_unlock(fn)
