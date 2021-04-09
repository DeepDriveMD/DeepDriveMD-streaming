import pickle
from OutlierDB import *


fn = 'OutlierDB.pickle'
my_lock(fn)

with open(fn, 'rb') as f:
    db = pickle.load(f)

my_unlock(fn)

print(len(db.sorted_index))

for k in db.sorted_index:
    print(db.dictionary[k])



