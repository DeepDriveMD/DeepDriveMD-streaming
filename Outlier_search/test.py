from OutlierDB import *
import pickle

f = open('OutlierDB.pickle','rb')
db = pickle.load(f)
f.close

db.is_consistent()
