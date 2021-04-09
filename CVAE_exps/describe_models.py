import glob
import numpy as np
import subprocess
from keras.models import load_model

models = glob.glob("cvae_runs_*/best.h5")

for model in models:
    a = load_model(model)
    a.summary()


