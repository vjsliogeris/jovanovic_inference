import pickle
import numpy as np


n_steps = 10000

env = pickle.load(open("PBN_jovanovic.pkl","rb"))

#env.reset(np.array([0,0,0,0,0,0,0,0,0,0,0,0,0]))
env.reset()

print("Starting state: {0}\n".format(env.render().astype(int)))

for i in range(n_steps):
    env.step()
    state = env.render().astype(int)
    print(state)
