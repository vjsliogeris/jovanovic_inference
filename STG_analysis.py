import numpy as np
import pickle
import gym
import networkx as nx
import matplotlib.pyplot as plt

#DO NOT PERFORM THIS WITH LARGE PBNs - WILL BREAK. THIS IS AROUND (i think) O(2^n), where n is the size of the PBN.

environment = pickle.load(open("PBN_jovanovic.pkl", "rb"))
STG = environment.render(mode='STG')

attr_gen = nx.attracting_components(STG)
for attr in attr_gen:
    print(attr)

