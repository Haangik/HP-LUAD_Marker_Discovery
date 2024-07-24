#https://github.com/FelSiq/DBCV?tab=readme-ov-file#multiprocessing
import numpy as np
import dbcv

#X = np.loadtxt("module_cluster_embed.txt")
#label = np.loadtxt("module_cluster_label.txt")

#score = dbcv.dbcv(X, label, noise_id=0)
#print(score)

def dbcv_execute(n):
  embedding_filename="module_cluster_embed_" + str(n) + ".txt"
  label_filename="module_cluster_label_" + str(n) + ".txt"
  X = np.loadtxt(embedding_filename)
  label = np.loadtxt(label_filename)
  score = dbcv.dbcv(X, label, noise_id=0)
  return score
