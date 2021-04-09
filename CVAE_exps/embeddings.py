import numpy as np
from sklearn.manifold import TSNE

import sys

ps = int(sys.argv[1])
components = int(sys.argv[2])

cvae_embeddings = np.load("projection.npy")
rmsd = np.load("rmsd.npy")
tsne = TSNE(n_components=components, n_jobs=ps)

tsne_embeddings = tsne.fit_transform(cvae_embeddings)

with open(f'tsne_embeddings_{components}.npy', 'wb') as f:
    np.save(f, tsne_embeddings)

