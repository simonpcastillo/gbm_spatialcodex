import tensorflow as tf
from tensorflow import keras
from keras import layers
import numpy as np
import pandas as pd
import os
import warnings
from sklearn.model_selection import train_test_split

warnings.filterwarnings("ignore")
np.random.seed(97)

data_dir = '/Users/spcastillo/Downloads/gbm_fromseadragon_20231109/output'
links = pd.read_csv(
    os.path.join(data_dir, "links_30umrand_all.csv"),
    sep=",",
    header=0,
    names=["target", "source"],
    )

features = pd.read_csv(
    os.path.join(data_dir, "features_all.csv"),
    sep=",",
    header=0,
    )
features = features.drop(columns=['phenotype', 'size','CD68', 'Nestin'])
features = features.rename(columns={"phenotype2": "phenotype", 'cell.id': 'cell_id'})
class_values = sorted(features["phenotype"].unique())
class_idx = {name: id for id, name in enumerate(class_values)}
cell_idx = {name: idx for idx, name in enumerate(sorted(features["cell_id"].unique()))}

features["cell_id"] = features["cell_id"].apply(lambda name: cell_idx[name])
links["source"] = links["source"].apply(lambda name: cell_idx[name])
links["target"] = links["target"].apply(lambda name: cell_idx[name])
features["phenotype"] = features["phenotype"].apply(lambda value: class_idx[value])

train_data, test_data = train_test_split(features, test_size=0.5, stratify=features[['phenotype']])

train_indices = train_data["cell_id"].to_numpy()
test_indices = test_data["cell_id"].to_numpy()

train_labels = train_data["phenotype"].to_numpy()
test_labels = test_data["phenotype"].to_numpy()

# Define graph, namely an edge tensor and a node feature tensor
edges = tf.convert_to_tensor(links[["target", "source"]])
x = np.asarray(features.sort_values("cell_id").iloc[:, 1:-1]).astype('float32')
node_states = tf.convert_to_tensor(x)

