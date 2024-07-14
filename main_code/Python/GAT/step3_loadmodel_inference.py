import datetime
import keras.optimizers.legacy
from keras import layers
import tensorflow as tf
from tensorflow import keras
import numpy as np
import pandas as pd
import os
import warnings
from sklearn.model_selection import train_test_split
import re

HIDDEN_UNITS = 120
NUM_HEADS = 8
NUM_LAYERS = 5
OUTPUT_DIM = 4 #len(class_values)
LEARNING_RATE = 3e-3

warnings.filterwarnings("ignore")
data_dir = '/Users/spcastillo/Downloads/gbm_fromseadragon_20231109/output_test'

for featuresdf in os.listdir(data_dir):
    if featuresdf.startswith("features"):
        roi_id = re.findall('_(.*).csv', featuresdf)[0]
        print(roi_id)
#       if not os.path.isfile('filename.txt'):
        links = pd.read_csv(
            os.path.join(data_dir, "links_1um_"+roi_id+".csv"),
            sep=",",
            header=0,
            names=["target", "source"],
            )

        features = pd.read_csv(
            os.path.join(data_dir, "features_"+roi_id+".csv"),
            sep=",",
            header=0,
            )


        ### DO NOT CHANGE
        features = features.drop(columns=['phenotype', 'size','CD68', 'Nestin'])
        features = features.rename(columns={"phenotype2": "phenotype", 'cell.id': 'cell_id'})
        #####
        class_values = sorted(features["phenotype"].unique())
        class_idx = {name: id for id, name in enumerate(class_values)}
        cell_idx = {name: idx for idx, name in enumerate(sorted(features["cell_id"].unique()))}

        features["cell_id"] = features["cell_id"].apply(lambda name: cell_idx[name])
        links["source"] = links["source"].apply(lambda name: cell_idx[name])
        links["target"] = links["target"].apply(lambda name: cell_idx[name])
        features["phenotype"] = features["phenotype"].apply(lambda value: class_idx[value])

        test_data = features #train_test_split(features, test_size=0.5, stratify=features[['phenotype']])

        # train_indices = train_data["cell_id"].to_numpy()
        test_indices = test_data["cell_id"].to_numpy()

        # train_labels = train_data["phenotype"].to_numpy()
        test_labels = test_data["phenotype"].to_numpy()

        # Define graph, namely an edge tensor and a node feature tensor
        edges = tf.convert_to_tensor(links[["target", "source"]])
        x = np.asarray(features.sort_values("cell_id").iloc[:, 1:-1]).astype('float32')
        node_states = tf.convert_to_tensor(x)


        # Load weights
        loss_fn = keras.losses.SparseCategoricalCrossentropy(from_logits=True)
        optimizer = keras.optimizers.legacy.Adam(LEARNING_RATE)
        accuracy_fn = keras.metrics.SparseCategoricalAccuracy(name="acc")

        # Build model
        gat_model = GraphAttentionNetwork(
            node_states, edges, HIDDEN_UNITS, NUM_HEADS, NUM_LAYERS, OUTPUT_DIM
        )

        # Compile model
        gat_model.compile(loss=loss_fn, optimizer=optimizer, metrics=[accuracy_fn])

        # Load weights
        gat_model.load_weights('ckpt/20240331-1824/')


        test_probs = gat_model.predict(x=test_indices)
        mapping = {v: k for (k, v) in class_idx.items()}

        for i, (probs, label) in enumerate(zip(test_probs[:1], test_labels[:1])):
            print(f"Example {i+1}: {mapping[label]}")
            for j, c in zip(probs, class_idx.keys()):
                print(f"test: {i+1}, GT: {mapping[label]}, pred: {c: <24}, prob:{j:7.3f}")
                #print(f"\tProbability of {c: <24} = {j*100:7.3f}%")
            print("---" * 20)

            d = []
            for i, (probs, label) in enumerate(zip(test_probs[:], test_labels[:])):
                for j, c in zip(probs, class_idx.keys()):

                    d.append(
                        {
                            'test': i,
                            'gt': mapping[label],
                            'gat': f"{c}",
                            "prob": f"{j: 7.3f}"

                        }
                    )

            d = pd.DataFrame(d)
            d.to_csv('/Users/spcastillo/Downloads/gbm_fromseadragon_20231109/predictions/spatialagnostic_1um/'+ roi_id +'.csv')
