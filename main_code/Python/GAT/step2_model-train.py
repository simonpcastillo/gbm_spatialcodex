import datetime
import keras.optimizers.legacy
import os

HIDDEN_UNITS = 120
NUM_HEADS = 8
NUM_LAYERS = 5
OUTPUT_DIM = len(class_values)

NUM_EPOCHS = 250
BATCH_SIZE = 1024
VALIDATION_SPLIT = 0.30
LEARNING_RATE = 3e-3
MOMENTUM = 0.9

log_id = datetime.datetime.now().strftime("%Y%m%d-%H%M")
log_dir = "logs/fit/" + log_id
ckpt_path = "ckpt/" + log_id + "/"

my_callbacks = [
    keras.callbacks.EarlyStopping(monitor="val_loss", min_delta=1e-5, patience=10, restore_best_weights=False),
    keras.callbacks.ModelCheckpoint(monitor="val_loss", mode= 'min', filepath=ckpt_path, save_weights_only=True, save_best_only=True,verbose=1),
    keras.callbacks.TensorBoard(log_dir=log_dir, histogram_freq=1),
]
loss_fn = keras.losses.SparseCategoricalCrossentropy(from_logits=True)
optimizer = keras.optimizers.legacy.SGD(LEARNING_RATE, momentum=MOMENTUM)
optimizer = keras.optimizers.legacy.Adam(LEARNING_RATE)
accuracy_fn = keras.metrics.SparseCategoricalAccuracy(name="acc")

# Build model
gat_model = GraphAttentionNetwork(
    node_states, edges, HIDDEN_UNITS, NUM_HEADS, NUM_LAYERS, OUTPUT_DIM
)

# Compile model
gat_model.compile(loss=loss_fn, optimizer=optimizer, metrics=[accuracy_fn])


gat_model.fit(
    x=train_indices,
    y=train_labels,
    validation_split=VALIDATION_SPLIT,
    batch_size=BATCH_SIZE,
    epochs=NUM_EPOCHS,
    callbacks=my_callbacks,
    verbose=2,
)


