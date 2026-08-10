import tensorflow as tf
def mse(y_true, y_pred):
    return tf.reduce_mean(tf.square(y_true - y_pred))
