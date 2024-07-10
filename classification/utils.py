import numpy as np
import xgboost as xgb
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA

TMP = 10

# Function to split the data into training and validation, given embeddings, labels and a split ratio
def split_train_test(embeddings, embedding_labels, TRAIN_RATIO=0.8):
    #import pdb; pdb.set_trace()
    # The data will be divided into 80% training 20% validation
    idx_train = int(len(embeddings) * TRAIN_RATIO)

    indices = np.arange(embeddings.shape[0])

    # Shuffle data
    np.random.shuffle(indices)

    # Separate the images and the labels
    embeddings_shuffle = embeddings[indices] 
    embedding_labels_shuffle = embedding_labels[indices]

    # Split to train and validation
    train_embeddings = embeddings_shuffle[:idx_train]
    train_labels = embedding_labels_shuffle[:idx_train]
    validation_embeddings = embeddings_shuffle[idx_train:]
    validation_labels = embedding_labels_shuffle[idx_train:]
    return train_embeddings, train_labels, validation_embeddings, validation_labels


# xgBoost classifier
def xgb_classify(train_embeddings, embed_train_labels_num, validation_embeddings, embed_val_labels_num, seed=None):
    # Initialize classifier
    xgb_classifier = xgb.XGBClassifier(n_jobs=-1, random_state=seed)
    # Fit
    xgb_classifier.fit(train_embeddings, embed_train_labels_num)
    # Predict
    preds = xgb_classifier.predict(validation_embeddings)
    # Score
    acc_score = accuracy_score(embed_val_labels_num, preds)
    # convert to float
    acc_score = float(acc_score)
    return acc_score*100, preds, xgb_classifier


# Function to plot the ROC curve
def plot_roc_curve(fpr, tpr, roc_auc, label_type_1, label_type_2, save_path=None):
    plt.rcParams.update({'font.size': 12})
    plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
             lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC ' + label_type_1 + ' ' + label_type_2)
    plt.legend(loc="lower right")
    if save_path is not None:
        plt.savefig(save_path, dpi=1200)
    plt.show()


# LDA classifier
def lda_classify(train_embeddings, embed_train_labels_num, validation_embeddings, embed_val_labels_num, seed=None):
    # Initialize classifier
    lda_classifier = LDA(n_components=1)
    # Fit
    lda_classifier.fit(train_embeddings, embed_train_labels_num)
    # Predict
    preds = lda_classifier.predict(validation_embeddings)
    # Score
    acc_score = accuracy_score(embed_val_labels_num, preds)
    return acc_score*100, preds, lda_classifier