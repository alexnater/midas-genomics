import argparse,time,sys
import numpy as np
import pandas as pd
from keras.models import model_from_json
from keras.models import Sequential, Model
from keras import optimizers
from keras.layers import Dense, Dropout, Activation, Flatten, Input
from keras.layers import Conv2D, MaxPooling2D, concatenate
from keras.utils import np_utils
from sklearn.model_selection import train_test_split
from keras.utils.layer_utils import convert_all_kernels_in_model
from keras.preprocessing.image import ImageDataGenerator
from keras.callbacks import EarlyStopping, ModelCheckpoint
import keras.backend as K
import fnmatch


def main(argsDict):
    numSubWins = argsDict['numSubWins']
    print("loading data now...")
    
    # import data from predictFile:
    x_df = pd.read_table(argsDict['predictFile'])
    testX = x_df[list(x_df)[4:]].as_matrix()
    nDims = int(testX.shape[1] / numSubWins)
    testX = testX.reshape(testX.shape[0], nDims, numSubWins, 1)
    
    # set up generator for normalization 
    validation_gen = ImageDataGenerator(
        featurewise_center=True,
        featurewise_std_normalization=True,
        horizontal_flip=False)
    validation_gen.fit(testX)
    
    # import model
    json_file = open(argsDict['modelStructure'], 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    model = model_from_json(loaded_model_json)
    
    # load weights into new model
    model.load_weights(argsDict['modelWeights'])
    print("Loaded model from disk")
    
    # get predictions
    preds = model.predict(validation_gen.standardize(testX))
    predictions = np.argmax(preds, axis=1)
    classDict = {0:'neutral', 1:'shared', 2:'divergent1', 3:'divergent2'}
    
    # output the predictions
    with open(argsDict['predictFileOutput'], 'w') as pred_file:
        pred_file.write('chrom\tclassifiedWinStart\tclassifiedWinEnd\tpredClass\tprob(neutral)\tprob(shared)\tprob(divergent1)\tprob(divergent2)\n')
        for index, row in x_df.iterrows():
            pred_file.write('{}\t{}\t{}\t{}\t{:f}\t{:f}\t{:f}\t{:f}\n'.format(int(row['chrom']), int(row['start']), int(row['end']), classDict[predictions[index]], 
                            preds[index][0], preds[index][1], preds[index][2], preds[index][3]))
    print("{} predictions complete".format(index+1))



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('modelStructure', help='path to CNN structure .json file')
    parser.add_argument('modelWeights', help='path to CNN weights .h5 file')
    parser.add_argument('predictFile', help='input file to predict')
    parser.add_argument('predictFileOutput', help='output file name')
    parser.add_argument('--numSubWins', type=int, help='number of subwindows that our chromosome is divided into (default = 11)', default=11)
    main(vars(parser.parse_args()))

