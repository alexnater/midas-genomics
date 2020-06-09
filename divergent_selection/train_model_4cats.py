import argparse,time,sys
import numpy as np
from keras.models import Sequential, Model
from keras import optimizers
from keras.layers import Dense, Dropout, Activation, Flatten, Input
from keras.layers import Conv2D, MaxPooling2D,concatenate
from keras.utils import np_utils
from sklearn.model_selection import train_test_split
from keras.utils.layer_utils import convert_all_kernels_in_model
from keras.preprocessing.image import ImageDataGenerator
from keras.callbacks import EarlyStopping,ModelCheckpoint
import keras.backend as K
import fnmatch


def same_length_arrays(array_list, nentries = None, insist=False):
    minnumber = min([len(x) for x in array_list])
    if nentries != None:
        if minnumber >= nentries:
            minnumber = nentries
        else:
            print("WARNING: minimum number of entries for training set ({}) is lower than requested number to keep ({})!".format(minnumber, nentries))
            if insist:
                minnumber = nentries
    print("Number of entries to keep for training sets: " + str(minnumber))
    for i in range(len(array_list)):
        np.random.shuffle(array_list[i])
        array_list[i] = array_list[i][:minnumber]



def main(argsDict):
    numSubWins = argsDict['numSubWins']
    trainingDir = argsDict['trainDir']
    testingDir = argsDict['testDir']
    equal = argsDict['equal']
    insist = argsDict['insist']
    keep_train = argsDict['keep_train']
    keep_test = argsDict['keep_test']
    epochOption = argsDict['epochs']
    patienceOption = argsDict['patience']
    outputModel = argsDict['outputModel']
    print("loading data now ...")
    
    cat0 = np.loadtxt(trainingDir + "neutral.fvec", skiprows=1)
    nDims = int(cat0.shape[1] / numSubWins)
    cat0 = np.reshape(cat0, (cat0.shape[0], nDims, numSubWins))
    cat1 = np.loadtxt(trainingDir + "shared.fvec", skiprows=1)
    cat1 = np.reshape(cat1, (cat1.shape[0], nDims, numSubWins))
    cat2 = np.loadtxt(trainingDir + "divergent1.fvec", skiprows=1)
    cat2 = np.reshape(cat2, (cat2.shape[0], nDims, numSubWins))
    cat3 = np.loadtxt(trainingDir + "divergent2.fvec", skiprows=1)
    cat3 = np.reshape(cat3, (cat3.shape[0], nDims, numSubWins))
    
    trainsets = [cat0, cat1, cat2, cat3]
    if equal: same_length_arrays(trainsets, keep_train, insist)
    both = np.concatenate(trainsets)
    both = both.reshape(both.shape[0], nDims, numSubWins, 1)
    y = np.concatenate([np.repeat(x, len(y)) for x, y in enumerate(trainsets)])
    
    if (trainingDir == testingDir):
        X_train, X_test, y_train, y_test = train_test_split(both, y, test_size=0.2)
    else:
        X_train = both
        y_train = y
        # read testing data:
        cat0 = np.loadtxt(testingDir + "neutral.fvec", skiprows=1)
        nDims = int(cat0.shape[1] / numSubWins)
        cat0 = np.reshape(cat0, (cat0.shape[0], nDims, numSubWins))
        cat1 = np.loadtxt(testingDir + "shared.fvec", skiprows=1)
        cat1 = np.reshape(cat1, (cat1.shape[0], nDims, numSubWins))
        cat2 = np.loadtxt(testingDir + "divergent1.fvec", skiprows=1)
        cat2 = np.reshape(cat2, (cat2.shape[0], nDims, numSubWins))
        cat3 = np.loadtxt(testingDir + "divergent2.fvec", skiprows=1)
        cat3 = np.reshape(cat3, (cat3.shape[0], nDims, numSubWins))

        trainsets = [cat0, cat1, cat2, cat3]
        if equal: same_length_arrays(trainsets, keep_test, insist)
        both2 = np.concatenate(trainsets)
        X_test = both2.reshape(both2.shape[0], nDims, numSubWins, 1)
        y_test = np.concatenate([np.repeat(x, len(y)) for x, y in enumerate(trainsets)])
    
    Y_train = np_utils.to_categorical(y_train, 4)
    Y_test = np_utils.to_categorical(y_test, 4)
    X_valid, X_test, Y_valid, Y_test = train_test_split(X_test, Y_test, test_size=0.5)
    
    datagen = ImageDataGenerator(
        featurewise_center=True,
        featurewise_std_normalization=True,
        horizontal_flip=True)
    
    validation_gen = ImageDataGenerator(
        featurewise_center=True,
        featurewise_std_normalization=True,
        horizontal_flip=False)
    
    test_gen = ImageDataGenerator(
        featurewise_center=True,
        featurewise_std_normalization=True,
        horizontal_flip=False)
    
    print("training set has %d examples" % X_train.shape[0])
    print("validation set has %d examples" % X_valid.shape[0])
    print("test set has %d examples" % X_test.shape[0])
    
    model_in = Input(X_train.shape[1:])
    h = Conv2D(128, 3, activation='relu',padding="same", name='conv1_1')(model_in)
    h = Conv2D(64, 3, activation='relu',padding="same", name='conv1_2')(h)
    h = MaxPooling2D(pool_size=3, name='pool1',padding="same")(h)
    h = Dropout(0.15, name='drop1')(h)
    h = Flatten(name='flaten1')(h)
    
    dh = Conv2D(128, 2, activation='relu',dilation_rate=[1,3],padding="same", name='dconv1_1')(model_in)
    dh = Conv2D(64, 2, activation='relu',dilation_rate=[1,3],padding="same", name='dconv1_2')(dh)
    dh = MaxPooling2D(pool_size=2, name='dpool1')(dh)
    dh = Dropout(0.15, name='ddrop1')(dh)
    dh = Flatten(name='dflaten1')(dh)
    
    dh1 = Conv2D(128, 2, activation='relu',dilation_rate=[1,4],padding="same", name='dconv4_1')(model_in)
    dh1 = Conv2D(64, 2, activation='relu',dilation_rate=[1,4],padding="same", name='dconv4_2')(dh1)
    dh1 = MaxPooling2D(pool_size=2, name='d1pool1')(dh1)
    dh1 = Dropout(0.15, name='d1drop1')(dh1)
    dh1 = Flatten(name='d1flaten1')(dh1)
    
    h =  concatenate([h, dh, dh1])
    h = Dense(512, name="512dense", activation='relu')(h)
    h = Dropout(0.2, name='drop7')(h)
    h = Dense(128, name="last_dense", activation='relu')(h)
    h = Dropout(0.1, name='drop8')(h)
    output = Dense(4, name="out_dense", activation='softmax')(h)
    model = Model(inputs=[model_in], outputs=[output])
    
    model.compile(loss='categorical_crossentropy',
                  optimizer='adam',
                  metrics=['accuracy'])
    
    # define early stopping callback:
    earlystop = EarlyStopping(monitor='val_accuracy', min_delta=0.001, patience=patienceOption, verbose=1, mode='auto')
    
    model_json = model.to_json()
    with open(outputModel + ".json", "w") as json_file:
        json_file.write(model_json)
    modWeightsFilepath = outputModel + ".weights.hdf5"
    checkpoint = ModelCheckpoint(modWeightsFilepath, monitor='val_accuracy', verbose=1, save_best_only=True, save_weights_only=True, mode='auto')

    callbacks_list = [earlystop, checkpoint]

    datagen.fit(X_train)
    validation_gen.fit(X_valid)
    test_gen.fit(X_test)
    start = time.perf_counter()
    model.fit_generator(datagen.flow(X_train, Y_train, batch_size=32), \
                        steps_per_epoch=len(X_train) / 32, epochs=epochOption, verbose=1, \
                        callbacks=callbacks_list, \
                        validation_data=validation_gen.flow(X_valid,Y_valid, batch_size=32), \
                        validation_steps=len(X_valid) / 32)
    score = model.evaluate_generator(test_gen.flow(X_test, Y_test, batch_size=32), len(Y_test) / 32)
    sys.stderr.write("total time spent fitting and evaluating: %f secs\n" %(time.perf_counter() - start))
    
    print("evaluation on test set:")
    print("model loss: %f" % score[0])
    print("model accuracy: %f" % score[1])
    
    # print scores to file:
    with open(outputModel + ".scores.txt", "w") as score_file:
        score_file.write(str(score[0]) + '\t' + str(score[1]))
    
    # print predictions of the test set:
    if not argsDict['predictFile'] == None:
        preds = model.predict(test_gen.standardize(X_test))
        predictions = np.argmax(preds, axis=1)
        classDict = {0:'neutral', 1:'shared', 2:'divergent1', 3:'divergent2'}
        trues = np.argmax(Y_test, axis=1)
        with open(argsDict['predictFile'], 'w') as pred_file:
            pred_file.write('trueClass\tpredClass\tprob(neutral)\tprob(shared)\tprob(divergent1)\tprob(divergent2)\n')
            for index in range(len(preds)):
                pred_file.write('{}\t{}\t{:f}\t{:f}\t{:f}\t{:f}\n'.format(classDict[trues[index]], classDict[predictions[index]],
                                preds[index][0], preds[index][1], preds[index][2], preds[index][3]))
        print("{} predictions complete".format(index+1))



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('trainDir', help='path to training set files')
    parser.add_argument('testDir', help='path to test set files, can be same as trainDir')
    parser.add_argument('outputModel', help='file name for output model, will create two files one with structure one with weights')
    parser.add_argument("--equal", dest='equal', action='store_true', help="shuffle and cut data sets to equal length", default=False)
    parser.add_argument("--insist", dest='insist', action='store_true', help="insist on number of entries to keep", default=False)
    parser.add_argument('--keep_train', type=int, help='number of entries to keep for each training set', default=None)
    parser.add_argument('--keep_test', type=int, help='number of entries to keep for each test set', default=None)
    parser.add_argument('--epochs', type=int, help='max epochs for training CNN (default = 100)', default=100)
    parser.add_argument('--patience', type=int, help='patience during training for early stop (default = 5)', default=5)
    parser.add_argument('--numSubWins', type=int, help='number of subwindows that our chromosome is divided into (default = 11)', default=11)
    parser.add_argument('--predictFile', type=str, help='file to write the predictions of the test data set', default=None)
    main(vars(parser.parse_args()))

