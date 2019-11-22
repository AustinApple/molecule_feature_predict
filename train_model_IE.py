import numpy as np
import pandas as pd
from rdkit.Chem import AllChem
from rdkit import Chem, DataStructs
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn import utils
from sklearn.metrics import mean_squared_error, mean_absolute_error 
from sklearn.ensemble import RandomForestRegressor
from sklearn.externals import joblib
from keras.models import Sequential
from keras.layers import Dense
from keras import backend as K
# K.tensorflow_backend._get_available_gpus()

from sklearn.model_selection import train_test_split
import argparse



def train(input_file, save_model):
    data = pd.read_csv(input_file)
    X_input = data.drop(['SMILES','EA','IE'],axis=1)
    y_EA_IE= data[['IE']]
    X_train, X_test, y_train, y_test = train_test_split(X_input.values, y_EA_IE.values, test_size=0.1, random_state=1)
    
    X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=0.1, random_state=1)
    # building model
    model = Sequential()
    model.add(Dense(output_dim=int(X_train.shape[1]/2), input_dim=X_train.shape[1],activation='relu'))
    model.add(Dense(output_dim=int(X_train.shape[1]/6),activation='relu'))
    model.add(Dense(output_dim=int(X_train.shape[1]/12),activation='relu'))
    model.add(Dense(output_dim=y_train.shape[1]))
    model.compile(loss='mae', optimizer='adam',metrics=['mae'],)
    print('Training -----------')
    model.fit(X_train, y_train, verbose=1, epochs=100, validation_data=[X_val, y_val])
    
    print("test testing set")
    print(model.evaluate(X_test, y_test))
    

# prediction = mlp.predict(X_test)
# print(mean_absolute_error(y_test, prediction))
# print(mean_squared_error(y_test, prediction))

    # save model
    model.save(save_model)	


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--data_file',
                        help='the file including smiles and corresponding IE and EA', default=None)
    parser.add_argument('-o', '--save_model',
                        help='the name of save model', default=None)
   
    
    
    args = vars(parser.parse_args())
    # train model 
    train(args['data_file'], args['save_model'])

    