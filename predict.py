import numpy as np
import pandas as pd
from keras.models import Sequential
from keras.layers import Dense
from keras import backend as K
from keras.models import load_model


# K.tensorflow_backend._get_available_gpus()

from sklearn.model_selection import train_test_split
import argparse



def predict(input_file, SYBYL, ECFP, ECFP_SYBYL, ECFP_num):
    
    data = pd.read_csv(input_file)
    data_feature = data.drop(["SMILES"], axis=1)
    
  

    if SYBYL:

        model_IE = load_model('model/SYBYL_IE.h5')
        model_EA = load_model('model/SYBYL_EA.h5')
        
        IE = model_IE.predict(data_feature.values)
        EA = model_EA.predict(data_feature.values)
        data['IE'] = IE
        data['EA'] = EA
        data[['SMILES','IE','EA']].to_csv("SYBYL_predict.csv", index=False)
    elif ECFP:

        model_IE = load_model('model/ECFP_IE.h5')
        model_EA = load_model('model/ECFP_EA.h5')

        IE = model_IE.predict(data_feature.values)
        EA = model_EA.predict(data_feature.values)
        
        data['IE'] = IE
        data['EA'] = EA
        data[['SMILES','IE','EA']].to_csv("ECFP_predict.csv", index=False)
        
    
    
    elif ECFP_SYBYL:

        model_IE = load_model('model/ECFP_SYBYL_IE.h5')
        model_EA = load_model('model/ECFP_SYBYL_EA.h5')
        
        IE = model_IE.predict(data_feature.values)
        EA = model_EA.predict(data_feature.values)
        
        data['IE'] = IE
        data['EA'] = EA
        data[['SMILES','IE','EA']].to_csv("ECFP_SYBYL_predict.csv", index=False)

    elif ECFP_num:
      
        model_IE = load_model('model/ECFP_num_IE.h5')
        model_EA = load_model('model/ECFP_num_EA.h5')
        
        IE = model_IE.predict(data_feature.values)
        EA = model_EA.predict(data_feature.values)

        data['IE'] = IE
        data['EA'] = EA
        data[['SMILES','IE','EA']].to_csv("ECFP_num_predict.csv", index=False)

        
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--data_file',
                        help='the file including smiles and its representation', default=None)
    parser.add_argument('-m', '--model',
                        help='predict model', default=None)
   
    
    
    args = vars(parser.parse_args())
    # train model 
    if args['model'] == 'ECFP':
        predict(input_file = args['data_file'], ECFP=True, SYBYL=None, ECFP_num=None, ECFP_SYBYL=None)
    elif args['model'] == 'SYBYL':
        predict(input_file = args['data_file'], ECFP=None, SYBYL=True, ECFP_num=None, ECFP_SYBYL=None)
    elif args['model'] == 'ECFP_SYBYL':
        predict(input_file = args['data_file'], ECFP=None, SYBYL=None, ECFP_num=None, ECFP_SYBYL=True)
    elif args['model'] == 'ECFP_num':
        predict(input_file = args['data_file'], ECFP=None, SYBYL=None, ECFP_num=True, ECFP_SYBYL=None)