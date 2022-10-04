#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import joblib
import numpy as np
import pandas as pd

def SLEder(path_to_model, path_to_mean, path_to_var, path_to_file):
    '''
    SLEDER is a programm for atherosclerosis plug prediction
    by metabolomic data from patients with SLE
    
    path_to_mode - model which will be used for prediction
    path_to_mean, path_to_var - scaling parametrs
    path_to_file - file from patients with target metabolites
    
    '''
    # Download model
    model = joblib.load(path_to_model)
    # load metabolites information
    box = []
    with open(path_to_file, 'r') as file:
        for line in file:
            box.append(line.strip())
    name, values = box[0], box[-1].split("\t")
    # transform matabolic values to using
    values = np.array(list(map(lambda x: float(x), 
                           list(map(lambda x: (".").join(x.split(",")), values)))))
    # load mean values for scaling
    common_mean = []
    with open(path_to_mean, 'r') as file:
        for line in file:
            common_mean.append(line.strip())
    common_mean = np.array(list(map(lambda x: float(x), common_mean)))
    
    # load varience values for scaling
    common_var = []
    with open(path_to_var, 'r') as file:
        for line in file:
            common_var.append(line.strip())
    common_var = np.array(list(map(lambda x: float(x), common_var)))
    
    # Scaling, prediction, prediction confidence
    values_sc = (values-common_mean)/np.sqrt(common_var)
    prediction = model.predict(np.reshape(values_sc, (1, 227)))
    confidence = round(max(model.predict_proba(np.reshape(values_sc, (1, 227)))[0]), 2)*100
    
    # Get prognosis
    # print(name, end='\n')
    # print("Plug+ ; Confidence: ", confidence, "%") if prediction[0] == 1 else print("Plug- ; Confidence: ", confidence, "%")
    if prediction[0]:
        answer = "PLUG POSITIVE"
    else:
        answer = "PLUG NEGATIVE"

    return name, answer, confidence
