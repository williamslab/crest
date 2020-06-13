#!/usr/bin/env python3
import sys
import os
import numpy as np
import csv
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import GridSearchCV
import pickle
from sklearn.base import BaseEstimator, ClassifierMixin
from argparse import ArgumentParser

# TODO: the kernel as parameter or fixed
class KDEClassifier(BaseEstimator, ClassifierMixin):
    """
    classifier based on KDE

    """
    def __init__(self, bandwidth=0.1, kernel='gaussian'):
        self.bandwidth = bandwidth
        self.kernel = kernel
        
    def fit(self, X, y,bandwidths):
        self.classes_ = np.sort(np.unique(y))
        training_sets = [X[y == yi] for yi in self.classes_]
        #
        params= [{'bandwidth': bandwidths, 'kernel': ['gaussian','tophat','exponential','linear']}]
        self.grid_ = GridSearchCV(KernelDensity(),params, cv=5)
        self.models_=[]
        for Xi in training_sets:
            self.grid_.fit(Xi)
            #print (self.grid_.best_estimator_.bandwidth)
            self.models_.append(KernelDensity(bandwidth=self.grid_.best_estimator_.bandwidth, kernel=self.grid_.best_estimator_.kernel).fit(Xi))
            
        
        #self.models_ = [KernelDensity(bandwidth=best_params[tuple(Xi)]).fit(Xi) for Xi in training_set]
        
        #self.models_ = [KernelDensity(bandwidth=self.bandwidth, kernel=self.kernel).fit(Xi)for Xi in training_sets]
                        
            self.logpriors_ = [np.log(Xi.shape[0]*1.0 / X.shape[0]) for Xi in training_sets]
        return self
        
    def predict_proba(self, X):
        logprobs = np.array([model.score_samples(X)
                             for model in self.models_]).T
        
        result = np.exp(logprobs + self.logpriors_)
        return result / result.sum(1, keepdims=True)
        
    def predict(self, X):
        return self.classes_[np.argmax(self.predict_proba(X), 1)]


def read_in_data(file, total_len):
    id1 = [ ]
    id2 = [ ]
    c = [ ]
    r1 = [ ]
    r2 = [ ]
    with open (file) as file:
        for lines in file:
            line = lines.split(',')
            id1.append(line[0])
            id2.append(line[1])
            c.append(float(line[2]))
            r1.append(float(line[3]))
            r2.append(float(line[4]))



    #id1, id2, c, r1,r2 = np.genfromtxt(file, delimiter=',', usecols=(0,1,2,3,4), unpack=True)
    #np.genfromtxt(StringIO(data), usecols=(0, -1))
    f1 = [min(r1[i],r2[i]) for i in range(len(c))]
    f2 = [max(r1[i],r2[i]) for i in range(len(c))]
    coverage = [c[i]/(total_len*2) for i in range(len(c))]
    d = [np.log(r1[i])-np.log(r2[i]) if r1[i]*r2[i]!=0 else 0 for i in range(len(c))]
    features = np.column_stack((id1, id2, coverage, f1, f2, d))
    
    return features

def train_KDE_model(start, end, inv, features, labels, model_file, bandwidths):
    '''
    use bins to store data 
    for loop once
    train model separately
    store one or more trained models?  
    '''
    
    #total_len = 3346.299

    total_num = features.shape[0]
    models = {}
    num_win = round((end - start)/inv)+1
    for index_win in range(num_win):
        if index_win < num_win - 1 :
            index_bin = [i for i in range(total_num) if ~np.isnan(float(features[i,2])) and float(features[i,3]) != 0 and start + index_win * inv <= float(features[i,2]) <= start + (index_win + 1) * inv]
              
        else:
            index_bin = [i for i in range(total_num) if ~np.isnan(float(features[i,2])) and float(features[i,3])!=0 and float(features[i,2]) >= start + index_win * inv ]


       # print(len(index))
       # print(features.shape)
        #features_bin = features[index,:]
       
        if len(index_bin) == 0:
            print("There is no enough data in this bin")
            continue
        if len(index_bin) == 1:
            X_train = features[index_bin,3:5].reshape(-1,1).astype(np.float)
        else:
            X_train = features[index_bin,3:5].astype(np.float)

        #X_train[:,0] = np.divide(X_train[:,0],X_train[:,1])
        Y_train = labels[index_bin]
        #bandwidths = np.logspace(-2, 1, 20)
        model.fit(X_train, Y_train,bandwidths)
        models[index_win] = model
        pickle.dump(models, open(model_file, 'wb'))
    return models 

# how to allow customize start/ end/ inv? 
def prediction_KDE(features, start, end, inv, class_models, direction_model,prior):

    total_num = features.shape[0]
    results = np.empty([total_num, 10], dtype=object)
    results[:,0:2] = features[:,0:2]
    dict_type = {0:'UN',1:'GP',2:'AV',3:'HS'}
    
    for i in range(total_num):
        if ~np.isnan(float(features[i,2])) and float(features[i,2]) >= start:
            coverage = float(features[i,2])
            if coverage >= end:
                model_index = int((end-start)/inv)
            else:
                model_index = int((coverage-start)/inv)
            X_pred = features[i,3:5].astype(np.float).reshape(1,-1)
        #X_pred[:,0] = np.divide(X_pred[:,0],X_pred[:,1])
            y_pred = class_models[model_index].predict(X_pred)
            pred_prob = class_models[model_index].predict_proba(X_pred)
            
            direction_pred = direction_model.predict(features[i,5].reshape(1,-1))
            pred_prob1 = direction_model.predict_proba(features[i,5].reshape(1,-1))

            

            results[i,2] = int(y_pred[0])
            results[i,3] = dict_type[results[i,2]]
            results[i,4:7] = (pred_prob[0] * prior)/ sum(pred_prob[0] * prior)
            results[i,7] = int(direction_pred[0])
            results[i,8:10] = pred_prob1[0]

        else: 
            results[i,2] = 0
            results[i,3] = dict_type[results[i,2]]
            results[i,4:7] = prior[:]
            results[i,7:10] = [0,0,0]
            
    return results


def main(args):
    
    features = read_in_data(args.input,args.total_len)

    if abs(sum(args.prior)-1) >0.01:
        print("The sum of prior probability is not 1, will renormalize it. ")
    prior = np.asarray(args.prior)/sum(args.prior)
    if args.train:
        #read in labels 
        with open(args.labels) as csvfile:
            labels = csv.reader(csvfile, delimiter=',')
        bandwidths = np.logspace(-2, -0.5, 30)
        test_model = train_KDE_model(args.start, args.end, args.inv, features, labels, args.models_type, bandwidths)
    else:

        test_model = pickle.load(open(args.models_type, 'rb'))

    direction_model = pickle.load(open(args.models_direction, 'rb'))

    results = prediction_KDE(features, args.start, args.end, args.inv, test_model, direction_model,prior)
    output = args.output + '.csv'
    with open(output, 'wb') as f:
        f.write(b'ID1,ID2,Class,Type,Prob_GP,Prob_AV,Prob_HS,Direction,Prob1,Prob2\n')
        np.savetxt(f, results, delimiter=',',fmt='%s')
      

if __name__ == '__main__':
    path = os.path.dirname(__file__)+'/'

    parser = ArgumentParser('CREST_relationships.py')
    
    parser.add_argument('-i', '--input',
                        type=str, required = True,
                        help='File with input data.')
    parser.add_argument('--total_len', 
                        type=float, default = 3536.5466,
                        help='The total length of genome in cM.')
    parser.add_argument('--models_type',
                        type=str, default = path + 'type_clf.pickle',
                        help='File of trained models to predict relationship types.')
    parser.add_argument('--models_direction', 
                        type=str, default = path + 'direction_clf.pickle',
                        help='File of trained models to predict directionality.')
    parser.add_argument('-o','--output',
                        type=str, default = 'relationships',
                        help='File to output results.')
    parser.add_argument('--start',
                        type=float, default = 0.025,
                        help='The minimum coverage rate to classify.')
    parser.add_argument('--end',
                        type=float, default = 0.20,
                        help='The coverage rate to start to merge models.')
    parser.add_argument('--inv',
                        type=float,default = 0.025,
                        help='The window size of coverage for each model.')
    parser.add_argument('--prior',
                        nargs = "*", type=float, default=[1/3,1/3,1/3],
                        help='Prior probability of three types.')
    parser.add_argument('--train',
                        help='To train new classifier with training data and labels.',
                        action="store_true")
    parser.add_argument('--labels',
                        type=str,
                        help='File of labels for trainning data.')
    
    versionFile = path + "version.h"
    if os.path.exists(versionFile):
        with open(versionFile) as f:
            for lines in f:
                line = lines.split( )
                if line[1] == "VERSION_NUMBER":
                    version = line[2].replace('"','')
                if line[1] == "RELEASE_DATE":
                    date = line[2].replace('"','') + ' ' + line[3] + ' ' + line[4].replace('"','')
        print("\nCREST  v" + version + "\n" + "Released " + date +"\n\n")
    else:
        print("Please download version.h file to get the version information.")

    args = parser.parse_args()
    main(args)
