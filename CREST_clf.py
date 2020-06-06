__author__ = 'qiaoqiao'
import sys
import os
import numpy as np
import csv
import math 
import random
from sklearn.neighbors.kde import KernelDensity
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import GridSearchCV
from scipy.stats import multivariate_normal
import pickle
from sklearn.model_selection import train_test_split 
from sklearn.metrics import precision_recall_curve
from yellowbrick.classifier import PrecisionRecallCurve
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
        
    def fit(self, X, y,bandwidths, prior = 0):
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
        if prior == 0:                 
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
    id1, id2, c, r1,r2 = np.genfromtxt(file, delimiter=',', usecols=(0,1,2,3,4), unpack=True)
    #np.genfromtxt(StringIO(data), usecols=(0, -1))
    f1 = [min(r1[i],r2[i]) for i in range(len(c))]
    f2 = [max(r1[i],r2[i]) for i in range(len(c))]
    coverage = [c/total_len for i in range(len(c))]
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
            index_bin = [i for i in range(total_num) if ~np.isnan(features[i,2]) and features[i,3] != 0 and start + index_win * inv <= features[i,2] <= start + (index_win + 1) * inv]
              
        else:
            index_bin = [i for i in range(total_num) if ~np.isnan(features[i,2]) and features[i,3]!=0 and features[i,2] >= start + index_win * inv ]


       # print(len(index))
       # print(features.shape)
        #features_bin = features[index,:]
       
        if len(index_bin) == 0:
            print("There is no enough data in this bin")
            continue
        if len(index_bin) == 1:
            X_train = features[index_bin,3:5].reshape(-1,1)
        else:
            X_train = features[index_bin,3:5]

        #X_train[:,0] = np.divide(X_train[:,0],X_train[:,1])
        Y_train = labels[index_bin]
        #bandwidths = np.logspace(-2, 1, 20)
        model.fit(X_train, Y_train,bandwidths)
        models[index_win] = model
        pickle.dump(models, open(model_file, 'wb'))
    return models 

# how to allow customize start/ end/ inv? 
def prediction_KDE(features, start, end, inv, class_models, direction_model,prior=0):

    total_num = features.shape[0]
    results = numpy.empty(total_num, 9)
    resutls[:,0:2] = features[:,0:2]
    
    for i in range(total_num):
        if ~np.isnan(features[i,2]) and float(features[i,2]) >= start:
            coverage = float(features[i,2])
            model_index = int((coverage-start)/inv)
            X_pred = features[i,3:5].reshape(-1,1)
        #X_pred[:,0] = np.divide(X_pred[:,0],X_pred[:,1])
            y_pred = class_models[model_index].predict(X_pred)
            pred_prob = class_models[model_index].predict_proba(X_pred)
            
            direction_pred = direction_model.predict(features[i,5].reshape(-1,1))
            pred_prob1 = direction_model.predict_proba(features[i,5].reshape(-1,1))
            results[i,2] = y_pred
            results[i,3:6] = pred_prob
            results[i,6] = direction_pred
            results[i,7:9] = pred_prob1

        else: 
            results[i,2] = 0
            if prior == 0:
                results[i,3:6] = [1/3, 1/3, 1/3]
            else:
                results[i,3:6] = prior[i,:]
            results[i,6:9] = [0,0,0]
            
    return results


def main(args):
    features = read_in_data(args.input)
    if args.train:
        #read in labels 
        with open(args.labels) as csvfile:
            labels = csv.reader(csvfile, delimiter=',')
        bandwidths = np.logspace(-2, -0.5, 30)

    test_model = pickle.load(open(args.models, 'rb'))
    direction_model = args.models_direction

    results = prediction_KDE(features, start, end, inv, test_models, direction_model)
    np.savetxt(args.output, results, delimiter=',')
      

if __name__ == '__main__':
    parser = ArgumentParser('CREST classifier to predict relationships')
    # parser.add_argument('--fblock-path',
    # type=bool,
    # default=0,
    #                    type=str,
    #                    help='path for fblocka.')
    parser.add_argument('--train',
                        help='To train new models with training data and labels.',
                        action="store_true")
    parser.add_argument('--start',
                        type=float, default = 0.025,
                        help='The minimum coverage to consider.')
    parser.add_argument('--end',
                        type=float, default = 0.2,
                        help='The coverage to start to merge models.')
    parser.add_argument('--inv',
                        type=float,default = 0.025,
                        help='The interval of coverage to train models.')
    parser.add_argument('--prior',
                        type=str, 
                        help='File to give prior probability.')
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help='Path for input data.')
    parser.add_argument('--models',
                        type=str, required=True,
                        help='Path for trained models to predict relationship types.')
    parser.add_argument('--models_direction',
                        type=str,
                        help='Path for trained models to predict directionality.')
    parser.add_argument('-o','--output',
                        type=str,
                        help='Path for output results.')
    parser.add_argument('--labels',
                        type=str,
                        help='File with labels for trainning data.')

    args = parser.parse_args()
    main(args)
