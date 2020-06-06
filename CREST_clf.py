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
    
    Parameters:
    bandwidth : float
    kernel : str
    """
    def __init__(self, bandwidth=1.0, kernel='gaussian'):
        self.bandwidth = bandwidth
        self.kernel = kernel
        
    def fit(self, X, y,bandwidths):
        self.classes_ = np.sort(np.unique(y))
        training_sets = [X[y == yi] for yi in self.classes_]
        
        params= [{'bandwidth': bandwidths, 'kernel': ['gaussian']}
         ]
        self.grid_ = GridSearchCV(KernelDensity(),params, cv=5)
        self.models_=[]
        for Xi in training_sets:
            self.grid_.fit(Xi)
            print (self.grid_.best_estimator_.bandwidth)
            self.models_.append(KernelDensity(bandwidth=self.grid_.best_estimator_.bandwidth).fit(Xi))
            
        
        #self.models_ = [KernelDensity(bandwidth=best_params[tuple(Xi)]).fit(Xi) for Xi in training_set]
        
        #self.models_ = [KernelDensity(bandwidth=self.bandwidth, kernel=self.kernel).fit(Xi)for Xi in training_sets]
                                     
        self.logpriors_ = [np.log(Xi.shape[0]*1.0 / X.shape[0])
                           for Xi in training_sets]
        return self
        
    def predict_proba(self, X):
        logprobs = np.array([model.score_samples(X)
                             for model in self.models_]).T
        result = np.exp(logprobs + self.logpriors_)
        return result / result.sum(1, keepdims=True)
        
    def predict(self, X):
        return self.classes_[np.argmax(self.predict_proba(X), 1)]


def read_in_data(file，total_len):
    id1, id2, c, r1,r2 = np.genfromtxt(file, delimiter=',', usecols=(0,1,2,3,4), unpack=True)
    #np.genfromtxt(StringIO(data), usecols=(0, -1))
    f1 = [min(r1[i],r2[i]) for i in range(len(y))]
    f2 = [max(r1[i],r2[i]) for i in range(len(y))]
    coverage = [c/total_len for i in range(len(y))]
    d = [np.log(r1[i])-np.log(r2[i]) if r1[i]*r2[i]!=0 else 0 for i in range(len(y))]
    features = np.column_stack((id1, id2, f1,f2,coverage,d))
    
    return features

def prediction_KDE(features, start, end, inv, class_models, direction_models):
    
    y_real = {}
    cm_dict = {}
    
    y_pred = {}
    pred_prob = {}

    total_num = len(coverage)
    bins = {}
    
    for coverage in range(start,end,inv):
        if coverage < end - inv:
            index_bin = [i for i in range(total_num) if ~np.isnan(features[i,2]) and  features[i,0]!=0 and coverage * 1.0 / 100 <= features[i,2] <= (coverage+inv)*1.0/100]
            #f1_bin = [features[i,0] for i in range(total_num) if features[i,2] >= coverage * 1.0 / 100 and features[i,0] <= (coverage+inv)*1.0/100]
            #f2_bin = [features[i,1] for i in range(total_num) if features[i,2] >= coverage * 1.0 / 100 and features[i,0] <= (coverage+inv)*1.0/100]
        else:
            index_bin = [i for i in range(total_num) if ~np.isnan(features[i,2]) and features[i,0]!=0 and features[i,2] >= coverage * 1.0 / 100 ]
            #f1_bin = [features[i,0] for i in range(total_num) if features[i,2] >= coverage * 1.0 / 100]
            #f2_bin = [features[i,1] for i in range(total_num) if features[i,2] >= coverage * 1.0 / 100]
        if len(index_bin)==0:
            continue
        index1 = np.random.choice([i for i in index_bin if features[i,3]==1],1500)
        index2 = np.random.choice([i for i in index_bin if features[i,3]==2],1500)
        index3 = np.random.choice([i for i in index_bin if features[i,3]==3],1500)
        
        index = np.concatenate((index1, index2, index3))
       
        
        X_pred = features[index,0:2]
        #X_pred[:,0] = np.divide(X_pred[:,0],X_pred[:,1])
        y_pred[coverage] = models[coverage].predict(X_pred)
        pred_prob[coverage] = models[coverage].predict_proba(X_pred)
        #model.fit(X_train, y_train)
    #if true:
    #    pred_prob[coverage] = 
    #    y_p[coverage] = p 
    
    #return cm_dict, pred_prob, y_p, fp_gp, fp_ac, fp_hs, fn_gp, fn_ac, fn_hs
    return id1,id2, y_pred, cm_dict, pred_prob

def write_results():
    
    return


def main(args):

    bandwidths = np.logspace(-2, -0.5, 30)

    features = read_in_data('csvfiles/train_test_10k.csv')
    test_model = pickle.load(open('models.pl', 'rb'))
    id1,id2,y_true, y_pred, cm, pred_prob = prediction_KDE(test_model,'train.csv',5,45,5)

    
	        write_results(ibd_union, os.path.join(args.out_dirpath, 'fc_'+str(i)+'_ibd.gz'))
	    

if __name__ == '__main__':
    parser = ArgumentParser('Compare HAPI and genotype.')
    # parser.add_argument('--fblock-path',
    # type=bool,
    # default=0,
    #                    type=str,
    #                    help='path for fblocka.')

    parser.add_argument('--fam-path',
                        type=str,
                        help='Path for fam file from pedigrees.')

    parser.add_argument('--out-dirpath',
                        type=str,
                        help='Path to output files.')

    args = parser.parse_args()
    main(args)
