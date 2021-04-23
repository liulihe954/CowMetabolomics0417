# project: p7
# submitter: lliu356
# partner: none
# hours: 4

import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures, OneHotEncoder
from sklearn.compose import make_column_transformer
import netaddr
import copy
import os

class UserPredictor():
    
    def __init__(self):
        self.model = Pipeline([('both',make_column_transformer((OneHotEncoder(),['badge']),
                                                               (OneHotEncoder(handle_unknown = 'ignore'),['region']),
                                                               (PolynomialFeatures(degree = 2,include_bias=False),['weighted_time']),
                                                               remainder = "passthrough")),
                               ('cls', LogisticRegression(fit_intercept=False))])       
        self.features_cols = ['region','badge','age','past_purchase_amt','weighted_time']
        self.y = 'y'
        
    def process_raw(self,users,logs, y = None):
        ### user:
        user_tmp = users #.copy()
        
        ### logs: match region
        ip2check_test = list(logs['ip_address'])
        ip_df = self.ip2location_load()
        match_region_out = self.ip_check(ip2check_test,ip_df)
        logs_add_region = logs #.copy()
        logs_add_region['region'] = match_region_out
        logs_match_region = logs_add_region[['id','region']].drop_duplicates()
        
        ### logs: match weigted time 
        page_type = list(logs['url_visited'])
        
        ### trnsform url (manually)
        url_transform = []
        #url_ref = {}
        for i in range(0,len(page_type)):
            tmp = page_type[i].replace('.html','').replace('/','')
            url_transform.append(tmp)
            #url_ref[tmp] = url_ref.setdefault(tmp,0) + 1
        
        # Any better idea than hardcoding these categories...?
        laptop = ['laptop']; office = ['tablet','keyboard','monitor','printer','desk']
        
        ### weight time
        page_type = list(logs['url_visited'])
        minutes_raw = list(logs['minutes_on_page'])
        minutes_weighted = [None] * len(minutes_raw)
        
        for i in range(0,len(page_type)):
            tmp = page_type[i].replace('.html','').replace('/','')
            if tmp in laptop:
                minutes_weighted[i] = minutes_raw[i] * 0.7
            elif tmp in office:
                minutes_weighted[i] = minutes_raw[i] * 0.2
            elif tmp == 'NotFound':
                minutes_weighted[i] = 0
        # extract        
        logs_match_weight = logs#.copy()
        logs_match_weight['weighted_time'] = minutes_weighted
        logs_match_weight = pd.DataFrame(logs_match_weight.groupby('id')['weighted_time'].sum())
        
        ### merge:
        user_wtime = pd.merge(user_tmp,logs_match_weight,on = 'id',how='left').fillna(0)
        out = pd.merge(user_wtime,logs_match_region,on = 'id',how='left').fillna('NotFound')
        if not y is None:
            out = pd.merge(out,y,how='left',on = 'id')
            out['y'] = out['y'].map({0:False,1:True})
            return out
        else:
            return out

    def fit(self, train_users, train_logs, train_y, crs_val = False):
        
        self.features_train = self.process_raw(train_users, train_logs, train_y)
        
        self.model.fit(self.features_train[self.features_cols],self.features_train[self.y])
        
        if crs_val:
            crs_val_scores  = cross_val_score(self.model,
                                              self.features_train[self.features_cols],
                                              self.features_train[self.y])
            print(f"AVG: {crs_val_scores.mean()}, STD: {crs_val_scores.std()}\n")
            
    def predict(self,test_users, test_logs):
        self.features_test = self.process_raw(test_users, test_logs)
        y_predicted = self.model.predict(self.features_test[self.features_cols])
                       
        return np.asarray(list(map(int, y_predicted)))
    
    def BinarySearch(self,arr,x):
        left, right = 0, len(arr) - 1
        while left < right:
            mid = left + (right - left) // 2
            if x - arr[mid] > arr[mid + 1] - x:
                left = mid + 1
            else:
                right = mid
        if arr[left-1] <= x and arr[left] > x:
            return left - 1
        else:
            return left
                       
    def ip2location_load(self):
        with open(os.path.join('data','ip2location.csv')) as f:
            out_raw = f.read()
        out = out_raw.split("\n")
        ip_raw = []
        for line in out:
            line_tmp = line.split(",")
            ip_raw.append(line_tmp)    
        ip_df = pd.DataFrame(ip_raw)
        ip_df.rename(columns = ip_df.iloc[0], inplace = True)
        ip_df.drop(ip_df.index[0],inplace = True)
        ip_df.drop(ip_df.index[-1],inplace = True) # something weird
        ip_df.sort_values(by = ["low"]) # sort for binary search
        ip_df.reset_index(drop = True,inplace = True)
        return ip_df
    
    def ip_check(self,ips,ip_df):
        search_sequence = list(map(int, ip_df["low"]))
        ip_match_out = []
        for ip in ips:    
            try:
                int_ip = int(netaddr.IPAddress(ip))
                index_tmp = self.BinarySearch(search_sequence,int_ip)
                matched_region = ip_df.iloc[index_tmp,3]
                ip_match_out.append(matched_region)  
            except:
                print("Execution halted...")
                break
        return ip_match_out