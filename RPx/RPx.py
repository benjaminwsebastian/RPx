import os, sys, datetime
import numpy as np
import scipy.stats as stats
import pandas as pd
import _thread

class detector(object):

    def __init__(self, file_name, period, n_cycles, **kwargs):
        self.config = {
            "EPS"                      : 1e-5,
            "USE_Z_SCORE"              : False,
            "NUM_PERMUTATIONS"         : 5000,
            "NUM_CORES"                : 1,
            "SHUFFLE_WITH_REPLACEMENT" : False,
            "MIN_RP24"                 : 0.0,
            "FILTER"                   : False
        }

        for kwarg in self.config:
            try:
                self.config[kwarg] = kwargs[kwarg]
            except:
                pass
                
        try:
            self.config["FILE_NAME"] = file_name
            self.config["PERIOD"] = period
            self.config["N_CYCLES"] = n_cycles
        except KeyError:
            raise KeyError("Missing required argument\n Required arguments: file_name, period, n_cycles")

    def __repr__(self):
        return '<File {}, Period {}, n Cycles {}, num cores {}, filter {}'.format(self.config["FILE_NAME"], self.config["PERIOD"], self.config["N_CYCLES"], self.config["NUM_CORES"], self.config["FILTER"])

    ####################
    ## PRE-PROCESSING ## 
    ####################

    def read_file(self):
        df = pd.read_csv(self.config["FILE_NAME"])
        df["detectable"] = df[1::].apply(detectable)
        df["zero"] = df[1::].apply(contains_zero)

        return df

    ###########################
    ## COMPUTE RPx FUNCTIONS ## Refactor, implement with pandas, and remove gene specific language
    ###########################

    def compute_phase(row):
        components = np.fft.fft(row)
        phi = np.angle(complex(components[2].real,components[2].imag))
        phase = -12/np.pi*phi
        if phase < 0: phase = 24 + phase

        return phase

    def compute_rp24(col):
        col = list(col)
        DFT_mod = (abs(np.fft.fft(col)) / len(col))**2
        # first step to getting the Power Spectral Density (PSD), 
        # the DFT_mod is just a modified Distcrete Fourier Transfor at this stage
    
        PSD_mod = [0] * len(DFT_mod)
        for i in range(1, int(len(DFT_mod) / 2) + 1):
            if i % len(DFT_mod/2) != 0:
                PSD_mod[i] = 2 * DFT_mod[i]
            else:
                PSD_mod[i] = DFT_mod[i]
                
        RP24 = np.log2(PSD_mod[self.config["PERIOD"]] / (sum(PSD_mod[ 0:int(len(PSD_mod)/2) ]) - PSD_mod[self.config["PERIOD"]]))

        return RP24

    def compute_pvalue(expression, PSD, RP24):
        # compute the p-value for a single gene
        if RP24 >= self.config["MIN_RP24"]:
            # first compute random shuffles of expression
            RP24s_shuffled = []    
            temp = list(expression)
            for i in range(self.config["NUM_PERMUTATIONS"]):
                if self.config["SHUFFLE_WITH_REPLACEMENT"]:
                    temp1 = np.random.choice(temp,len(temp))
                    RP24s_shuffled.append(compute_rp24(temp1))
                else:
                    np.random.shuffle(temp)
                    RP24s_shuffled.append(compute_rp24(temp))
            # next compute the p-value
            if self.config["USE_Z_SCORE"]:
                # Calculate p-value using a zscore and normal distribution
                z_score = (RP24-np.mean(RP24s_shuffled))/np.std(RP24s_shuffled)
                return stats.norm.sf(z_score)
            else:
                #Calculate empirical p-value
                count = 0.0        
                for value in RP24s_shuffled:
                    if value >= RP24:
                        count+= 1.0
                return float(count)/len(RP24s_shuffled)
        else:
            # if RP24 negative, then don't test RP24 significance
            return "NOTEST"

    def preprocess_pvalues(expression, detectable):
        expression_cleaned = {}
        for gene in expression:
            if detectable[gene][1] != 0:
                expression_cleaned[gene] = expression[gene]
        return expression_cleaned

    def compute_qvalues(pvalues):
        pvalues_unsorted = []
        for id in pvalues:
            if pvalues[id] != "NOTEST":
                pvalues_unsorted.append((id,pvalues[id]))

        # Sort list by p value
        pvalues_sorted = sorted(pvalues_unsorted, key=lambda p:p[1])
    
        # calculate q values from sorted list, going backwards so q values can be adjusted if needed
        qprev = 1.0
        qvalues = {}
        for i in reversed(range(len(pvalues_sorted))):
            id, pvalue = pvalues_sorted[i]
            r=float(i+1.0)
            qtemp=pValue*len(pvalues_sorted)/r
            qvalue=min(qtemp,qprev)
            qprev=qvalue
            qvalues[id] = qvalue

        for id in pvalues:
            if pvalues[id] == "NOTEST":
                qvalues[id] = "NOTEST"

        return qvalues
        
    ###############
    ## DETECTION ## Implement this with pandas
    ############### 
    '''
    def detect(self):
        for gene in expression_cleaned:
            PSDs[gene] = computePSD(expression_cleaned[gene])
            RP24s[gene] = compute_rp24(PSDs[gene])
            RP24_Pvalues[gene] = compute_pvalue(expression_cleaned[gene],PSDs[gene],RP24s[gene],params)
            phases[gene] = compute_phase(expression[gene])
            medians[gene] = np.median(expression[gene])

            RP24_Qvalues = computeQvalues(RP24_Pvalues)
    '''
    ######################
    ## HELPER FUNCTIONS ##
    ######################

    def detectable(values):
        if np.median(values) >= 1.0 and (max(values) + self.config["EPS"]) / (min(values) + self.config["EPS"]) > 1.5:
            return 1
        else:
            return 0

    def contains_zero(values):
        if 0 in values:
            return 0
        else:
            return 1
            
    
