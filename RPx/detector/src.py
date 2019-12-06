import os
import sys
import datetime
import numpy as np
import scipy.stats as stats
import pandas as pd
import multiprocessing

class detector(object):

    def __init__(self, file_name, period, n_cycles, **kwargs):
        self.config = {
            "EPS"                      : 1e-5,
            "STAGGER"                  : False,
            "USE_Z_SCORE"              : False,
            "NUM_PERMUTATIONS"         : 5000,
            "NUM_PROCS"                : 1,
            "SHUFFLE_WITH_REPLACEMENT" : False,
            "MIN_RP24"                 : 0.0,
            "FILTER_DETECTABLE"        : False,
            "FILTER_ZERO"              : False,
            "LEN_SIGNALS"              : None,
        }

        for kwarg in self.config:
            try:
                self.config[kwarg] = kwargs[kwarg.lower()]
            except:
                pass
                
        try:
            self.config["FILE_NAME"] = file_name
            self.config["PERIOD"] = period
            self.config["N_CYCLES"] = n_cycles
        except KeyError:
            raise KeyError("Missing required argument\n Required arguments: file_name, period, n_cycles\n note that n_cycles is the same as number of replacates in circadian research")

    def __repr__(self):
        return '<File {}, Period: {}, n_Cycles: {}, num_procs: {}, filter_detectable: {}, filter_zero: {}>'.format(self.config["FILE_NAME"], self.config["PERIOD"], self.config["N_CYCLES"], self.config["NUM_PROCS"], self.config["FILTER_DETECTABLE"], self.config["FILTER_ZERO"])

    ###########
    ## INPUT ##
    ###########

    def read_file(self, delim = None):
        if delim not in (None, "tab"):
            raise KeyError("delimiter not recognized, please use 'tab' if tab delimited wanted\n default: comma separated")
        elif delim == "tab":
            delim = "\t"
        else:
            delim = ","

        signals = pd.read_csv(self.config["FILE_NAME"], delimiter = delim)

        if self.config["STAGGER"]:
            signals = self._correct_stagger(signals)

        self.config["LEN_SIGNALS"] = len(signals.columns.tolist()) - 1

        signals["detectable"] = signals.iloc[:, 1:(self.config["LEN_SIGNALS"] + 1)].apply(self._detectable, axis = 1)
        signals["zero"] = signals.iloc[:, 1:(self.config["LEN_SIGNALS"] + 1)].apply(self._contains_zero, axis = 1)

        if self.config["FILTER_DETECTABLE"]:
            signals = self.filter_detectable(signals)
            
        if self.config["FILTER_ZERO"]:
            signals = self.filter_zero(signals)

        return signals

    ####################
    ## PRE-PROCESSING ##
    ####################
    
    def filter_zero(self, signals):
        has_zero = signals['zero'] == 1

        signals_cleaned = signals[has_zero]
        return signals_cleaned

    def filter_detectable(self, signals):
        is_detectable = signals['detectable'] == 1

        signals_cleaned = signals[is_detectable]
        return signals_cleaned

    ###########################
    ## COMPUTE RPx FUNCTIONS ## Refactor, implement with pandas, and remove gene specific language
    ###########################

    def compute_phase(self, signal):
        signal = signal.tolist()
        
        components = np.fft.fft(signal)
        phi = np.angle(complex(components[2].real, components[2].imag))
        phase = -12 / np.pi*phi
        if phase < 0: phase = 24 + phase

        return phase

    def compute_rpx(self, signal):
        if not isinstance(signal, list):
            signal = signal.tolist()
        
        DFT_mod = (abs(np.fft.fft(signal)) / len(signal)) ** 2
        # first step to getting the Power Spectral Density (PSD), 
        # the DFT_mod is just a modified Distcrete Fourier Transfor at this stage
    
        PSD_mod = [0] * len(DFT_mod)
        for i in range(1, int(len(DFT_mod) / 2) + 1):
            if i % int(len(DFT_mod) / 2) != 0:
                PSD_mod[i] = 2 * DFT_mod[i]
            else:
                PSD_mod[i] = DFT_mod[i]
                
        RPx = np.log2((PSD_mod[self.config["PERIOD"] * self.config["N_CYCLES"]] + self.config["EPS"]) / (sum(PSD_mod[ 0:int(len(PSD_mod) / 2 + 1) ]) - PSD_mod[self.config["PERIOD"] * self.config["N_CYCLES"]] + self.config["EPS"]))

        return RPx

    def compute_p_value(self, signal_rpx):
        # compute the p-value for a single gene
        signal = signal_rpx.tolist()
        RPx = signal.pop(len(signal_rpx) - 1)
        
        if RPx >= self.config["MIN_RP24"]:
            # first compute random shuffles of expression
            # RPxs_shuffled = []
            queue = multiprocessing.Queue()
            jobs = []
            for i in range(self.config["NUM_PROCS"]):
                process = multiprocessing.Process(target=self._compute_permutations, args=(signal, queue,))
                jobs.append(process)
                process.start()

            for j in jobs:
                j.join()
                
            RPxs_shuffled = []
            while not queue.empty():
                RPxs_shuffled.append(queue.get())

            # next compute the p-value
            if self.config["USE_Z_SCORE"]:
                # Calculate p-value using a zscore and normal distribution
                z_score = (RPx - np.mean(RPxs_shuffled)) / np.std(RPxs_shuffled)
                return stats.norm.sf(z_score)
            else:
                #Calculate empirical p-value
                count = 0.0
                for value in RPxs_shuffled:
                    if value >= RPx:
                        count += 1.0
                return float(count) / len(RPxs_shuffled)
        else:
            # if RP24 less than min defined value (default is 0), then don't test RP24 significance
            return "NOTEST"

    def compute_q_values(self, p_values):
        p_values = p_values.tolist()
        q_values_final = [0] * len(p_values)
        p_values_unsorted = []
        for i, p_value in enumerate(p_values):
            if p_value != "NOTEST":
                p_values_unsorted.append((i, p_value))

        # Sort list by p value
        p_values_sorted = sorted(p_values_unsorted, key=lambda p:p[1])
    
        # calculate q values from sorted list, going backwards so q values can be adjusted if needed
        q_prev = 1.0
        q_values = []
        for i in reversed(range(len(p_values_sorted))):
            j, p_value = p_values_sorted[i]
            r = float(i+1.0)
            q_temp = p_value * len(p_values_sorted) / r
            q_value = min(q_temp, q_prev)
            q_prev = q_value
            q_values.append((j, q_value))

        for i in range(len(q_values)):
            j, q_value = q_values[i]
            q_values_final[j] = q_value

        for i in range(len(p_values)):
            if p_values[i] == "NOTEST":
                q_values_final[i] = "NOTEST"

        return q_values_final
        
    ###############
    ## DETECTION ## Implement this with pandas
    ###############

    def detect(self, signals):
        signals["RPx"] = signals.iloc[:, 1:(self.config["LEN_SIGNALS"] + 1)].apply(self.compute_rpx, axis = 1)

        indexes = list(range(1, self.config["LEN_SIGNALS"] + 1)) + [len(signals.columns.tolist()) - 1] # separated for clarity - just the signal indexes and the index of the RPx
        signals["p_values"] = signals.iloc[:, indexes].apply(self.compute_p_value, axis = 1)

        signals["q_values"] = self.compute_q_values(signals["p_values"])

        signals["phase"] = signals.iloc[:, 1:(self.config["LEN_SIGNALS"] + 1)].apply(self.compute_phase, axis = 1)
        signals["mean"] = signals.iloc[:, 1:(self.config["LEN_SIGNALS"] + 1)].apply(np.mean, axis = 1)
        return signals

    ###########
    ## UTILS ##
    ###########

    def _correct_stagger(self, df):
        cols = df.columns.tolist()
        new_order = []
        new_order.append(cols.pop(0))
        for i in range(self.config["N_CYCLES"]):
            new_order = new_order + cols[i::self.config["N_CYCLES"]]
        return df[new_order]

    def _detectable(self, signal):
        signal = signal.tolist()

        if np.median(signal) >= 1.0 and (max(signal) + self.config["EPS"]) / (min(signal) + self.config["EPS"]) > 1.5:
            return 1
        else:
            return 0

    def _contains_zero(self, signal):
        if 0 in signal.tolist():
            return 0
        else:
            return 1

    def _compute_permutations(self, signal, queue):
        for i in range(int(self.config["NUM_PERMUTATIONS"] / self.config["NUM_PROCS"])):
            if self.config["SHUFFLE_WITH_REPLACEMENT"]:
                shuffled = np.random.choice(signal, len(signal))
                queue.put(self.compute_rpx(shuffled))
            else:
                np.random.shuffle(signal)
                queue.put(self.compute_rpx(signal))
