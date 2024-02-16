#!/usr/bin/env python3

#######################
### Library Imports ###
#######################
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from Bio import BiopythonDeprecationWarning
warnings.simplefilter(action='ignore', category=BiopythonDeprecationWarning)
from pathlib import Path
from colabfold.download import download_alphafold_params
from colabfold.utils import setup_logging
from colabfold.batch import get_queries, run, set_model_type
from colabfold.plot import plot_msa_v2
from colabfold.colabfold import plot_protein
import matplotlib.pyplot as plt

############################
### Function Definitions ###
############################
display_images = False
def check(folder):
    if os.path.exists(folder):
        return False
    else:
        return True

def iterative_foldername(jobname):
    # check if directory with jobname exists
    if not check(jobname):
        n = 0
        while not check(f"{jobname}_{n}"): 
            n += 1
        jobname = f"{jobname}_{n}"
    os.makedirs(jobname, exist_ok=True)
    return jobname

def input_features_callback(input_features):
    if display_images:
        print("Attempting to display images")
        plot_msa_v2(input_features)
        plt.show()
        plt.close()

def prediction_callback(protein_obj, length, prediction_result, input_features, mode):
    model_name, relaxed = mode
    if all([not relaxed, display_images]):
        print("Attempting to display images")
        fig = plot_protein(protein_obj, Ls=length, dpi=150)
        plt.show()
        plt.close()

def ProcessJobSettings(request):
    settings = {}
    settings["SEQUENCE_1"] = request.form.get('SEQUENCE_1').upper()
    settings["SEQUENCE_2"] = request.form.get('SEQUENCE_2').upper()
    settings["MSA_MODE"] = request.form.get('MSA_MODE')
    settings["PAIR_MODE"] = request.form.get('PAIR_MODE')
    settings["MODEL_TYPE"] = request.form.get('MODEL_TYPE')
    settings["NUM_RECYCLES"] = request.form.get('NUM_RECYCLES')
    settings["RECYCLE_EARLY_STOP_TOLERANCE"] = request.form.get('RECYCLE_EARLY_STOP_TOLERANCE')
    settings["RELAX_MAX_ITERATIONS"] = request.form.get('RELAX_MAX_ITERATIONS')
    settings["PAIRING_STRATEGY"] = request.form.get('PAIRING_STRATEGY')
    settings["MAX_MSA"] = request.form.get('MAX_MSA')
    settings["NUM_SEEDS"] = int(request.form.get('NUM_SEEDS'))
    settings["USE_DROPOUT"] = bool(request.form.get('USE_DROPOUT'))
    settings["SAVE_ALL"] = bool(request.form.get('SAVE_ALL'))
    settings["SAVE_RECYCLES"] = bool(request.form.get('SAVE_RECYCLES'))
    settings["DISPLAY_IMAGES"] = bool(request.form.get('DISPLAY_IMAGES'))
    settings["DPI"] = int(request.form.get('DPI'))
    return settings
