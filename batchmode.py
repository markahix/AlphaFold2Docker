#!/usr/bin/env python3
#######################
### LIBRARY IMPORTS ###
#######################
# No point in continuing if the libraries aren't all available!
# from google.colab import files
import os
import numpy as np
import re
import hashlib
import random
from pathlib import Path
import matplotlib.pyplot as plt
import sys
import glob
import warnings
from sys import version_info
python_version = f"{version_info.major}.{version_info.minor}"
warnings.simplefilter(action='ignore', category=FutureWarning)

from Bio import BiopythonDeprecationWarning
warnings.simplefilter(action='ignore', category=BiopythonDeprecationWarning)

from colabfold.download import download_alphafold_params, default_data_dir
from colabfold.utils import setup_logging
from colabfold.batch import get_queries, run, set_model_type
from colabfold.plot import plot_msa_v2
from colabfold.colabfold import plot_protein

def sequence_to_structure_alphafold(query_sequence="",
                                    jobname="",
                                    num_relax=0,
                                    template_mode=None,
                                    msa_mode = "mmseqs2_uniref_env",
                                    pair_mode = "unpaired_paired",
                                    model_type = "auto",
                                    num_recycles = "auto",
                                    recycle_early_stop_tolerance = "auto",
                                    pairing_strategy = "greedy",
                                    max_msa = "auto",
                                    num_seeds = 1,
                                    save_all = False,
                                    save_recycles = False,
                                    use_dropout = False,
                                    dpi = 300,
                                    display_images = True):
    num_recycles = None if num_recycles == "auto" else int(num_recycles)
    recycle_early_stop_tolerance = None if recycle_early_stop_tolerance == "auto" else float(recycle_early_stop_tolerance)
    if max_msa == "auto": 
        max_msa = None
    ############################
    ### FUNCTION DEFINITIONS ###
    ############################
    def add_hash(x,y):
        return x+"_"+hashlib.sha1(y.encode()).hexdigest()[:5]

    # check if directory with jobname exists
    def check(folder):
        if os.path.exists(folder):
            return False
        else:
            return True
        
    def input_features_callback(input_features):
        if display_images:
            plot_msa_v2(input_features)
            plt.show()
            plt.close()

    def prediction_callback(protein_obj, length,
                            prediction_result, input_features, mode):
        model_name, relaxed = mode
        if not relaxed:
            if display_images:
                fig = plot_protein(protein_obj, Ls=length, dpi=150)
                plt.show()
                plt.close()

    # remove whitespaces
    query_sequence = "".join(query_sequence.split())
    basejobname = "".join(jobname.split())
    basejobname = re.sub(r'\W+', '', basejobname)
    jobname = add_hash(basejobname, query_sequence)

    if not check(jobname):
        n = 0
        while not check(f"{jobname}_{n}"): n += 1
        jobname = f"{jobname}_{n}"

    # make directory to save results
    os.makedirs(jobname, exist_ok=True)

    # save queries
    queries_path = os.path.join(jobname, f"{jobname}.csv")
    with open(queries_path, "w") as text_file:
        text_file.write(f"id,sequence\n{jobname},{query_sequence}")

    if template_mode == "pdb100":
        use_templates = True
        custom_template_path = None
    else:
        custom_template_path = None
        use_templates = False

    # decide which a3m to use
    if "mmseqs2" in msa_mode:
        a3m_file = os.path.join(jobname,f"{jobname}.a3m")

    else:
        a3m_file = os.path.join(jobname,f"{jobname}.single_sequence.a3m")
        with open(a3m_file, "w") as text_file:
            text_file.write(">1\n%s" % query_sequence)

    result_dir = jobname
    log_filename = os.path.join(jobname,"log.txt")
    if not os.path.isfile(log_filename) or 'logging_setup' not in globals():
        setup_logging(Path(log_filename))
        logging_setup = True

    queries, is_complex = get_queries(queries_path)
    model_type = set_model_type(is_complex, model_type)

    if "multimer" in model_type and max_msa is not None:
        use_cluster_profile = False
    else:
        use_cluster_profile = True

    download_alphafold_params(model_type, Path("."))  # Can I do this for each model type option above and then just have it done at the start of the service?
    results = run(
        queries=queries,
        result_dir=result_dir,
        use_templates=use_templates,
        custom_template_path=custom_template_path,
        num_relax=num_relax,
        msa_mode=msa_mode,
        model_type=model_type,
        num_models=5,
        num_recycles=num_recycles,
        recycle_early_stop_tolerance=recycle_early_stop_tolerance,
        num_seeds=num_seeds,
        use_dropout=use_dropout,
        model_order=[1,2,3,4,5],
        is_complex=is_complex,
        data_dir=Path("."),
        keep_existing_results=False,
        rank_by="auto",
        pair_mode=pair_mode,
        pairing_strategy=pairing_strategy,
        stop_at_score=float(100),
        prediction_callback=prediction_callback,
        dpi=dpi,
        zip_results=False,
        save_all=save_all,
        max_msa=max_msa,
        use_cluster_profile=use_cluster_profile,
        input_features_callback=input_features_callback,
        save_recycles=save_recycles,
    )
    results_zip = f"{jobname}.result.zip"
    os.system(f"zip -r {results_zip} {jobname}")
    if glob.glob(results_zip):
        os.system(f"rm -r {jobname}")