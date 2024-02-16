### Base Requirements for a Flask application.
from utilities.file_management import *
from utilities.server_functions import *
from utilities.application_functions import *

### Basic startup, sets a secret key and some standard folders
app = Flask(__name__,static_folder='static')
app.config['SECRET_KEY'] = "01189998819991197253"
app.config['UPLOAD_FOLDER'] = "/uploads/"
app.config['DATABASE_FOLDER'] = "/database/"
app.config['TEMPLATES_AUTO_RELOAD'] = True


### Starting Page
@app.route('/',methods=['GET', 'POST'])
def start_page():
    if request.method=="POST":
        jobid = random_job_identifier()
        session["jobid"] = jobid
        CURRENT_JOBS[jobid] = ProcessJobSettings(request)
        return render_template("running.html")
    return render_template("main.html")

@app.route('/console_output')
def console_output():
    """
    Run AlphaFold2 Job, report back as it progresses
    """
    jobid = session["jobid"]
    settings = CURRENT_JOBS[jobid]
    def inner():
        JobName                      = jobid
        Sequence_1                   = settings["SEQUENCE_1"]
        Sequence_2                   = settings["SEQUENCE_2"]
        msa_mode                     = settings["MSA_MODE"]
        pair_mode                    = settings["PAIR_MODE"]
        model_type                   = settings["MODEL_TYPE"]
        num_recycles                 = settings["NUM_RECYCLES"]
        recycle_early_stop_tolerance = settings["RECYCLE_EARLY_STOP_TOLERANCE"]
        relax_max_iterations         = settings["RELAX_MAX_ITERATIONS"]
        pairing_strategy             = settings["PAIRING_STRATEGY"]
        max_msa                      = settings["MAX_MSA"]
        num_seeds                    = settings["NUM_SEEDS"]
        use_dropout                  = settings["USE_DROPOUT"]
        save_all                     = settings["SAVE_ALL"]
        save_recycles                = settings["SAVE_RECYCLES"]
        dpi                          = settings["DPI"]
        display_images               = settings["DISPLAY_IMAGES"]
        
        Sequence_1 = "".join(Sequence_1.split())
        Sequence_2 = "".join(Sequence_2.split())

        # Validate Sequences:
        valid_letters = set(["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"])
        if not all([set(Sequence_1).issubset(valid_letters),set(Sequence_2).issubset(valid_letters)]):
            return f'<p style="color:#03e7d4"><i>ERROR:  Improper Amino Acid sequence given.  Check for bad one-letter codes.</i></p>'
            
        if all([Sequence_1.strip() != "",Sequence_2.strip() != "" ]):
            query_sequence = Sequence_1 + ":" + Sequence_2
        else:
            query_sequence = Sequence_1 + Sequence_2
        if query_sequence.strip() == "":
            return f'<p style="color:#03e7d4"><i>ERROR:  Empty sequence given.</i></p>'
        yield f'<p style="color:#03e7d4"><b>Complete Query Sequence: </b><br>{query_sequence}</p>'

        
        ######################################
        # Generate Fresh Job Name for Folder #
        ######################################
        jobname = iterative_foldername(JobName)
        yield f'<p style="color:#03e7d4">Creating job directory for job ID: {JobName}<br></p>'
        
        ################
        # Save Queries #
        ################
        queries_path = os.path.join(jobname, "queries.csv")
        with open(queries_path, "w") as text:
            text.write(f"id,sequence\n{jobname},{query_sequence}")
        
        custom_template_path = None
        use_templates = True
        # decide which a3m to use
        if "mmseqs2" in msa_mode:
            a3m_file = os.path.join(jobname,f"{jobname}.a3m")
        
        else:
            a3m_file = os.path.join(jobname,f"{jobname}.single_sequence.a3m")
            with open(a3m_file, "w") as text_file:
                text_file.write(">1\n%s" % query_sequence)

        num_recycles = None if num_recycles == "auto" else int(num_recycles)
        recycle_early_stop_tolerance = None if recycle_early_stop_tolerance == "auto" else float(recycle_early_stop_tolerance)
        if max_msa == "auto": 
            max_msa = None
        
        result_dir = jobname
        log_filename = os.path.join(jobname,"log.txt")
        if not os.path.isfile(log_filename) or 'logging_setup' not in globals():
            setup_logging(Path(log_filename))
        yield '<p style="color:#03e7d4">Processing job settings...</p>'
        
        queries, is_complex = get_queries(queries_path)
        model_type = set_model_type(is_complex, model_type)
        
        if "multimer" in model_type and max_msa is not None:
            use_cluster_profile = False
        else:
            use_cluster_profile = True
        yield '<p style="color:#03e7d4">Downloading AlphaFold parameters.  This may take several minutes.</p>'
        download_alphafold_params(model_type, Path("."))
        yield '<p style="color:#03e7d4">AlphaFold parameters downloaded.  Beginning model generation.  This may take some time.</p>'
        results = run(
            queries=queries,
            result_dir=result_dir,
            use_templates=use_templates,
            custom_template_path=custom_template_path,
            num_relax=0,
            msa_mode=msa_mode,
            model_type=model_type,
            num_models=5,
            num_recycles=num_recycles,
            relax_max_iterations=relax_max_iterations,
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
            user_agent="colabfold/google-colab-main",
        )
        yield '<p style="color:#03e7d4">Job complete!  Zipping results.</p>'
        results_zip = f"/app/download/{jobname}.result.zip"
        os.system(f"zip -r /app/download/{jobname}.results.zip {jobname}")
        yield "<p style='color:#03e7d4'><b>Job Complete!</b></p>"
        yield "<p style='color:#03e7d4'><a href=" + f"download/{jobname}.results.zip" + ">Results (zipped)</a></p>"
    return Response(inner(), mimetype='text/html')

@app.route('/download/<path:filename>', methods=['GET', 'POST'])
def download(filename):
    return send_from_directory('download/', filename, as_attachment=True)

@app.route('/results')
def results():
    return render_template("results.html",jobid=session["jobid"])