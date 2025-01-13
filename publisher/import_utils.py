
import logging
from datetime import datetime
import pandas as pd
import numpy as np
from model_import_actions import model_import_actions
from sys import stderr

from dotenv import load_dotenv
import os
from natsort import natsorted

data_issue_logger = {}
output_logger = None
error_logger = None

load_dotenv()

chunk_size = int(os.environ.get("CHUNK_SIZE"))
fail_fast = os.environ.get("FAIL_FAST") == "true" or os.environ.get("FAIL_FAST") == "True"
verbose = os.environ.get("VERBOSE") == "true" or os.environ.get("VERBOSE") == "True"
print('verbose is', verbose)

def get_job_dir( db_name, assembly, dry_run, tsv_dir, start_at_model, start_at_file):
        
    current_dir = os.path.dirname(os.path.realpath(__file__))
    jobs_dir = os.path.abspath(os.path.join(current_dir, "jobs"))    
    
    date_yyyy_mm_dd = datetime.now().strftime("%Y-%m-%d")
    start_at_info = f"{'-'+start_at_model if start_at_model else '' }{'-'+start_at_file+'-' if start_at_file else ''}"
    desired_job_dir_name = f"{date_yyyy_mm_dd}-{tsv_dir}-v{assembly}-into-{db_name}{start_at_info}{'-dryrun' if dry_run else ''}_1";
    
    if os.path.exists(os.path.join(jobs_dir, desired_job_dir_name)):
        dir_exists = True
        while dir_exists:
            number_at_end = int(desired_job_dir_name.split("_")[-1])
            up_to_underscore = "_".join(desired_job_dir_name.split("_")[:-1])
            desired_job_dir_name = f"{up_to_underscore}_{str(number_at_end+1)}"
            dir_exists = os.path.exists(os.path.join(jobs_dir, desired_job_dir_name))

    return os.path.join(jobs_dir, desired_job_dir_name)
            
        
def inspectTSV(file):
    total_rows = 0
    separator = "\t"
    small_read = pd.read_csv(file, sep="\t", nrows=4, header=None)

    num_columns = small_read.shape[1]
    if num_columns == 1 and "gene" not in file:
        separator = " "
        small_read = pd.read_csv(file, sep=" ", nrows=4, header=None)
        if small_read.shape[1] == 1:
            print("Could not determine separator for " + file)
            quit()

    columns = [col.lower() for col in small_read.values.tolist()[0]]
# Read the first chunk to infer data types 
    first_chunk = pd.read_csv(file, sep=separator, chunksize=chunk_size, header=0) 
    type_dict = dict(first_chunk.get_chunk(1).dtypes) 
    first_chunk.close()

    for chunk in pd.read_csv(file, sep="\t", chunksize=chunk_size):
        total_rows += len(chunk)
        
#    print("types", type_dict)

    return {
        "total_rows": total_rows,
        "num_columns": num_columns,
        "columns": columns,
        "types": type_dict,
        "separator": separator
#        "chromosome"
    }

def readTSV(file, info, dtype={}):
#    df = pd.read_csv(file, sep=info["separator"], dtype=dtype, na_values=["NA"], keep_default_na=False)
    
    df = pd.read_csv(file, sep=info["separator"], dtype=dtype)

    df.rename(columns={"All_info$variant": "variant"}, inplace=True)
    df.columns = [col.lower() for col in df.columns]
    
    df.replace(np.nan, None, inplace=True)
    df.replace(".", None, inplace=True)
    
    return df

def setup_loggers(job_dir):
    global data_issue_logger,output_logger, error_logger
    model_names = list(model_import_actions.keys())
    model_names.append("severities")
    for model_name in model_names:
        a_logger = logging.getLogger(model_name)
        a_logger.setLevel(logging.WARNING)
        a_logger.propagate = False        
        a_logger_handler = logging.FileHandler(os.path.join(job_dir,f"warnings-{model_name}.log"))
        a_logger_handler.setLevel(logging.WARNING)
        a_logger.addHandler(a_logger_handler)
        data_issue_logger[model_name] = a_logger
        
    output_logger = logging.getLogger("output")
    output_logger.setLevel(logging.INFO)
    output_logger.propagate = False

    output_logger_handler = logging.FileHandler(os.path.join("./",job_dir,"output.log"))
    output_logger_handler.setLevel(logging.INFO)
    output_logger.addHandler(output_logger_handler)

    output_logger.info(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    
    error_logger = logging.getLogger("error")
    error_logger.setLevel(logging.ERROR)
    error_logger.propagate = False
    error_logger_handler = logging.FileHandler(os.path.join(job_dir,"error.log"))
    error_logger_handler.setLevel(logging.ERROR)
    error_logger.addHandler(error_logger_handler)
    
    print("logging to", os.path.join(job_dir,"error.log"), os.path.join(job_dir,"output.log"))

def log_data_issue(s, model=None):
    if model is not None:
        data_issue_logger[model].warning(s)
    if (fail_fast):
        print(s)
        exit()
    if (verbose):
        print(s)
def log_output(s):
    output_logger.info(s)
    if (verbose):
        print(s)
def log_error(s):
    error_logger.error(s)
    if (fail_fast):
        print(s)
        exit()
    if (verbose):
        print(s)

def report_counts(counts):
    percent_success = "N/A"
    if counts["success"] + counts["fail"] > 0:
        percent_success = (
            100
            * (counts["success"])
            / (counts["rowcount"])
        )
    log_output(
        f"{percent_success}% published. "+
        f"success={counts['success']} "+
        f"fail={counts['fail']} "+
        f"rowcount={counts['rowcount']} "+
        f"missingrefs={counts['missingRef']} "+
        f"inserts={counts['inserted']} "+
        f"duplicates={counts['duplicate']} "+
        f"updates={counts['updated']} "+
        f"successfulchunks={counts['successful_chunks']} "+
        f"failchunks={counts['fail_chunks']}"  
    )


def chunks(l, n):
    """Yield successive n-sized chunks from list l."""
    for i in range(0, len(l), n):
        yield l[i : i + n]