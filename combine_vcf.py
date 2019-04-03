#!/usr/bin/env python3

import os
import time
import datetime
import subprocess
import logging
import multiprocessing
import itertools
import pandas as pd
import numpy as np
from scriptpy import process_args, report_function_name

# get files
@report_function_name
def get_files(folder):
  logging.info('Combining files in {folder}'.format(**locals()))
  vcfs = [os.path.join(folder, x) for x in os.listdir(folder) if 'dose' in x]
  infos = [os.path.join(folder, x) for x in os.listdir(folder) if 'info' in x]
  return vcfs, infos

# combine dose via vcf-concat
@report_function_name
def start_dose(vcfs):
  logging.info("Combining Dose Files")
  proc_dose = subprocess.Popen('vcf-concat {} > combined.dose.txt'.format(' '.join(vcfs)), shell = True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  return proc_dose  


# dose is the bottle neck (bigger file)
@report_function_name
def wait_for_dose(proc_dose):
  proc_dose.wait()
  out, err = proc_dose.communicate()

  if len(out) > 0:
    logging.info(out.decode('utf-8'))
    
  if len(err) > 0:
    logging.error(err.decode('utf-8'))

    
# combine info manually
@report_function_name
def run_info(infos):
  logging.info("Combining the Info Files")
  subprocess.check_call('zcat {} > combined.info.txt'.format(infos[0]), shell = True)
  for s in infos[1:]:
    subprocess.check_call('zcat {} | awk "NR>1" >> combined.info.txt'.format(s), shell = True)


# QC: check that rows are the same
@report_function_name
def run_qc():
  logging.info("Quality Check")
  a = int(subprocess.check_output("sed -n '$=' combined.info.txt", shell = True).decode('utf-8').split('\n')[0])
  b = int(subprocess.check_output("grep -v '^#' combined.dose.txt | sed -n '$='", shell = True).decode('utf-8').split('\n')[0]) + 1 
  assert a == b, "Mismatched dose and info files"

  logging.info("SNPs found (dose):{a}".format(**locals()))
  logging.info("SNPs found (info):{b}".format(**locals()))


def run(folder):
  vcfs, info = get_files(folder)
  
  proc_dose = start_dose(vcfs)
  run_info(info)
  wait_for_dose(proc_dose)
  
  run_qc()
  

if __name__ == "__main__":

  import argparse
  parser = argparse.ArgumentParser(description = "COMBINE VCF FILES")
  parser.add_argument("--folder", help="Folder location of vcf data", required = True)
  parser.add_argument("--log", default = None, help="Prefix of log file [%(default)s]")
  args = process_args(parser)

  try:
    run(args.folder)
  except Exception as e:
    logging.exception("Combine VCF failed")
    raise

  logging.info("Complete")
