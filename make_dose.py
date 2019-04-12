#!/usr/bin/env python3

import os
import time
import datetime
import subprocess
import logging
import pandas as pd
import numpy as np
from scriptpy import process_args, assert_dim

# combine VCF, generate dosage file, add the SNP column to dosage file
def run(folder, out):

  logging.info("Combining individual VCF files")
  this_script = os.path.join(os.path.dirname(os.path.realpath(__file__)), "combine_vcf.py")
  subprocess.check_call("python3 %(this_script)s --folder %(folder)s" % locals(), shell = True)

  logging.info("Generating dosage file")
  subprocess.check_call("vcftools --vcf combined.dose.txt --extract-FORMAT-info DS --out combined", shell = True)

  logging.info("Making final dosage file")
  info = pd.read_csv("combined.info.txt", sep = "\t", usecols = [0,1,2]) # cols [0,2) refers to SNP, REF, ALT (info file)
  info.columns = ['SNP', 'REF0', 'ALT1']

  # maintain consistency between combine_vcf.py and vcftools
  nr = int(subprocess.check_output("grep '^#' combined.dose.txt | sed -n '$='", shell = True).decode('utf-8').split('\n')[0]) - 1 # remove all comment rows except header
  dose_other = pd.read_csv("combined.dose.txt", sep = "\t", usecols = [3,4], skiprows = range(nr)) # cols [3,4] refers to REF, ALT (dose file)
  dose = pd.read_csv("combined.DS.FORMAT", sep = "\t")

  dose['SNP'] = ['%s:%s' % x for x in zip(dose.CHROM, dose.POS)]
  dose['REF0'] = dose_other['REF']
  dose['ALT1'] = dose_other['ALT']

  # ensures that the order of info and new dose are the same (by placing info first in merge)
  m = pd.merge(info, dose, on = ['SNP', 'REF0', 'ALT1'])

  @assert_dim(dose.shape)
  def check_shape(x):
    return x
  check_shape(m)

  logging.info("Writing out {}".format(out))
  m.to_csv(out, sep = "\t", index = False)


# process args and run program
if __name__ == "__main__":

  import argparse
  parser = argparse.ArgumentParser(description = "COMBINES VCF FILES AND MAKES DOSAGES")
  parser.add_argument("--folder", help = "Folder that contains all files", required = True)
  parser.add_argument("--out", default = "final_dose.txt", help = "Name of output file [%(default)s]")
  parser.add_argument("--log", default = None, help="Prefix of log file [%(default)s]")
  args = process_args(parser)

  try:
    run(folder = args.folder, out = args.out)
  except Exception as e:
    logging.exception("Making dosages failed")
    raise

  logging.info("Complete")

