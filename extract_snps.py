#!/usr/bin/env python3

import re
import os
import time
import datetime
import subprocess
import logging
import multiprocessing
from itertools import repeat
import pandas as pd
import numpy as np
from scriptpy import process_args, ensure_requisite_folders


#' checks existence of a file
def check_existence(path):
  if not os.path.isfile(path):
    raise FileNotFoundError("{path} does not exist".format(**locals()))


#' makes headers for the filtered files; these are saved to outfile
def make_headers(in_dose, in_info, out_dose, out_info):
  dose_proc = subprocess.Popen("zcat {in_dose} | head -100 | grep \# > {out_dose}".format(**locals()), shell = True)
  info_proc = subprocess.Popen("zcat {in_info} | head -1 > {out_info}".format(**locals()), shell = True)
  return dose_proc, info_proc


#' makes dose command to extract SNP via tabix; these are saved to outfile.tmp
def make_dose_cmd(snps, in_file, out_file):
  template = "tabix -p vcf -f {in_file} {c}:{p}-{p} >> {out_file}.tmp"
  pos = [x.split(':') for x in snps]
  cmds = [None]*len(pos)
  for i,v in enumerate(pos):
    c,p = v
    cmds[i] = template.format(**locals())
  return cmds


#' makes info command to extract SNP; these are saved to outfile.tmp
def make_info_cmd(snps, in_file, out_file):
  filters = "||".join(['$1 == "{}"'.format(x) for x in snps])
  cmd = "zcat %(in_file)s | awk '(%(filters)s){print}' > %(out_file)s.tmp" % locals()
  return cmd
  

#' extracts SNP from dose/info using tabix; these are saved to outfile.tmp
def extract(snps, in_dose, in_info, out_dose, out_info):  
  info_cmd = make_info_cmd(snps = snps, in_file = in_info, out_file = out_info)
  dose_cmds = make_dose_cmd(snps = snps, in_file = in_dose, out_file = out_dose)
  
  info_extract = subprocess.Popen(info_cmd, shell = True)
  dose_extract = [subprocess.Popen(cmd, shell = True) for cmd in dose_cmds]
  return dose_extract, info_extract
  
  
# implements extraction for 1 chr
def process_1(c, snp_list, vcf_dir, out_dir):

  in_dose = '{vcf_dir}/chr{c}.dose.vcf.gz'.format(**locals())
  in_info = '{vcf_dir}/chr{c}.info.gz'.format(**locals())
  out_dose = '{out_dir}/chr{c}.dose.vcf.filter'.format(**locals())
  out_info = '{out_dir}/chr{c}.info.filter'.format(**locals())
  
  check_existence(snp_list)
  check_existence(in_info)
  check_existence(in_dose)
  
  snp_file = pd.read_csv(snp_list, header = None, sep = "\s+")
  
  # confirm snp_list format: 3 cols, id format 'chr:pos'
  # format of alleles not specified, but if not correct format will result in no results
  assert snp_file.shape[1] >= 3, "The SNP file's first three columns should contain id, effect allele, noneffect allele, etc. Please check again"
  assert all(snp_file.iloc[:,0].str.contains("\d+:\d+")), "Format of the id column is invalid, should have the format 'chr:pos'"

  snps = snp_file.iloc[:, 0].astype(str)
  snps = [x for x in snps if x.startswith(c + ':')]
  
  if len(snps) > 0:
    # extract headers - should be fast
    # saved to outfile
    dose_header, info_header = make_headers(in_dose, in_info, out_dose, out_info)
    info_header.wait()  
    dose_header.wait()
  
    # extract the rest - dose extract fast w/ tabix, info extract may be bottleneck
    # saved to outfile.tmp
    dose_extract, info_extract = extract(snps, in_dose, in_info, out_dose, out_info)
    time.sleep(30)
    info_extract.wait()
    for p in dose_extract:
      p.wait()

    # remove duplicates (based on allele); final result moves outfile.tmp -> outfile
    check_lines = int(subprocess.check_output("wc -l {out_info}.tmp".format(**locals()), shell=True).decode("utf-8").split()[0])
    if check_lines > 0:
      directory = os.path.dirname( os.path.realpath(__file__) )
      rm_dup_cmd = "Rscript {directory}/rm_duplicates.R --prefix {c} --vcf_dose {out_dose}.tmp --vcf_info {out_info}.tmp --snp_file_name {snp_list} --output_dose {out_dose} --output_info {out_info}".format(**locals())
      subprocess.check_call(rm_dup_cmd, shell = True)
  
      # finish up
      trash = subprocess.check_call("rm -f {out_dose}.tmp {out_info}.tmp".format(**locals()), shell = True)


#' runs extraction across chr in parallel
def run(snp_list, vcf_dir, out_dir):
  
  ensure_requisite_folders(out_dir)
  chrs = list(set([re.findall("\d+", x)[0] for x in os.listdir(vcf_dir) if x.startswith("chr")]))
  logging.info("Found {} files to extract SNPs from".format(len(chrs)))
  assert len(chrs) > 0, "Could not find any dose files in {vcf_dir}, all files should start with 'chr'".format(**locals())
  
  logging.info("Initiating")
  with multiprocessing.Pool(processes = len(chrs)) as pool:
    pool.starmap(process_1, zip(chrs, repeat(snp_list), repeat(vcf_dir), repeat(out_dir)))
  
  #subprocess.check_call("gzip {out_dir}/*".format(**locals()), shell = True)


if __name__ == "__main__":

  import argparse
  parser = argparse.ArgumentParser(description = "EXTRACT SNPS FROM VCF FILES")
  parser.add_argument("--snp_list", help = "File [columns id (format: 'chr:pos'), effect allele, noneffect allele, etc] containing the SNPs to be extracted, no header", required = True)
  parser.add_argument("--vcf_dir", help = "Directory of VCF files", required = True)
  parser.add_argument("--out_dir", help = "Directory to save filtered VCF files", required = True)
  parser.add_argument("--log", default = None, help="Name of log config file [%(default)s]")
  args = process_args(parser)

  try:
    run(snp_list = args.snp_list, vcf_dir = args.vcf_dir, out_dir = args.out_dir)
  except:
    logging.exception("SNP extraction failed")
    raise

  logging.info("Complete")
