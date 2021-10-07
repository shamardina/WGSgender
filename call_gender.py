from __future__ import print_function
import subprocess
import datetime
import sys
import os
import csv
import re
import numpy as np
import multiprocessing
import shutil
import collections
import shlex
import argparse

Dataset = collections.namedtuple("Dataset", ["filename", "mtime", "header", "data", "updated"])


def read_conf(filename):
  """Reading and setting up configuration parameters"""
  conf = {}
  with open(filename, "r") as f:
    for l in f:
      if "=" in l:
        line = l.rstrip().split("=")
        k = line[0].strip()
        v = re.sub("#.+$", "", line[1].strip()).strip()
        conf[k] = v
  conf["H_threshold"] = float(conf["H_threshold"])
  conf["nonPAR_start"] = int(conf["nonPAR_start"])
  conf["nonPAR_end"] = int(conf["nonPAR_end"])
  conf["NoSDs"] = int(conf["NoSDs"])
  conf["cores"] = int(conf["cores"])
  conf["autosomes"] = map(lambda x: conf["chrom_prefix"] + str(x), range(1, 23))
  conf["chromosomes"] = conf["autosomes"] + [conf["chrom_prefix"] + "X", conf["chrom_prefix"] + "Y"]
  conf["HHH_file"] = os.path.join(conf["output"], conf["name"], "hethomhemi.csv")
  conf["RC_file"] = os.path.join(conf["output"], conf["name"], "read_counts.csv")
  conf["chrom_lengths"] = {}
  with open(conf["chrom_lengths_file"], "r") as f:
    for l in f:
      line = l.rstrip().split("\t")
      if not line[0].startswith(conf["chrom_prefix"]):  # be kind and add "chr" if needed
        line[0] = conf["chrom_prefix"] + line[0]
      conf["chrom_lengths"][line[0]] = int(line[1])
  return conf


def read_manifest(filename):
  """Reading the manifest file. Checking uniqueness of samples (warning only if non-unique).
  Checking 'Included' is 0 or 1, 'Declared_gender' is F, M or NA (exiting if not)"""
  manifest = {}
  with open(filename, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
      if row["ID"] in manifest:
        print("WARNING: sample {} is duplicated in the manifest file {}. Will use the last entry and continue.".format(
            row["ID"], filename))
      if row["Included"] not in {"0", "1"}:
        sys.exit("ERROR: 'Included' should be 0 or 1 in the manifest, sample {} has '{}' in {}, exiting".format(
            row["ID"], row["Included"], filename))
      if row["Declared_gender"] not in {"F", "M", "NA"}:
        sys.exit("ERROR: 'Declared_gender' should be F, M or NA in the manifest, sample {} has '{}' in {}, exiting".format(
            row["ID"], row["Declared_gender"], filename))
      manifest[row["ID"]] = row
  return manifest


def read_data(filename):
  """Reading the saved data from previous runs"""
  data = {}
  with open(filename, "r") as f:
    header = f.readline().rstrip().split("\t")
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
      data[row[0]] = row
  return Dataset(
      filename=filename,
      mtime=datetime.datetime.fromtimestamp(os.path.getmtime(filename)),
      header=header,
      data=data,
      updated=[False])


def write_new(newdataset):
  """Rewriting the results if anything is updated. Previous results saved with .old"""
  if not newdataset.updated[0]:
    print(newdataset.filename + " is up to date")
    return
  shutil.move(newdataset.filename, newdataset.filename + ".old")
  with open(newdataset.filename, "w") as f:
    writer = csv.writer(f, delimiter="\t", lineterminator="\n")
    writer.writerow(newdataset.header)
    for s in sorted(newdataset.data.iterkeys()):
      writer.writerow(newdataset.data[s])


def init(conf):
  """Initialising directory structure and output files;
  checking if samtools and bcftools are available"""
  if not os.path.exists(conf["output"]):
    os.mkdir(conf["output"], 0775)
  output = os.path.join(conf["output"], conf["name"])
  if not os.path.exists(output):
    os.mkdir(output, 0775)
  FNULL = open(os.devnull, 'w')
  print("Checking samtools is available: ", end="")
  if subprocess.call(["which", "samtools"], stdout=FNULL, stderr=subprocess.STDOUT) != 0:
    sys.exit("\nERROR: samtools is required for this script, exiting")
  print("Done")
  print("Checking bcftools is available: ", end="")
  if subprocess.call(["which", "bcftools"], stdout=FNULL, stderr=subprocess.STDOUT) != 0:
    sys.exit("\nERROR: bcftools is required for this script, exiting")
  print("Done")
  FNULL.close()
  if not os.path.exists(conf["HHH_file"]):
    with open(conf["HHH_file"], "w") as f:
      f.write("\t".join(["ID", "filename", "het", "hom", "hemi"]) + "\n")
  if not os.path.exists(conf["RC_file"]):
    with open(conf["RC_file"], "w") as f:
      f.write("\t".join(["ID", "filename"] + conf["chromosomes"]) + "\n")


def count_reads(filename, conf):
  """Getting read counts from BAM file with samtools idxstats"""
  print(filename)
  RC = {}
  output = subprocess.check_output(["samtools", "idxstats", filename]).strip()
  lines = output.split("\n")
  for line in lines:
    a = line.split("\t")
    RC[a[0]] = a[2]
  return [RC[k] for k in conf["chromosomes"]]


def count_hethomhemi(filename, conf):
  """Getting numbers of heterozygous, homozygous and hemizygous passing SNPs
  in the non-PAR of X from VCF file with bcftools query"""
  print(filename)
  command = "bcftools query -i 'FILTER==\"PASS\" && TYPE==\"snp\"' -r {}:{}-{} -f '[%GT]\n' {}".format(
      conf["chrom_prefix"] + "X", conf["nonPAR_start"], conf["nonPAR_end"], filename)
  output = subprocess.check_output(shlex.split(command)).strip()
  counts = {"0/1": 0, "1/1": 0, "1": 0}
  for gt in output.split("\n"):
    if gt not in counts:
      counts[gt] = 0
    counts[gt] += 1
  return [counts[i] for i in ["0/1", "1/1", "1"]]


def count_star(args):
  """Wrapper for calling count functions in multiprocessing pool.map"""
  func = args[0]
  a = args[1]
  return func(*a)


def needs_updating(sample, filename, current_data, date):
  """Checking if sample's info needs updating. Update if:
     - if it is a new sample;
     - filename has changed for old sample;
     - file modification time has changed since last results were generated"""
  if sample not in current_data:
    return True
  if filename != current_data[sample][1]:
    return True
  if datetime.datetime.fromtimestamp(os.path.getmtime(filename)) > date:
    return True
  return False


def drop_removed(dataset, manifest):
  """Remove samples from results if they are removed form the manifest"""
  to_drop = []
  for sample in dataset.data:
    if sample not in manifest:
      to_drop.append(sample)
  if len(to_drop) > 0:
    for sample in to_drop:
      print("Removed outdated sample: " + sample)
      del dataset.data[sample]
    dataset.updated[:] = [True]


def update(dataset, manifest, filetype, count_func, conf):
  """Updating previous results"""
  to_update = []
  for sample in manifest:
    if needs_updating(sample, manifest[sample][filetype], dataset.data, dataset.mtime):
      to_update.append(sample)
  if len(to_update) > 0:
    pool = multiprocessing.Pool(processes=conf["cores"])
    tasks = [(count_func, (manifest[i][filetype], conf)) for i in to_update]
    updated_results = pool.map(count_star, tasks)
    for i, sample in enumerate(to_update):
      dataset.data[sample] = [sample, manifest[sample][filetype]] + updated_results[i]
    dataset.updated[0] = True
  drop_removed(dataset, manifest)


def find_H_ratio(dataset):
  """Calculating H_ratio = #het/(#hom + #hemi)"""
  H_ratio = {}
  het_field, hom_field, hemi_field = [dataset.header.index(i) for i in ["het", "hom", "hemi"]]
  for sample in dataset.data:
    if dataset.data[sample][hom_field] == 0 and dataset.data[sample][hemi_field] == 0:
      H_ratio[sample] = "NA"
    else:
      H_ratio[sample] = float(
          dataset.data[sample][het_field]) / (float(dataset.data[sample][hom_field]) + float(dataset.data[sample][hemi_field]))
  return H_ratio


def find_cov_ratios(dataset, conf):
  """Calculating ratios of normalised X and Y coverage to median of the normalised autosomal coverage"""
  cov_ratios = {}
  for sample in dataset.data:
    RC_values = {
        chrom: float(dataset.data[sample][dataset.header.index(chrom)]) / conf["chrom_lengths"][chrom]
        for chrom in conf["chromosomes"]
    }
    auto_median = np.median([RC_values[chrom] for chrom in conf["autosomes"]])
    if auto_median == 0:
      X_auto = "NA"
      Y_auto = "NA"
    else:
      X_auto = RC_values[conf["chrom_prefix"] + "X"] / auto_median
      Y_auto = RC_values[conf["chrom_prefix"] + "Y"] / auto_median
    cov_ratios[sample] = {"XAutoRatio": X_auto, "YAutoRatio": Y_auto}
  return cov_ratios


def make_call(rc, hhh, manifest, conf):
  """Getting the input and calling set_gender_and_flags()"""
  rv = {}
  H_ratio = find_H_ratio(hhh)
  cov_ratios = find_cov_ratios(rc, conf)
  males = []
  females = []
  for sample, H in H_ratio.iteritems():
    if manifest[sample]["Included"] != "1":  # use only included samples for statistical estimations
      continue
    if H < conf["H_threshold"]:
      males.append(sample)
    else:
      females.append(sample)
  mean_male_X = np.mean([cov_ratios[sample]["XAutoRatio"] for sample in males])
  mean_male_Y = np.mean([cov_ratios[sample]["YAutoRatio"] for sample in males])
  mean_female_X = np.mean([cov_ratios[sample]["XAutoRatio"] for sample in females])
  mean_female_Y = np.mean([cov_ratios[sample]["YAutoRatio"] for sample in females])
  dev_male_X = conf["NoSDs"] * np.std([cov_ratios[sample]["XAutoRatio"] for sample in males], ddof=1)
  dev_male_Y = conf["NoSDs"] * np.std([cov_ratios[sample]["YAutoRatio"] for sample in males], ddof=1)
  dev_female_X = conf["NoSDs"] * np.std([cov_ratios[sample]["XAutoRatio"] for sample in females], ddof=1)
  dev_female_Y = conf["NoSDs"] * np.std([cov_ratios[sample]["YAutoRatio"] for sample in females], ddof=1)
  means = (mean_male_X, mean_male_Y, mean_female_X, mean_female_Y, dev_male_X, dev_male_Y, dev_female_X, dev_female_Y)
  for sample in manifest:
    rv[sample] = set_gender_and_flags(cov_ratios[sample]["XAutoRatio"], cov_ratios[sample]["YAutoRatio"], H_ratio[sample],
                                      manifest[sample]["Declared_gender"], means, conf["H_threshold"])
    rv[sample]["ID"] = sample
  return rv


def in_gates(X, Y, Xmin, Xmax, Ymin, Ymax):
  """Checking if a point is located inside given boundaries"""
  if (Xmin < X < Xmax) and (Ymin < Y < Ymax):
    return True
  return False


def set_gender_and_flags(X, Y, H, declared_gender, means, H_threshold):
  """Inferring gender, sex karyotype; setting flags for outliers"""
  mean_male_X, mean_male_Y, mean_female_X, mean_female_Y, dev_male_X, dev_male_Y, dev_female_X, dev_female_Y = means
  rv = {"XAutoRatio": X, "YAutoRatio": Y, "Hratio": H, "Declared_gender": declared_gender}
  if in_gates(X, Y, mean_male_X - dev_male_X, mean_male_X + dev_male_X, mean_male_Y - dev_male_Y,
              mean_male_Y + dev_male_Y):  # normal XY
    rv.update({"Genotype_gender": "M", "Karyotype": "XY", "Discordant_Hratio": int(H != 0), "Outside_thresholds": 0})
  elif in_gates(X, Y, mean_female_X - dev_female_X, mean_female_X + dev_female_X, mean_female_Y - dev_female_Y,
                mean_female_Y + dev_female_Y):  # normal XX
    rv.update({"Genotype_gender": "F", "Karyotype": "XX", "Discordant_Hratio": int(H < H_threshold), "Outside_thresholds": 0})
  elif in_gates(X, Y, mean_female_X - 0.5 * dev_female_X, mean_female_X + 0.5 * dev_female_X, mean_male_Y - 0.5 * dev_male_Y,
                mean_male_Y + 0.5 * dev_male_Y):  # inferred XXY
    rv.update({"Genotype_gender": "M", "Karyotype": "XXY", "Discordant_Hratio": 0, "Outside_thresholds": 1})
  elif in_gates(X, Y, 1.5 * mean_female_X - 0.5 * dev_female_X, 1.5 * mean_female_X + 0.5 * dev_female_X,
                mean_female_Y - dev_female_Y, mean_female_Y + dev_female_Y):  # inferred XXX
    rv.update({"Genotype_gender": "F", "Karyotype": "XXX", "Discordant_Hratio": 0, "Outside_thresholds": 1})
  elif mean_male_X - dev_male_X < X < mean_male_X + dev_male_X and Y >= mean_male_Y + dev_male_Y:  # inferred XYY
    rv.update({"Genotype_gender": "M", "Karyotype": "XYY", "Discordant_Hratio": int(H != 0), "Outside_thresholds": 1})
  elif H < H_threshold and Y <= mean_male_Y - dev_male_Y and mean_male_X - dev_male_X < X < mean_male_X + dev_male_X:  # inferred X0
    rv.update({"Genotype_gender": "F", "Karyotype": "X0", "Discordant_Hratio": int(H != 0), "Outside_thresholds": 1})
  elif H >= H_threshold and in_gates(
      X, Y, mean_female_X / 2.0, mean_female_X - 0.9 * dev_female_X, mean_female_Y - 0.5 * dev_female_Y,
      mean_female_Y + 0.5 * dev_female_Y
  ):  # probably partial X0 or XX with large deletions; greedy (and safe here) approach to include X right hand side boundary
    rv.update({"Genotype_gender": "F", "Karyotype": "XX", "Discordant_Hratio": 0, "Outside_thresholds": 1})
  else:
    rv.update({"Genotype_gender": "NA", "Karyotype": "NA", "Discordant_Hratio": 0, "Outside_thresholds": 1})
  rv["Gender_mismatch"] = int(declared_gender != "NA" and declared_gender != rv["Genotype_gender"])
  if rv["Discordant_Hratio"] == 1 or rv["Outside_thresholds"] == 1 or rv["Gender_mismatch"] == 1:
    rv["Flag"] = 1
  else:
    rv["Flag"] = 0
  return rv


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument(
      "-c",
      "--conf",
      help="provide path to the configuration file [default.conf]",
      default="default.conf",
      metavar="FILENAME",
      dest="conf_file")
  parser.add_argument(
      "-p",
      "--process-only",
      help="skip the counting step and use the output of the previous runs instead",
      action="store_true",
      dest="process_only")
  args = parser.parse_args()
  CONF = read_conf(args.conf_file)
  if not args.process_only:
    init(CONF)
  MANIFEST = read_manifest(CONF["manifest"])
  HHH = read_data(CONF["HHH_file"])
  RC = read_data(CONF["RC_file"])
  if not args.process_only:
    update(HHH, MANIFEST, "VCF", count_hethomhemi, CONF)
    write_new(HHH)
    update(RC, MANIFEST, "BAM", count_reads, CONF)
    write_new(RC)
  result = make_call(RC, HHH, MANIFEST, CONF)
  header = [
      "ID", "XAutoRatio", "YAutoRatio", "Hratio", "Declared_gender", "Genotype_gender", "Gender_mismatch", "Discordant_Hratio",
      "Outside_thresholds", "Flag", "Karyotype"
  ]
  with open(os.path.join(CONF["output"], CONF["name"], "gender.csv"), "w") as f:
    f.write("\t".join(header) + "\n")
    for sample in sorted(result.keys()):
      f.write("\t".join([str(result[sample][k]) for k in header]) + "\n")


if __name__ == "__main__":
  main()
