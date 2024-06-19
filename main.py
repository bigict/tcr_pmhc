"""Tools for preprocess TCR_pMHC data"""
import os
import re
from collections import defaultdict
import csv
from dataclasses import dataclass
import functools
import multiprocessing as mp
import random
from urllib.parse import urlparse, parse_qs

from profold2.data.parsers import parse_fasta

_seq_index_pattern = "(\\d+)-(\\d+)"


def decompose_pid(pid, return_domain=False):
  k = pid.find("/")
  if k != -1:
    pid, domains = pid[:k], pid[k + 1:]
  else:
    domains = None

  k = pid.rfind("_")
  if k != -1:
    pid, chain = pid[:k], pid[k + 1:]
  else:
    chain = None

  if return_domain:
    return pid, chain, domains
  return pid, chain


def seq_index_split(text):
  for s in text.split(","):
    r = re.match(_seq_index_pattern, s)
    assert r
    yield tuple(map(int, r.group(1, 2)))


def seq_index_join(seq_index):
  return ",".join(f"{i}-{j}" for i, j in seq_index)


@dataclass
class DBUri:
  path: str
  chain_idx: str = "chain.idx"
  mapping_idx: str = "mapping.idx"
  attr_idx: str = "attr.idx"


def parse_db_uri(db_uri):
  o = urlparse(db_uri)
  chain_idx, mapping_idx, attr_idx = "chain.idx", "mapping.idx", "attr.idx"
  if o.query:
    attrs = parse_qs(o.query)
    if "chain_idx" in attrs:
      chain_idx = attrs["chain_idx"][0]
    if "mapping_idx" in attrs:
      mapping_idx = attrs["mapping_idx"][0]
    if "attr_idx" in attrs:
      chain_idx = attrs["attr_idx"][0]
  return DBUri(o.path, chain_idx, mapping_idx, attr_idx=attr_idx)


def read_mapping_idx(data_uri, mapping_dict=None):
  if mapping_dict is None:
    mapping_dict = {}
  mapping_idx_path = os.path.join(data_uri.path, "mapping.idx")
  if os.path.exists(mapping_idx_path):
    with open(mapping_idx_path, "r") as f:
      for line in filter(lambda x: x, map(lambda x: x.strip(), f)):
        k, v = line.split()
        mapping_dict[v] = k
  return mapping_dict


def read_chain_idx(data_uri, chain_dict=None):
  if chain_dict is None:
    chain_dict = {}
  chain_idx_path = os.path.join(data_uri.path, data_uri.chain_idx)
  if os.path.exists(chain_idx_path):
    with open(chain_idx_path, "r") as f:
      for line in filter(lambda x: x, map(lambda x: x.strip(), f)):
        k, *v = line.split()
        chain_dict[k] = v
  return chain_dict


def read_attrs_idx(data_uri, attr_dict=None):
  if attr_dict is None:
    attr_dict = {}
  attr_idx_path = os.path.join(data_uri.path, data_uri.attr_idx)
  if os.path.exists(attr_idx_path):
    with open(attr_idx_path, "r") as f:
      for line in filter(lambda x: x, map(lambda x: x.strip(), f)):
        k, *v = line.split()
        attr_dict[k] = " ".join(v)
  return attr_dict


def _read_fasta(data_dir, pid):
  with open(os.path.join(data_dir, "fasta", f"{pid}.fasta"), "r") as f:
    fasta_str = f.read()
  sequences, _ = parse_fasta(fasta_str)
  assert len(sequences) == 1
  return pid, sequences[0]


def read_fasta(data_uri, mapping_idx):
  fasta_list = {}
  with mp.Pool(processes=int(os.environ.get("NUM_PROCESSES", 16))) as p:
    f = functools.partial(_read_fasta, data_uri.path)
    for pid, sequence in p.imap(f,
                                set(mapping_idx.values()),
                                chunksize=int(os.environ.get("CHUNKSIZE",
                                                             256))):
      fasta_list[sequence] = pid
  return fasta_list


def align_peptide_main(args):  # pylint: disable=redefined-outer-name
  assert args.db
  seq_db, desc_db = [], []
  for db in args.db:
    with open(db, "r") as f:
      fasta_str = f.read()
    seqs, descs = parse_fasta(fasta_str)
    seq_db += seqs
    desc_db += descs
  desc_db = [desc.split()[0] for desc in desc_db]
  assert len(seq_db) == len(desc_db)
  seq_len = defaultdict(list)
  for i, seq in enumerate(seq_db):
    seq_len[len(seq)].append(i)

  for fasta_file in args.files:
    if args.verbose:
      print(f"process {fasta_file} ...")
    with open(fasta_file, "r") as f:
      fasta_string = f.read()
    sequences, descriptions = parse_fasta(fasta_string)
    assert len(sequences) == 1
    assert len(sequences) == len(descriptions)

    pid, _ = os.path.splitext(os.path.basename(fasta_file))
    output_path = os.path.join(args.output, pid, "msas")
    os.makedirs(output_path, exist_ok=True)
    with open(os.path.join(output_path, "bfd_uniclust_hits.a3m"), "w") as f:
      f.write(f">{descriptions[0]}\n")
      f.write(f"{sequences[0]}\n")

      k = len(sequences[0])
      if k in seq_len:
        for i in seq_len[k]:
          f.write(f">{desc_db[i]}/1-{k}\n")
          f.write(f"{seq_db[i]}\n")


def align_peptide_add_argument(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument("files", type=str, nargs="*", help="list of fasta files")
  parser.add_argument("-o",
                      "--output",
                      type=str,
                      default=".",
                      help="output dir, default=\".\"")
  parser.add_argument("--db",
                      type=str,
                      default=None,
                      nargs="+",
                      help="peptide fasta db, default=None")
  return parser


def read_a3m(data_uri, mapping_idx, pdb_id):
  pdb_id = mapping_idx.get(pdb_id, pdb_id)
  a3m_path = os.path.join(data_uri.path, "a3m", pdb_id, "msas", f"{pdb_id}.a3m")
  with open(a3m_path, "r") as f:
    a3m_string = f.read()
  sequences, descriptions = parse_fasta(a3m_string)
  return sequences, descriptions


def align_a3m(a3m_data, mapping_dict, align_dict, **kwargs):
  a3m_seqs, a3m_descs = a3m_data
  for seq, desc in zip(a3m_seqs[1:], a3m_descs[1:]):
    pid, chain, domains = decompose_pid(desc.split()[0], return_domain=True)
    domains = list(seq_index_split(domains))

    pdb_id = f"{pid}_{chain}" if chain else pid
    pdb_id_list = set([pdb_id]) | set(mapping_dict.get(pdb_id, []))
    for pdb_id in pdb_id_list:
      pid, chain = decompose_pid(pdb_id)  # pylint: disable=unbalanced-tuple-unpacking
      align_dict[pid].append((chain, domains, seq, desc, kwargs))
  return align_dict


def align_complex(item, **kwargs):
  target_pid, target_chain_list = item

  def _seq_at_i(a3m_data, i):
    seqs, descs = a3m_data
    return seqs[i], descs[i]

  # retrieve aligned chains
  a3m_list, a3m_dict = [], defaultdict(list)
  for chain in target_chain_list:
    pid = f"{target_pid}_{chain}" if chain else target_pid
    a3m_data = read_a3m(kwargs["target_uri"], kwargs["target_mapping_idx"], pid)
    a3m_dict = align_a3m(a3m_data,
                         kwargs["db_mapping_dict"],
                         a3m_dict,
                         target_chain=chain)
    a3m_list.append(a3m_data)

  # db_mapping_idx = kwargs["db_mapping_idx"]
  db_chain_idx = kwargs["db_chain_idx"]

  def _is_aligned(pid, chain_list):
    return pid in db_chain_idx and len(db_chain_idx[pid]) == len(chain_list)

  # filter complex with all chains aligned
  a3m_dict = {k: v for k, v in a3m_dict.items() if _is_aligned(k, v)}

  # realign the complex: iterate each target chain
  new_a3m_list = []

  # add target
  target_seq, domains = "", []
  n = 1
  for i, chain in enumerate(target_chain_list):
    seq, _ = _seq_at_i(a3m_list[i], 0)
    target_seq += seq
    domains.append((n, n + len(seq) - 1))
    n += len(seq) + 100
  domains = seq_index_join(domains)

  target_desc = f">{target_pid} domains:{domains}"
  new_a3m_list.append(target_desc)
  new_a3m_list.append(target_seq)

  def _repl(m):
    return "*" * len(m.group(0))

  # hit chains
  db_attr_idx = kwargs["db_attr_idx"]
  for pid, chain_list in a3m_dict.items():
    hit_desc = f">{pid} chains:{','.join(c for c, *_ in chain_list)}"
    if pid in db_attr_idx:
      hit_desc = f"{hit_desc} {db_attr_idx[pid]}"
    new_a3m_list.append(hit_desc)

    new_hit_seq = ""
    for i, target_chain in enumerate(target_chain_list):
      seq, _ = _seq_at_i(a3m_list[i], 0)
      hit_seq_at_i = "*" * len(seq)
      for _, _, hit_seq, _, attrs in chain_list:
        if attrs["target_chain"] == target_chain:
          hit_seq_at_i = hit_seq
          hit_seq_at_i = re.sub("^[-]+", _repl, hit_seq_at_i)
          hit_seq_at_i = re.sub("[-]+$", _repl, hit_seq_at_i)
          break
      new_hit_seq += hit_seq_at_i
    new_a3m_list.append(new_hit_seq)

  return item, new_a3m_list


def align_complex_main(args):  # pylint: disable=redefined-outer-name
  db_mapping_idx, db_chain_idx, db_attr_idx = {}, {}, {}
  for db_uri in args.db_uri:
    db_uri = parse_db_uri(db_uri)
    db_mapping_idx = read_mapping_idx(db_uri, db_mapping_idx)
    db_chain_idx = read_chain_idx(db_uri, db_chain_idx)
    db_attr_idx = read_attrs_idx(db_uri, db_attr_idx)
  db_mapping_dict = defaultdict(list)
  for k, v in db_mapping_idx.items():
    db_mapping_dict[v].append(k)

  target_uri = parse_db_uri(args.target_uri)
  target_mapping_idx = read_mapping_idx(target_uri)
  target_chain_idx = read_chain_idx(target_uri)

  with mp.Pool(processes=args.processes) as p:
    f = functools.partial(align_complex,
                          target_uri=target_uri,
                          target_mapping_idx=target_mapping_idx,
                          db_chain_idx=db_chain_idx,
                          db_mapping_idx=db_mapping_idx,
                          db_attr_idx=db_attr_idx,
                          db_mapping_dict=db_mapping_dict)
    for (pid, chain_list), a3m_list in p.imap(f,
                                              target_chain_idx.items(),
                                              chunksize=args.chunksize):
      output_path = os.path.join(args.output, pid, "msas")
      os.makedirs(output_path, exist_ok=True)
      with open(os.path.join(output_path, f"{pid}.a3m"), "w") as f:
        f.write("\n".join(a3m_list))
      print(f"{pid}\t{len(a3m_list)/2}\t{chain_list}")


def align_complex_add_argument(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument("files", type=str, nargs="*", help="list of fasta files")
  parser.add_argument("-o",
                      "--output",
                      type=str,
                      default=".",
                      help="output dir, default=\".\"")
  parser.add_argument("--db_uri",
                      type=str,
                      default=None,
                      nargs="+",
                      help="db uri.")
  parser.add_argument("--target_uri",
                      type=str,
                      default=None,
                      help="target uri.")
  parser.add_argument("--processes",
                      type=int,
                      default=None,
                      help="num of processes.")
  parser.add_argument("--chunksize", type=int, default=1, help="chunksize.")
  return parser


def csv_to_fasta_main(args):  # pylint: disable=redefined-outer-name
  os.makedirs(args.output, exist_ok=True)

  print(f"load {args.target_uri} ...")
  target_uri = parse_db_uri(args.target_uri)
  mapping_idx = read_mapping_idx(target_uri)
  fasta_idx = read_fasta(target_uri, mapping_idx)

  def cell_check(c):
    return c != "" and c.find("nan") == -1

  def cell_write(c, pid):
    if c not in fasta_idx:
      fasta_path = os.path.join(args.output, "fasta")
      os.makedirs(fasta_path, exist_ok=True)
      with open(os.path.join(fasta_path, f"{pid}.fasta"), "w") as f:
        f.write(f">{pid}\n")
        f.write(c)
      mapping_idx[pid] = pid
      fasta_idx[c] = pid
    else:
      mapping_idx[pid] = fasta_idx[c]

  print(f"process {args.csv_file} ...")
  with open(args.csv_file, "r") as f:
    reader = csv.DictReader(f)

    for i, row in enumerate(reader, start=args.start_idx):
      for key, chain in (("Antigen", "P"), ("MHC_str", "M"), ("a_seq", "A"),
                         ("b_seq", "B")):
        if key in row:
          if cell_check(row[key]):
            cell_write(row[key], f"{args.pid_prefix}{i}_{chain}")

  print("write mapping.idx ...")
  with open(os.path.join(args.output, "mapping.idx"), "w") as f:
    for v, k in mapping_idx.items():
      f.write(f"{k}\t{v}\n")


def csv_to_fasta_add_argument(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument("-o",
                      "--output",
                      type=str,
                      default=".",
                      help="output dir.")
  parser.add_argument("--target_uri", type=str, default=".", help="target dir.")
  parser.add_argument("--start_idx",
                      type=int,
                      default=0,
                      help="start index for each protein.")
  parser.add_argument("--pid_prefix",
                      type=str,
                      default="tcr_pmhc_",
                      help="pid prefix.")
  parser.add_argument("csv_file", type=str, default=None, help="csv file")
  return parser


def create_negative_main(args):  # pylint: disable=redefined-outer-name
  random.seed()

  mhc_seq_dict = {}

  peptides = set()
  tcr_mhc_rows = defaultdict(set)
  print(f"process {args.csv_file} ...")
  with open(args.csv_file, "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
      tcr_mhc_rows[(row["a_seq"], row["b_seq"], row["MHC_str"])].add(
          row["Antigen"])
      peptides.add(row["Antigen"])

  print("write {args.output} ...")
  with open(args.output, "w") as f:
    writer = csv.DictWriter(f, ("Antigen", "a_seq", "b_seq", "MHC_str"))
    writer.writeheader()
    for i, (row, antigens) in enumerate(tcr_mhc_rows.items()):
      a_seq, b_seq, mhc_str = row
      if args.verbose:
        print(f"{i}\t{len(antigens)}/{len(peptides)}")
      negatives = list(peptides - antigens)
      if negatives:
        random.shuffle(negatives)
        for antigen in negatives[:max(1, int(len(antigens) * args.amplify))]:
          writer.writerow({
              "Antigen": antigen,
              "a_seq": a_seq,
              "b_seq": b_seq,
              "MHC_str": mhc_str
          })


def create_negative_add_argument(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument("-o",
                      "--output",
                      type=str,
                      default=".",
                      help="output dir.")
  parser.add_argument("-n",
                      "--amplify",
                      type=float,
                      default=1.0,
                      help="amplify.")
  parser.add_argument("csv_file", type=str, default=None, help="csv file")
  return parser


def mhc_preprocess_main(args):  # pylint: disable=redefined-outer-name
  mhc_seq_dict = {}

  print(f"load {args.mhc_seq_file} ...")
  with open(args.mhc_seq_file, "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
      name = re.split(r"[*:]", row["name"])

      for idx in (3, 2):
        k = "".join(name[:idx])
        if k not in mhc_seq_dict:
          mhc_seq_dict[k] = row["sqe"]

  mhc_rows = []
  print(f"process {args.csv_file} ...")
  with open(args.csv_file, "r") as f:
    reader = csv.DictReader(f)

    for row in reader:
      allele = row["Allele"]
      if allele in mhc_seq_dict:
        mhc_rows.append({
            "Antigen": row["Peptide"],
            "MHC_str": mhc_seq_dict[allele]
        })
      elif args.verbose:
        print(f"{allele} not found")

  print("write {args.output} ...")
  with open(args.output, "w") as f:
    writer = csv.DictWriter(f, ("Antigen", "a_seq", "b_seq", "MHC_str"))
    writer.writeheader()
    for row in mhc_rows:
      writer.writerow(row)


def mhc_preprocess_add_argument(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument("-o",
                      "--output",
                      type=str,
                      default=".",
                      help="output dir.")
  parser.add_argument("--mhc_seq_file",
                      type=str,
                      default=None,
                      help="mhc sequence file.")
  parser.add_argument("csv_file", type=str, default=None, help="csv file")
  return parser


def split_data_main(args):  # pylint: disable=redefined-outer-name
  random.seed()

  def mapping_filter(k):
    _, chain = decompose_pid(k)  # pylint: disable=unbalanced-tuple-unpacking
    return chain == args.cluster_chain

  print(f"load {args.target_uri} ...")
  target_uri = parse_db_uri(args.target_uri)
  mapping_idx = read_mapping_idx(target_uri)
  mapping_idx = {k: v for k, v in mapping_idx.items() if mapping_filter(k)}
  fasta_idx = read_fasta(target_uri, mapping_idx)
  fasta_idx = {v: k for k, v in fasta_idx.items()}  # pid -> fasta
  chain_idx = read_chain_idx(target_uri)

  cluster_idx = {}
  with open(args.cluster_csv_file, "r") as f:
    for row in csv.DictReader(f):
      cluster_idx[row["Antigen"]] = row["Cluster"]

  # split by cluster id
  cluster_list = list(set(cluster_idx.values()))
  random.shuffle(cluster_list)
  print(f"Total clusters: {len(cluster_list)}")
  test_cluster_list = set(
      cluster_list[:int(len(cluster_list) * args.test_ratio)])
  print(f"Clusters for test: {len(test_cluster_list)}")
  with open(os.path.join(target_uri.path, "cluster.idx"), "w") as f:
    for seq, c in cluster_idx.items():
      if c in test_cluster_list:
        f.write(f"{seq},{c},test\n")
      else:
        f.write(f"{seq},{c},train\n")

  for pid, chain_list in chain_idx.items():
    label = "train"
    if args.cluster_chain in chain_list:
      k = f"{pid}_{args.cluster_chain}"
      k = mapping_idx.get(k, k)
      c = fasta_idx[k]
      if c in cluster_idx and cluster_idx[c] in test_cluster_list:
        label = f"test\t{cluster_idx[c]}"
    print(f"split_data\t{pid} {' '.join(chain_list)}\t{label}")


def split_data_add_argument(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument("--target_uri", type=str, default=".", help="target uri.")
  parser.add_argument("--test_ratio",
                      type=float,
                      default=0.2,
                      help="test set ratio.")
  parser.add_argument("--cluster_chain",
                      type=str,
                      default="P",
                      help="cluster by this chain only.")
  parser.add_argument("cluster_csv_file",
                      type=str,
                      default=None,
                      help="cluster csv file")
  return parser


if __name__ == "__main__":
  import argparse

  commands = {
      "align_peptide": (align_peptide_main, align_peptide_add_argument),
      "align_complex": (align_complex_main, align_complex_add_argument),
      "csv_to_fasta": (csv_to_fasta_main, csv_to_fasta_add_argument),
      "create_negatives": (create_negative_main, create_negative_add_argument),
      "mhc_preprocess": (mhc_preprocess_main, mhc_preprocess_add_argument),
      "split_data": (split_data_main, split_data_add_argument),
  }

  formatter_class = argparse.ArgumentDefaultsHelpFormatter
  parser = argparse.ArgumentParser(formatter_class=formatter_class)

  sub_parsers = parser.add_subparsers(dest="command", required=True)
  for cmd, (_, add_argument) in commands.items():
    cmd_parser = sub_parsers.add_parser(cmd, formatter_class=formatter_class)
    cmd_parser = add_argument(cmd_parser)
    cmd_parser.add_argument("-v",
                            "--verbose",
                            action="store_true",
                            help="verbose")

  args = parser.parse_args()
  main, _ = commands[args.command]
  main(args)
