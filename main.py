"""Tools for preprocess TCR_pMHC data"""
import os
import re
from collections import defaultdict
import csv
from dataclasses import dataclass
import functools
import json
import multiprocessing as mp
from urllib.parse import urlparse, parse_qs

import click
import sqlitedict

from profold2.data.parsers import parse_fasta, parse_a3m
from profold2.data.utils import compose_pid, decompose_pid, seq_index_join, seq_index_split
from profold2.utils import timing

import task


class DictObject(object):
  def __init__(self, **args):
    self.__dict__.update(args)


@dataclass
class DBUri:
  path: str
  chain_idx: str = "chain.idx"
  mapping_idx: str = "mapping.idx"
  attr_idx: str = "attr.idx"
  a3m_dir: str = "a3m"
  append: bool = True


def parse_db_uri(db_uri):
  o = urlparse(db_uri)
  chain_idx, mapping_idx, attr_idx = "chain.idx", "mapping.idx", "attr.idx"
  a3m_dir = "a3m"
  append = True
  if o.query:
    attrs = parse_qs(o.query)
    if "chain_idx" in attrs:
      chain_idx = attrs["chain_idx"][-1]
    if "mapping_idx" in attrs:
      mapping_idx = attrs["mapping_idx"][-1]
    if "attr_idx" in attrs:
      attr_idx = attrs["attr_idx"][-1]
    if "a3m_dir" in attrs:
      a3m_dir = attrs["a3m_dir"][-1]
    if "append" in attrs:
      append = bool(int(attrs["append"][-1]))
  return DBUri(
      o.path, chain_idx, mapping_idx, attr_idx=attr_idx, a3m_dir=a3m_dir, append=append
  )


def _db_uri_abs_path(data_uri, data_idx):
  if os.path.isabs(data_idx):
    return data_idx
  return os.path.join(data_uri.path, data_idx)


def read_mapping_idx(data_uri, mapping_dict=None):
  if mapping_dict is None:
    mapping_dict = {}
  # mapping_idx_path = os.path.join(data_uri.path, data_uri.mapping_idx)
  mapping_idx_path = _db_uri_abs_path(data_uri, data_uri.mapping_idx)
  if data_uri.append and os.path.exists(mapping_idx_path):
    with open(mapping_idx_path, "r") as f:
      for line in filter(lambda x: x, map(lambda x: x.strip(), f)):
        k, v = line.split()
        mapping_dict[v] = k
  return mapping_dict


def read_chain_idx(data_uri, chain_dict=None):
  if chain_dict is None:
    chain_dict = {}
  # chain_idx_path = os.path.join(data_uri.path, data_uri.chain_idx)
  chain_idx_path = _db_uri_abs_path(data_uri, data_uri.chain_idx)
  if data_uri.append and os.path.exists(chain_idx_path):
    with open(chain_idx_path, "r") as f:
      for line in filter(lambda x: x, map(lambda x: x.strip(), f)):
        k, *v = line.split()
        chain_dict[k] = v
  return chain_dict


def read_attrs_idx(data_uri, attr_dict=None):
  if attr_dict is None:
    attr_dict = {}
  # attr_idx_path = os.path.join(data_uri.path, data_uri.attr_idx)
  attr_idx_path = _db_uri_abs_path(data_uri, data_uri.attr_idx)
  if data_uri.append and os.path.exists(attr_idx_path):
    with open(attr_idx_path, "r") as f:
      for line in filter(lambda x: x, map(lambda x: x.strip(), f)):
        k, *v = line.split()
        try:
          attr_dict[k] = json.loads(" ".join(v))
        except:
          pass
  return attr_dict


def read_fasta(data_uri, **kwargs):
  return sqlitedict.open(os.path.join(data_uri.path, "fasta.db"), **kwargs)


def create_shared_obj(**kwargs):
  context = mp.get_context()
  manager = context.Manager()
  shared_obj = manager.Namespace()
  for k, v in kwargs.items():
    setattr(shared_obj, k, v)
  return shared_obj


@click.group()
def main():
  pass


@main.command("peptide_align")
@click.argument("fasta_file", type=click.Path(), nargs=-1)
@click.option("-o", "--output_dir", type=str, default=".", help="output dir.")
@click.option("--target_db", type=click.Path(), multiple=True, help="db file.")
@click.option("-v", "--verbose", is_flag=True, help="verbose output.")
def peptide_align(**args):  # pylint: disable=redefined-outer-name
  args = DictObject(**args)

  assert args.target_db
  seq_db, desc_db = [], []
  for db in args.target_db:
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

  for fasta_file in args.fasta_file:
    if args.verbose:
      print(f"process {fasta_file} ...")
    with open(fasta_file, "r") as f:
      fasta_string = f.read()
    sequences, descriptions = parse_fasta(fasta_string)
    assert len(sequences) == 1
    assert len(sequences) == len(descriptions)

    pid, _ = os.path.splitext(os.path.basename(fasta_file))
    output_path = os.path.join(args.output_dir, pid, "msas")
    os.makedirs(output_path, exist_ok=True)
    with open(os.path.join(output_path, "bfd_uniclust_hits.a3m"), "w") as f:
      f.write(f">{descriptions[0]}\n")
      f.write(f"{sequences[0]}\n")

      k = len(sequences[0])
      if k in seq_len:
        for i in seq_len[k]:
          f.write(f">{desc_db[i]}/1-{k}\n")
          f.write(f"{seq_db[i]}\n")


def read_a3m(data_uri, mapping_idx, pdb_id):
  pdb_id = mapping_idx.get(pdb_id, pdb_id)
  a3m_path = os.path.join(
      _db_uri_abs_path(data_uri, data_uri.a3m_dir), pdb_id, "msas", f"{pdb_id}.a3m"
  )
  with open(a3m_path, "r") as f:
    a3m_string = f.read()
  sequences, descriptions = parse_fasta(a3m_string)
  return sequences, descriptions


def complex_align_read_a3m(a3m_data, mapping_dict, align_dict, **kwargs):
  a3m_seqs, a3m_descs = a3m_data

  align_data_dict, align_chain_dict = align_dict
  for seq, desc in zip(a3m_seqs[1:], a3m_descs[1:]):
    pid, chain, domains = decompose_pid(desc.split()[0], return_domain=True)
    if domains:
      domains = list(seq_index_split(domains))
    else:
      domains = []

    pdb_id = compose_pid(pid, chain)
    align_data_dict[pdb_id] = (domains, seq, desc, kwargs)

    pdb_id_list = set([pdb_id]) | set(mapping_dict.get(pdb_id, []))
    for pdb_id in pdb_id_list:
      pid, chain = decompose_pid(pdb_id)  # pylint: disable=unbalanced-tuple-unpacking
      # align_dict[pid].append((chain, domains, seq, desc, kwargs))
      align_chain_dict[pid].append(chain)
  # return align_dict
  return align_data_dict, align_chain_dict


_db_mapping_idx, _db_chain_idx, _db_attr_idx = {}, {}, {}
_db_mapping_dict = defaultdict(list)


def complex_align_init(*db_uri_list):
  for db_uri in db_uri_list:
    db_uri = parse_db_uri(db_uri)
    _db_mapping_idx = read_mapping_idx(db_uri, _db_mapping_idx)
    _db_chain_idx = read_chain_idx(db_uri, _db_chain_idx)
    _db_attr_idx = read_attrs_idx(db_uri, _db_attr_idx)
  for k, v in _db_mapping_idx.items():
    _db_mapping_dict[v].append(k)


def complex_align_func(output_dir, target_uri, target_mapping_idx, item):
  target_pid, target_chain_list = item

  def _seq_at_i(a3m_data, i):
    seqs, descs = a3m_data
    return seqs[i], descs[i]

  # retrieve aligned chains
  a3m_list, a3m_dict = [], ({}, defaultdict(list))
  with timing(f"read_a3m_dict ({target_pid})", print_fn=print):
    for chain in target_chain_list:
      pid = f"{target_pid}_{chain}" if chain else target_pid
      a3m_data = read_a3m(target_uri, target_mapping_idx, pid)
      with timing(f"align_a3m ({target_pid}_{chain})", print_fn=print):
        a3m_dict = complex_align_read_a3m(
            a3m_data, _db_mapping_dict, a3m_dict, target_chain=chain
        )
      a3m_list.append(a3m_data)

  def _is_aligned(pid, chain_list):
    if pid in _db_attr_idx:
      return pid in _db_chain_idx and len(_db_chain_idx[pid]) == len(chain_list)
    return False

  align_data_dict, a3m_dict = a3m_dict
  # filter complex with all chains aligned
  with timing(f"filter a3m_dict ({target_pid})", print_fn=print):
    new_a3m_dict = {k: v for k, v in a3m_dict.items() if _is_aligned(k, v)}

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
  for pid, chain_list in new_a3m_dict.items():
    hit_desc = f">{pid} chains:{','.join(c for c, *_ in chain_list)}"
    if pid in _db_attr_idx:
      hit_desc = f"{hit_desc} {_db_attr_idx[pid]}"
    new_a3m_list.append(hit_desc)

    new_hit_seq = ""
    for i, target_chain in enumerate(target_chain_list):
      seq, _ = _seq_at_i(a3m_list[i], 0)
      hit_seq_at_i = "*" * len(seq)
      # for _, _, hit_seq, _, attrs in chain_list:
      for chain in chain_list:
        pdb_id = f"{pid}_{chain}" if chain else pid
        assert pdb_id in _db_mapping_idx, (target_pid, pdb_id)
        if pdb_id in align_data_dict:
          _, hit_seq, _, attrs = align_data_dict[pdb_id]
        else:
          assert _db_mapping_idx[pdb_id] in align_data_dict
          _, hit_seq, _, attrs = align_data_dict[_db_mapping_idx[pdb_id]]
        if attrs["target_chain"] == target_chain:
          hit_seq_at_i = hit_seq
          hit_seq_at_i = re.sub("^[-]+", _repl, hit_seq_at_i)
          hit_seq_at_i = re.sub("[-]+$", _repl, hit_seq_at_i)
          break
      new_hit_seq += hit_seq_at_i
    new_a3m_list.append(new_hit_seq)

  output_path = os.path.join(output_dir, target_pid, "msas")
  os.makedirs(output_path, exist_ok=True)
  with open(os.path.join(output_path, f"{target_pid}.a3m"), "w") as f:
    f.write("\n".join(new_a3m_list))

  return item, len(new_a3m_list)


@main.command("complex_align")
@click.argument("fasta_file", type=click.Path(), nargs=-1)
@click.option("-o", "--output_dir", type=click.Path(), default=".", help="output dir.")
@click.option("--query_uri", type=str, default=".", help="query dir.")
@click.option("--target_uri", type=str, multiple=True, help="target dir.")
@click.option("--processes", type=int, default=None, help="num of processes.")
@click.option("--chunksize", type=int, default=2, help="chunksize.")
@click.option("-v", "--verbose", is_flag=True, help="verbose output.")
def complex_align(**args):  # pylint: disable=redefined-outer-name
  args = DictObject(**args)

  query_uri = parse_db_uri(args.query_uri)
  query_mapping_idx = read_mapping_idx(query_uri)
  query_chain_idx = read_chain_idx(query_uri)

  with mp.Pool(
      processes=args.processes, initializer=complex_align_init, initargs=args.target_uri
  ) as p:
    f = functools.partial(
        complex_align_func, args.output_dir, query_uri, query_mapping_idx
    )
    for (pid, chain_list), a3m_list in p.imap_unordered(
        f, query_chain_idx.items(), chunksize=args.chunksize
    ):
      print(f"{pid}\t{a3m_list/2}\t{chain_list}")


@main.command("csv_to_fasta")
@click.argument("csv_file", type=click.Path(), nargs=-1)
@click.option("--target_uri", type=str, default=".", help="target dir.")
@click.option("--start_idx", type=int, default=0, help="start index for each protein.")
@click.option("--pid_prefix", type=str, default="tcr_pmhc_", help="pid prefix.")
@click.option("--default_y", type=float, default=None, help="default label.")
@click.option("-v", "--verbose", is_flag=True, help="verbose output.")
def csv_to_fasta(**args):
  args = DictObject(**args)

  print(f"load {args.target_uri} ...")
  target_uri = parse_db_uri(args.target_uri)
  os.makedirs(target_uri.path, exist_ok=True)
  mapping_idx = read_mapping_idx(target_uri)
  attr_idx = read_attrs_idx(target_uri)

  with sqlitedict.open(
      os.path.join(target_uri.path, "fasta.db"), autocommit=False
  ) as fasta_db:
    fasta_idx = {seq: pid for pid, seq in fasta_db.items()}

    def cell_check(c):
      return c != "" and c.find("nan") == -1

    def cell_write(c, pid):
      if c not in fasta_idx:
        fasta_db[pid] = c
        mapping_idx[pid] = pid
        fasta_idx[c] = pid
      else:
        mapping_idx[pid] = fasta_idx[c]

    for csv_idx, csv_file in enumerate(args.csv_file):
      print(f"process [{csv_idx}] {csv_file} ...")
      with open(csv_file, "r") as f:
        reader = csv.DictReader(f)

        for i, row in enumerate(reader, start=args.start_idx):
          pdb_id = f"{args.pid_prefix}{csv_idx}_{i}"

          label = None
          if "y" in row:
            label = float(row["y"])
          elif "label" in row:
            label = float(row["label"])
          elif args.default_y is not None:
            label = args.default_y

          chain_list = set()
          for key, chain in (
              ("Antigen", "P"),
              ("Peptide", "P"),
              ("MHC_str", "M"),
              ("a_seq", "A"),
              ("b_seq", "B"),
              ("tcrb", "B"),
              ("TCRA", "A"),
              ("TCRB", "B")
          ):
            if key in row:
              if cell_check(row[key]):
                cell_write(row[key], f"{pdb_id}_{chain}")
                chain_list.add(chain)

          label, label_mask = task.make_label(label, chain_list)

          attr_idx[pdb_id] = {"label": label, "label_mask": label_mask}
          for key in ("HLA", "Allele"):
            if key in row:
              if pdb_id in attr_idx:
                attr_idx[pdb_id]["MHC"] = row[key]
              else:
                attr_idx[pdb_id] = {"MHC": row[key]}
              break

    fasta_db.commit()

  if args.verbose:
    print(f"write {target_uri.mapping_idx} ...")
  with open(os.path.join(target_uri.path, target_uri.mapping_idx), "w") as f:
    for v, k in mapping_idx.items():
      f.write(f"{k}\t{v}\n")
  if args.verbose:
    print(f"write {target_uri.attr_idx} ...")
  with open(os.path.join(target_uri.path, target_uri.attr_idx), "w") as f:
    for k, v in attr_idx.items():
      v = json.dumps(v)
      f.write(f"{k}\t{v}\n")


@main.command("sampling_weight")
@click.option("--target_uri", type=str, default=".", help="target uri.")
@click.option("--pid_topk", type=int, default=1, help="pid topk.")
@click.option("-v", "--verbose", is_flag=True, help="verbose output.")
def sampling_weight(**args):  # pylint: disable=redefined-outer-name
  args = DictObject(**args)

  print(f"load {args.target_uri} ...")
  target_uri = parse_db_uri(args.target_uri)
  chain_idx = read_chain_idx(target_uri)
  mapping_idx = read_mapping_idx(target_uri)
  attr_idx = read_attrs_idx(target_uri)

  new_mapping_dict = defaultdict(list)
  for pid, chain_list in chain_idx.items():
    if pid in attr_idx:
      if "P" in chain_list and "M" in chain_list:
        new_mapping_dict[(mapping_idx[f"{pid}_P"], mapping_idx[f"{pid}_M"])].append(pid)

  for (pid_p, pid_m), new_pid_list in new_mapping_dict.items():
    print(f"new_pid_list_count: {pid_p} {pid_m} {len(new_pid_list)}")

    m = len(new_pid_list)
    if 0 < args.pid_topk < m:
      m = args.pid_topk

    for i, new_pid in enumerate(new_pid_list):
      weight = m / len(new_pid_list)
      attr_idx[new_pid].update(weight=weight)

  print(f"write {target_uri.attr_idx} ...")
  with open(os.path.join(target_uri.path, target_uri.attr_idx), "w") as f:
    for k, v in attr_idx.items():
      v = json.dumps(v)
      f.write(f"{k}\t{v}\n")


@main.command("a3m_filter")
@click.argument("fasta_file", type=click.Path(), nargs=-1)
@click.option("-o", "--output_dir", type=click.Path(), default=".", help="output dir.")
@click.option(
    "-t",
    "--aligned_ratio_threshold",
    type=float,
    default=0.5,
    help="aligned ratio threshold."
)
@click.option("--trim_gap", is_flag=True, help="trim gap.")
@click.option("-v", "--verbose", is_flag=True, help="verbose output.")
def a3m_filter(**args):
  args = DictObject(**args)

  def _aligned_ratio(msa, query=None, n=None, trim_gap=True):
    i, j = 0, len(msa)
    if trim_gap:
      while i < len(msa) and msa[i] == "-":
        i += 1
      while j > 0 and msa[j - 1] == "-":
        j -= 1
    if not n:
      n = len(msa[i:j])
    if query:
      c = 0
      while i < j:
        if query[i] != msa[i]:
          c += 1
        i += 1
    else:
      c = msa[i:j].count("-")
    r = 1. - c / n
    return r

  for fasta_file in args.fasta_file:
    print(f"processing {fasta_file} ...")
    pid, _ = os.path.splitext(os.path.basename(fasta_file))
    a3m_file = os.path.join(args.output_dir, pid, "msas", f"{pid}.a3m")
    with open(a3m_file, "r") as f:
      a3m_string = f.read()

    sequences, descriptions = parse_fasta(a3m_string)
    aligned_seqs, _ = parse_a3m(a3m_string)

    with open(a3m_file, "w") as f:
      n = 0
      for i, seq in enumerate(aligned_seqs):
        msa_aligned_ratio = _aligned_ratio(
            seq, query=sequences[0], trim_gap=args.trim_gap
        )
        if msa_aligned_ratio >= args.aligned_ratio_threshold:
          f.write(f">{descriptions[i]} aligned_ratio:{msa_aligned_ratio}\n")
          f.write(f"{sequences[i]}\n")
          n += 1
      print(f"filtering {fasta_file} {n}/{len(sequences)}")


@main.command("fasta_extract")
@click.option("--target_uri", type=str, default=".", help="target dir.")
@click.option("--chain", type=str, multiple=True, help="chain.")
@click.option("-v", "--verbose", is_flag=True, help="verbose output.")
def fasta_extract(**args):
  args = DictObject(**args)

  target_uri = parse_db_uri(args.target_uri)
  with sqlitedict.open(
      os.path.join(target_uri.path, "fasta.db"), autocommit=False
  ) as fasta_db:
    for pid, seq in fasta_db.items():
      k = pid.rfind("_")
      assert k != -1
      chain = pid[k + 1:]
      if chain in args.chain:
        print(f">{pid}\n{seq}")


@main.command("attr_update_weight_and_task")
@click.argument("attr_idx", type=click.Path())
@click.option("--weight", type=float, default=1.0, help="weight.")
@click.option("-v", "--verbose", is_flag=True, help="verbose output.")
def attr_update_weight_and_task(**args):
  args = DictObject(**args)

  task_def = task.make_def()
  with open(args.attr_idx, "r") as f:
    for line in filter(lambda x: x, map(lambda x: x.strip(), f)):
      k, *v = line.split()
      v = json.loads(" ".join(v))

      v["weight"] = args.weight
      v["task_def"] = task_def

      v = json.dumps(v)
      print(f"{k}\t{v}")


if __name__ == "__main__":
  main()
