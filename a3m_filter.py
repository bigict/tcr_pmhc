import os

from profold2.data.parsers import parse_fasta, parse_a3m


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

def main(args):
  for fasta_file in args.files:
    print(f"processing {fasta_file} ...")
    pid, _ = os.path.splitext(os.path.basename(fasta_file))
    a3m_file = os.path.join(args.output, pid, 'msas', f'{pid}.a3m')
    with open(a3m_file, "r") as f:
      a3m_string = f.read()

    sequences, descriptions = parse_fasta(a3m_string)
    aligned_seqs, _ = parse_a3m(a3m_string)

    with open(a3m_file, "w") as f:
      n = 0
      for i, seq in enumerate(aligned_seqs):
        msa_aligned_ratio = _aligned_ratio(seq, query=sequences[0], trim_gap=args.trim_gap)
        if msa_aligned_ratio >= args.aligned_ratio_threshold:
          f.write(f">{descriptions[i]} aligned_ratio:{msa_aligned_ratio}\n")
          f.write(f"{sequences[i]}\n")
          n += 1
      print(f"filtering {fasta_file} {n}/{len(sequences)}")

if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument('files', type=str, nargs='*',
                      help='list of fasta files')
  parser.add_argument('-o', '--output', type=str, default='.',
                      help='output dir, default=\'.\'')
  parser.add_argument('-t', '--aligned_ratio_threshold', type=float, default=0.5,
                      help='aligned ratio threshold.')
  parser.add_argument('--trim_gap', action='store_true',
                      help='trim gap.')
  parser.add_argument('-v', '--verbose', action='store_true',
                      help='verbose')
  args = parser.parse_args()

  main(args)
