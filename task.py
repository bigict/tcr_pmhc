"""task def"""
from collections import defaultdict

task_mapping = {
    0: [(all, ["P", "M"])],                     # pMHC
    1: [(all, ["P"]), (any, ["A", "B"])],       # pTCR
    2: [(all, ["P", "M"]), (any, ["A", "B"])]   # TCR_pMHC
}

task_num = len(task_mapping.keys())

task_name_list = ["pMHC", "pTCR", "TCR_pMHC"]

def make_def():
  task_def = defaultdict(list)
  for task_idx, op_list in task_mapping.items():
    for _, chains in op_list:
      for chain in chains:
        task_def[chain].append(task_idx)
  return task_def


def _task_label_mask(task_idx, chain_list):
  for op, chains in task_mapping[task_idx]:
    if not op(chain in chain_list for chain in chains):
      return False
  return True


def make_label(label, chain_list):
  label_mask = [_task_label_mask(task_idx, chain_list) for task_idx in range(task_num)]

  # FIX: pMHC=1 for all TCR_pMHCs, task[0] is pMHC
  assert label is not None
  if label > 0:
    label = list(map(lambda x: float(x) * label, label_mask))
  elif any(label_mask[1:]):
    label = [float(label_mask[0]) * 1.0] + list(
        map(lambda x: float(x) * label, label_mask[1:])
    )
  else:
    assert not any(label_mask[1:])
    label = [label] + [0.0] * (task_num - 1)

  return label, label_mask
