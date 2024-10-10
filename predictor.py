import os
import logging

import torch
from torch.utils.data.distributed import DistributedSampler

from profold2.command.worker import main
from profold2.data.dataset import _make_var_features
from profold2.model.head import HeaderBuilder


def _load_models(rank, args):
  def _location_split(model_location):
    k = model_location.find("=")
    if k != -1:
      return model_location.split("=", 1)
    model_name = os.path.basename(model_location)
    model_name, _ = os.path.splitext(model_name)
    return model_name, model_location

  for i, model_location in enumerate(args.refs):
    model_name, model_location = _location_split(model_location)
    logging.info("Load model [%d/%d] %s from %s",
        i, len(refs), model_name, model_location)
    pth = torch.load(model_location)
    yield model_name, pth

def _create_dataloader(rank, args):
  pass

def predict(rank, args):
  dim_single, _dim_pairwise = 384, 128
  headers = [("fitness", {"mask": "-", "pooling": "mean"}, {})]

  headers = HeaderBuilder.build((_dim_single, _dim_pairwise), headers)

  model_runners = dict(_load_models(rank, args))
  # constract r.headers from refs
  print(headers)


def add_arguments(parser):  # pylint: disable=redefined-outer-name
  parser.add_argument("--refs", type=str, nargs="+",
                      help="list of refs files")
  parser.add_argument("--vars", type=str, default=None,
                      help="aligned ratio threshold.")


if __name__ == "__main__":
  import argparse

  formatter_class = argparse.ArgumentDefaultsHelpFormatter
  parser = argparse.ArgumentParser(formatter_class=formatter_class)

  # init distributed env
  parser.add_argument("--nnodes", type=int, default=None,
      help="number of nodes.")
  parser.add_argument("--node_rank", type=int, default=0,
      help="rank of the node.")
  parser.add_argument("--local_rank", type=int, default=None,
      help="local rank of xpu.")
  parser.add_argument("--init_method", type=str,
      default="file:///tmp/profold2-tcr_pmhc.dist",
      help="method to initialize the process group")

  # output dir
  parser.add_argument("-o", "--prefix", type=str, default=".",
      help="prefix of out directory.")
  add_arguments(parser)
  # verbose
  parser.add_argument("-v", "--verbose", action="store_true", help="verbose")

  args = parser.parse_args()
  main(args, predict)
