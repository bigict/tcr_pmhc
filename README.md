![header](assets/fig1.png)

# ProTCR

[![Conventional Commits](https://img.shields.io/badge/Conventional%20Commits-1.0.0-yellow.svg)](https://conventionalcommits.org)

This package provides an implementation of the inference pipeline of [ProTCR](https://github.com/bigict/tcr_pmhc). 

## Installation

1.  Clone this repository and `cd` into it.
  ```bash
  $git clone https://github.com/bigict/tcr_pmhc.git
  $cd tcr_pmhc
  $git submodule update --init --remote
  ```
2. Create conda an environment and activate it.
  ```bash
  $conda create -n tcr python=3.11
  $conda acrivate tcr
  ```

  You need to install the [ProFOLD2](https://github.com/bigict/ProFOLD2) dependencies
  ```bash
  bash profold2/install_env.sh
  ```

  And then
  ```
  pip install -r requirements.txt 
  ```

## Running ProTCR

1.  Inference
  ```bash
  $python main.py predict --models [MODEL_NAME1:]MODEL_FILE1 [MODEL_NAME2:]MODEL_FILE2
  ```
  
  Just like `train`, you can run
  ```bash
  $python main.py predict --help
  ```
  
2.  Train a model
  Create the 5-fold dataset first
  ```bash
  bash data/train_data_fold_5.sh
  ```

  ```bash
  $./env/bin/python main.py train --prefix=OUTPUT_DIR
  ```
  
  There are a lot of parameters, you can run
    
  ```bash
  $./env/bin/python main.py train -h
  ```
  
  for further help.
  
  `ProFOLD2` logs it's metrics to [TensorBoard](https://www.tensorflow.org/tensorboard). You can run
  
  ```bash
  $tensorboard --logdir=OUTPUT_DIR
  ```
  
  Then open http://localhost:6006 in you browser.
