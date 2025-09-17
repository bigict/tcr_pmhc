![header](assets/fig1.png)

# ProTCR

![GitHub License](https://img.shields.io/github/license/bigict/ProTCR)
[![Conventional Commits](https://img.shields.io/badge/Conventional%20Commits-1.0.0-yellow.svg)](https://conventionalcommits.org)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/bigict/ProTCR/pylint.yml)
![GitHub Release](https://img.shields.io/github/v/release/bigict/ProTCR)


This package provides an implementation of the inference pipeline of [ProTCR](https://github.com/bigict/ProTCR). 

## Requirements

* [Python3.9+](https://www.python.org)
* [ProFOLD2](https://github.com/bigict/ProFOLD2)
* [Dependencies](https://github.com/bigict/tcr_pmhc/network/dependencies)

## Installation

1.  Clone this repository and `cd` into it.

  ```bash
  git clone https://github.com/bigict/ProTCR.git
  cd ProTCR
  git submodule update --init --remote
  ```

2. Create conda an environment and activate it.

  ```bash
  conda create -n tcr python=3.11
  conda activate tcr
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

1. Model Parameters

   The [ProTCR](https://github.com/bigict/ProTCR) model parameters are hosted on [Hugging Face](https://huggingface.co/bigict/ProTCR).

   ```bash
   huggingface-cli download --local-dir params bigict/ProTCR
   
   ```

2. Inference

  Once you have installed ProTCR, you can test your setup using e.g. the following input `csv` file named `tcr_pmhc_input.csv`

  ```csv
  Antigen,TCRB,TCRA,MHC_str
  ALSKGVHFV,VSQHPSWVICKSGTSVKIECRSLDFQATTMFWYRQFPKQSLMLMATSNEGSKATYEQGVEKDKFLINHASLTLSTLTVTSAHPEDSSFYICSARGSSGRAEYTQYFGPGTRLTVLE,,SHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTD
  ```

  ```bash
  ./predict.sh --output_dir [OUTPUT_DIR] tcr_pmhc_input.csv
  ```
  
  You can run
  ```bash
  ./predict.sh --help
  ```
  for further help.

## Known Issues

  Please [create an issue](https://github.com/bigict/ProTCR/issues/new) if it is not already listed in the [issues tracker](https://github.com/bigict/ProTCR/issues)

