# optimalCytosolicDensity_wholeCellModel

## A Model to Investigate the Optimal Cytosolic Density of a Cell

This repository stores the python source code of the 'whole cell model' that investigates how a bacterial cell optimizes its cytosolic density.

## Full Documentation of the Model

Pang, T. Y., & Lercher, M. J. (2020). Optimal density of bacterial cells. Cold Spring Harbor Laboratory. https://doi.org/10.1101/2020.11.18.388744

See Fig. 2 of the manuscript for an illustration of the model.

## Description of the Files

The python3 script encodes the model; it requires the numpy and scipy python package. The file 'save_command_lines.txt' contains the shell commands to run the script; it contains the random initial conditions fed to the python script.

### Inputs of python script

For example:
```shell 
python3 script_wholeCellModel.py rho_ratio s_ext Nrepeat Nm input_s input_p input_T input_M input_R
```
* rho_ratio: occupancy of cell, i.e. ratio of cytosolic volume occupied by dry mass
* s_ext: concentration substrate s in the environment; unit: µM
* input_s: initial substrate concentration; unit: log(copy-number-per-cubic-micron)
* input_p: initial precursor concentration; unit: log(copy-number-per-cubic-micron)
* input_T: initial transporter T per cell; unit: log(copy-number-per-cell)
* input_M: initial enzyme M concentration; unit: log(copy-number-per-cubic-micron)
* input_R: initial ribosome R concentration; unit: log(copy-number-per-cubic-micron)
* Nrepeat: number of attempts to solve the problem
* Nm: number of steps in the metabolic pathway

### Outputs of python script

The script outputs a table 16 columns; the number of rows is Nrepeat, i.e. the number of attempts to solve for the optimal growth rate given the circumstances.
The different columns correspond to: index, substrate concentration in the environment, occupancy rho, optimal growth rate found in the current attempt, K_M^* of metabolic reaction, K_M^* of ribosomal reaction, dummy output, log-concentraion** of substrate s, log-concentration** of precursor p, log-copy-number of transporter T per cell, log-concentration** of metabolic enzyme M, log-concentration** of ribosome R, volume fracton of substrate s, volume fracton of prevursor p, dummay output, volume fracton of metabolic enzyme M, volume fracton of ribosome R.

**unit of concentraion: copy number per cubic micron

The code was tested on a machine with Debian GNU/Linux 10 (buster), installed with python 3.7.3, numpy 1.21.5, and scipy 1.7.3.
