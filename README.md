# optimalCytosolicDensity_wholeCellModel

## A Model to Investigate the Optimal Cytosolic Density of a Cell

This repository stores the python source code of the 'whole cell model' that investigates how a bacterial cell optimizes its cytosolic density.

## Full Documentation of the Model

Pang, T. Y., & Lercher, M. J. (2020). Optimal density of bacterial cells. Cold Spring Harbor Laboratory. https://doi.org/10.1101/2020.11.18.388744

See Fig. 2 of the manuscript for an illustration of the model.

## Description of the Files

The python3 script encodes the model; it requires the numpy and scipy python package. The file 'save_command_lines.txt' contains the shell commands to run the script; it contains the random initial conditions fed to the python script.

### Inputs of python script

The script takes 9 command line arguments as input. For example:
```shell 
python3 script_wholeCellModel.py rho_ratio s_ext Nrepeat Nm input_s input_p input_T input_M input_R
```
1. rho_ratio: occupancy of cell, i.e. ratio of cytosolic volume occupied by dry mass
2. s_ext: concentration substrate s in the environment; unit: µM
3. Nrepeat: number of attempts to solve the problem
4. Nm: number of steps in the metabolic pathway
5. input_s: initial substrate concentration; unit: log(copy-number-per-cubic-micron)
6. input_p: initial precursor concentration; unit: log(copy-number-per-cubic-micron)
7. input_T: initial transporter T per cell; unit: log(copy-number-per-cell)
8. input_M: initial enzyme M concentration; unit: log(copy-number-per-cubic-micron)
9. input_R: initial ribosome R concentration; unit: log(copy-number-per-cubic-micron)

### Outputs of python script

The script prints out a table with 17 columns; the number of rows is Nrepeat, i.e. the number of attempts to solve for the optimal growth rate given the circumstances. The different columns correspond to: 
1. row index
2. substrate concentration in the environment
3. occupancy rho
4. optimal growth rate found in the current attempt
5. K<sub>M</sub><sup>*</sup> of metabolic reaction
6. K<sub>M</sub><sup>*</sup> of ribosomal reaction
7. dummy output
8. log(concentraion<sup>‡</sup>) of substrate s
9. log(concentraion<sup>‡</sup>) of precursor p
10. log(copy-number) of transporter T per cell
11. log(concentraion<sup>‡</sup>) of metabolic enzyme M
12. log(concentraion<sup>‡</sup>) of ribosome R
13. volume fracton of substrate s
14. volume fracton of prevursor p
15. dummay output
16. volume fracton of metabolic enzyme M
17. volume fracton of ribosome R

<sup>‡</sup>unit of concentraion: copy number per cubic micron

The code was tested on a machine with Debian GNU/Linux 10 (buster), installed with python 3.7.3, numpy 1.21.5, and scipy 1.7.3.
