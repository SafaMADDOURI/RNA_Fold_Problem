
# Creation of an objective function for the RNA folding problem

For a given ribonucleotide chain, the RNA folding problem consists in finding the native fold among the astronomically large number of possible conformations. The native fold is the one with the lowest Gibbs free energy, the objective function should be an estimator of this energy.

For this purpose, three puthon scripts were implemented.
## Script1: Training script
The objective of the training script is to:
-Compute the interatomic distances from a given dataset of PDB files.Only C3’ atoms, “intrachain” distances and residues separated by at least 3 positions on the sequence are considered.
-Compute the observed frequencies: 10 x 21 distance intervals (0 to 21 A)
-Compute the reference frequency (=the "XX" pair)
-Compute the log-ration of the 2 frequencies.

The result of this script is a csv file containing a table representing the for each pair the corresponding distances(distances going from 0 to 20 ==> 21 lines for the distances: 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)

PS: The distances go until 20, because we are going to need the distance 20 in the prediction script to do the interpolation when the distance is between 19 and 20.

## Script2:Plotting
In this scriptn the interaction profiles are plotted: The scores as a function of the distance.

## Script3:Predecting
This script is used to evaluate RNA conformation. It is partially similar to the first one, as it will compute all the distances for a given structure (same thresholds: 20 Å and i, i+4). For each distance, a scoring value will be computed, using a linear interpolation formula.
By summing all these scores, the script will calculate the estimated Gibbs free energy of the evaluated RNA conformation.

## Code execution
### Script1:
Clone the git in your computer:
```
https://github.com/SafaMADDOURI/RNA_Fold_Problem.git
```
Run the script 1: 
```
$ python3 Training.py <pdb files>
```
The training script creates one csv file, containing the scores at different distances for each pair.
Input the pdb files directly and consicutevely
### Script2: Plotting the scores for each pair
```
$ python3 Plot.py <csv file> -o <Output directory>
```
### Script3: Prediction of Gibbs free energy
```
$ python3 Predict.py <pdb file>
```
## Dependencies
numpy and pandas must be installed to run the scripts
```
pip install numpy
pip install pandas
```
