# Introduction
----------------------------
this is a python program to calculate atomic-wavefunction using Restricted Hatree-Fock (RHF) method with Gaussian basis. So it is less accurate when calculating open-shell atoms compairing to close-shell atoms.
# contents
1. [Hatree-Fock Equation](#HFE)
2. [The structure of the program](#SP) 
3. [The output of the program](#OP) 

## <span id="HFE">Hatree-Fock Equation<span>
The most simple atom is Hydrogen! We can easily solve the schrodinger equation of hydrogen to get the exact wavefunction of hydrogen atom.
## <span id="SP">The structure of the program<span>
This program is roughly divided into **Four Parts** which are written as four funtions.

The first one is the *overlap integral* 
```python
def overlap(basis)
```
the result of the function is a matrix named S

The second one is *electrons' knetic energy and colum potential by nucleus* 
```python
def simpleh(basis, znucl)
```
The third one is *interaction between electrons* 
```python
def interaction(basis)
```
The last one is the **self-consistent** function which is used to solve HF equation and get the answer
```python
def scf(basis, znucl, K)
```

## <span id="OP">The output of the program<span>
You will get the eigenvalues and eigenvectors respectively. Using these eigenvectorswith Gausian basis, you can get the wavefunctions of each electron
