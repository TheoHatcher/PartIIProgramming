# Ex 1 : Solving Huckel systems

## Specification
For an arbitrary molecule, specified by its connectivity between adjacent atoms, calculate and print
the Hückel energies and degeneracies of its π-system. You may assume that the only relevant atoms
are all of the same element, are sp2 hybridized, and that adjacent atoms are all the same distance
from each other.

## Run Instructions
N/A - Simply run the file and complete required user input fields.

## Results
Summary of possible outputs from Exercise1.py.
```
Energies of butadiene:
energy: α -1.61803β ---- g: 1
energy: α -0.61803β ---- g: 1
energy: α +0.61803β ---- g: 1
energy: α +1.61803β ---- g: 1


Energies of benzene:
energy: α -2.0β ---- g: 1
energy: α -1.0β ---- g: 2
energy: α +1.0β ---- g: 2
energy: α +2.0β ---- g: 1


Energies of tetrahedral layout:
energy: α -1.0β ---- g: 3
energy: α +3.0β ---- g: 1


Energies of hexahedral layout:
energy: α -3.0β ---- g: 1
energy: α -1.0β ---- g: 3
energy: α +1.0β ---- g: 3
energy: α +3.0β ---- g: 1


Energies of dodecahedron layout:
energy: α -2.23607β ---- g: 3
energy: α -2.0β ---- g: 4
energy: α -0.0β ---- g: 4
energy: α +1.0β ---- g: 5
energy: α +2.23607β ---- g: 3
energy: α +3.0β ---- g: 1


Energies of buckminster fullerene:
energy: α -2.61803β ---- g: 3
energy: α -2.56155β ---- g: 4
energy: α -2.0β ---- g: 4
energy: α -1.61803β ---- g: 5
energy: α -1.43828β ---- g: 3
energy: α -1.30278β ---- g: 5
energy: α -0.38197β ---- g: 3
energy: α -0.13856β ---- g: 3
energy: α +0.61803β ---- g: 5
energy: α +1.0β ---- g: 9
energy: α +1.56155β ---- g: 4
energy: α +1.82025β ---- g: 3
energy: α +2.30278β ---- g: 5
energy: α +2.7566β ---- g: 3
energy: α +3.0β ---- g: 1

```


## Code
The general layout of the code involves taking a user input. Calling the correct function, to generate a connectivity matrix for the specified situation. Then the matrix is diagonalised to find the energies, followed by outputting in a human readable way. The 'platonic_solids.py' contains the connectivity matrices for
the platonic solids, and buckminster fullerene.

'Exercise1.py' contains 5 functions. The code functions by requesting an input from the user. This input is sanitized and the match-case loop handles the different options. Finally if the input matches no valid options, then it repeats notifying the user of the error.
```python
if __name__ == "__main__":
    while True:
        user_input = input("Calculate the energies of (Linear, Cyclic,\
Tetrahedron, Hexahedron, Dodecahedron, Buckminster):")

        match user_input.lower().strip():

            case "linear":
                try:
                    input_number = int(input("Size of linear chain:"))
                except ValueError:
                    continue
                linear_huckel_matrix = construct_linear_polyene_huckel(n=input_number)
                calculate_energies(linear_huckel_matrix)

            case "cyclic":
                try:
                    input_number = int(input("Size of cyclic system:"))
                except ValueError:
                    continue
                cyclic_huckel_matrix = construct_cyclic_polyene_huckel(n=input_number)
                calculate_energies(cyclic_huckel_matrix)

            case "tetrahedron":
                calculate_energies(tetrahedron_connectivity_matrix)
            case "hexahedron":
                calculate_energies(hexahedron_connectivity_matrix)
            case "dodecahedron":
                calculate_energies(dodecahedron_connectivity_matrix)
            case "buckminster":
                calculate_energies(buckminster_connectivity_matrix)
            case "exit":
                exit()
            case _:
                print("Error - Invalid input")
```


The function 'calculate_energies' calls the following two functions, one to calculate the eigenvalues of a connectivity matrix. The second to display them as energy levels.

```python
def solve_huckel_matrix_for_energies(huckel_matrix : np.ndarray) -> np.ndarray[float]:
    """Computes the eigenvalues of a connectivity matrix. Corresponding to the energies
    of a huckel matrix.

    :param huckel_matrix: A connectivity matrix/huckel hamiltonian.
    :type huckel_matrix: np.ndarray
    :return eigenvalues: A list of the eigenvalues of the connectivity matrix. Corresponds
    to the energies in a huckel system.
    :rtype: np.ndarray
    """
    eigenvalues = -np.sort(linalg.eig(huckel_matrix).eigenvalues)[::-1]
    return eigenvalues
```
```python
def output_energies(energies : np.ndarray) -> None:
    r"""Function prints the energies in a human readable output. Adds detail
    of the alpha +/- beta.
    ...
    :param energies: An array containing a list of all the energies.
    :type energies: np.ndarray
    """
    energies = [np.round(energy, 5) for energy in energies]

    print("\n")
    for energy_level in np.unique(energies):
        print(f"energy: \u03B1 {energy_level.real:+}\u03B2 ---- g:\
 {np.count_nonzero(energies == energy_level)}")

    print("\n")
```


The last two functions are responsible for generation of the matrices for the linear and cyclic polyenes.

```python
def construct_linear_polyene_huckel(n : int) -> np.ndarray[float]:
    """Constructs the huckel hamiltonian (equivalent to the connectivity matrix) of
    a linear polyene.

    :param n: Number of carbons/interacting pi orbitals in the linear polyene
    :type n: int
    :param print_energies: controls whether the energies are printed, defaults to True
    :type print_energies: bool, optional
    :return: A tuple containing the matrix and energies of the huckel system.
    :rtype: np.ndarray[float]
    """
    assert n >= 2

    matrix = np.diagflat([-1 for i in range(n-1)], -1) +\
             np.diagflat([-1 for i in range(n-1)],  1)

    return matrix
```

```python
def construct_cyclic_polyene_huckel(n : int) -> np.ndarray[float]:
    """Constructs the huckel hamiltonian (equivalent to the connectivity matrix) of
    a cyclic polyene.

    :param n: Number of carbons/interacting pi orbitals in the cyclic polyene
    :type n: int
    :param print_energies: controls whether the energies are printed, defaults to True
    :type print_energies: bool, optional
    :return: A tuple containing the matrix and energies of the huckel system.
    :rtype: np.ndarray[float]
    """
    assert n >= 3

    matrix = np.diagflat([-1 for i in range(n-1)], -1) +\
             np.diagflat([-1 for i in range(n-1)],  1)

    matrix[n-1, 0] = matrix[0, n-1] = -1

    return matrix
```
