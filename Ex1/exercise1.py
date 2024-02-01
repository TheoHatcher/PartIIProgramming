"""Exercise1.py - diagonalisation of huckel matrices and displaying the
corresponding energy levels.
"""
import numpy.linalg as linalg
import numpy as np
from platonic_solids import tetrahedron_connectivity_matrix,\
                            hexahedron_connectivity_matrix,\
                            dodecahedron_connectivity_matrix,\
                            buckminster_connectivity_matrix


def calculate_energies(huckel_matrix : np.ndarray) -> None:
    """Diagonalises the matrix given and parses the eigenvalues
    for displaying as energy levels.

    :param huckel_matrix: Huckel matrix being diagonalised and eigenvalues displayed.
    :type huckel_matrix: np.ndarray
    """
    energy_values = solve_huckel_matrix_for_energies(huckel_matrix)
    output_energies(energy_values)



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
        print(f"energy: \u03B1 {energy_level.real:+}\u03B2 ---- g: \
{np.count_nonzero(energies == energy_level)}")
    print("\n")


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





if __name__ == "__main__":
    while True:
        user_input = input("Calculate the energies of (Linear, Cyclic, \
Tetrahedron, Hexahedron, Dodecahedron, Buckminster or exit):")

        match user_input.lower().strip():

            case "linear":
                try:
                    input_number = int(input("Size of linear chain:"))
                except ValueError:
                    print("Error - Invalid Chain Length")
                    continue
                linear_huckel_matrix = construct_linear_polyene_huckel(n=input_number)
                calculate_energies(linear_huckel_matrix)

            case "cyclic":
                try:
                    input_number = int(input("Size of cyclic system:"))
                except ValueError:
                    print("Error - Invalid Cycle Size")
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
                print("Error - Invalid Input")
