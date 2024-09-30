import numpy as np

# Define base point and perturbation percent at the top of the script
C10 = 0.004
C01 = 0.001
C11 = 0.001
h   = 10.0

# 1% perturbation
perturbation_percentage = 0.01
base_point = [C10, C01, C11, h]
base_point = np.array(base_point)


def generate_doe_matrix(base_point, perturbation_percentage):
    """
    Generates the DOE matrix for N parameters.
    
    Parameters:
    base_point (list): The base point (theta values) for calculation, a list of N parameters.
    perturbation_percentage (float): The percentage by which each parameter will be perturbed (e.g., 0.01 for 1%).
    
    Returns:
    doe_matrix (list of lists): DOE matrix containing perturbed parameter sets.
    """
    N = len(base_point)  # Number of parameters
    N_Element = int(2*N + 4*N*(N - 1)/2)

    # Create DOE matrix with perturbations
    doe_matrix = np.ones((N_Element, N))

    for i in range(1, N+1):
        # Location for diagonal term
        Entry = (i-1)*2 + (2*N - i)*(i-1)*2
        
        # Handle diagonal terms
        doe_matrix[Entry, i-1] = 1.0 + perturbation_percentage
        doe_matrix[Entry+1, i-1] = 1.0 - perturbation_percentage
        
        # Handle off-diagonal terms
        for j in range(i+1, N+1):
            perturb_index = int(Entry + 2 + 4*(j - i - 1))
            
            if perturb_index + 3 >= N_Element:
                raise IndexError(f"Calculated perturb_index {perturb_index} exceeds DOE matrix size {N_Element}")
            
            doe_matrix[perturb_index, i-1] = 1.0 + perturbation_percentage
            doe_matrix[perturb_index, j-1] = 1.0 + perturbation_percentage
            
            doe_matrix[perturb_index + 1, i-1] = 1.0 + perturbation_percentage
            doe_matrix[perturb_index + 1, j-1] = 1.0 - perturbation_percentage
            
            doe_matrix[perturb_index + 2, i-1] = 1.0 - perturbation_percentage
            doe_matrix[perturb_index + 2, j-1] = 1.0 + perturbation_percentage
            
            doe_matrix[perturb_index + 3, i-1] = 1.0 - perturbation_percentage
            doe_matrix[perturb_index + 3, j-1] = 1.0 - perturbation_percentage
    # Multiply base point value
    

    return doe_matrix

if __name__ == "__main__":
    # Generate DOE matrix
    doe_matrix = generate_doe_matrix(base_point, perturbation_percentage)
    doe_matrix = doe_matrix * base_point
    
    # Save DOE matrix to a text file
    file_name = "DOE_matrix.txt"  # Specify the output file name
    np.savetxt(file_name, doe_matrix, fmt='%.6e', delimiter=",")
