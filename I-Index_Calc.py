import numpy as np
from numpy import linalg as LA

# Define base point and perturbation percent at the top of the script
C10 = 0.004
C01 = 0.001
C11 = 0.001
h   = 10.0

# 1% perturbation
perturbation_percentage = 0.01
base_point = [C10, C01, C11, h]
base_point = np.array(base_point)

N = len(base_point) # Number of parameters
N_Element = int(2*N + 4*N*(N - 1)/2)
# Load W values as text file
file_name = "W_Values.txt"
W = np.loadtxt(file_name)

H = np.zeros((N,N))
# Calculate Hessian Matrix
for i in range(1, N+1):
    # Location for diagonal term
    Entry = (i-1)*2 + (2*N - i)*(i-1)*2
    theta_i = base_point[i-1]
    Delta_i = theta_i*perturbation_percentage
    H[i-1,i-1] = (theta_i*theta_i)*( W[Entry] + W[Entry+1] )/(Delta_i*Delta_i)  
    # Handle off-diagonal terms
    for j in range(i+1, N+1):
        perturb_index = int(Entry + 2 + 4*(j - i - 1))
        theta_j = base_point[j-1]
        Delta_j = theta_j*perturbation_percentage
        H[i-1,j-1] = (theta_i*theta_j/(4.0*Delta_i*Delta_j))*( W[perturb_index] - W[perturb_index+1] - W[perturb_index+2] + W[perturb_index+3] )
        H[j-1,i-1] = H[i-1,j-1]

EigenValues, EigenVector = LA.eig(H);
Lambda_Max = max(EigenValues); Lambda_Min = min(EigenValues);
I_index = np.log10(Lambda_Max/Lambda_Min)
I=np.array([I_index])
print(I_index)
np.savetxt('I_index.txt', I, delimiter=',')