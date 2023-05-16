
from scipy.integrate import quad
import numpy as np
# Define the function to integrate
def f2(x):
    return 1/np.sqrt(x)

# Perform the integration
f2result, f2error = quad(f2, 0, 1)

# Print the result
print("Result:", f2result)
print("Error:", f2error)
