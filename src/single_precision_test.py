import numpy as np

# Answer to problem 1a

machinePrecision = np.finfo(np.single).eps
print(f"The machine precision value for single precision is {machinePrecision}")