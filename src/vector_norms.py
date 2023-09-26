import numpy as np

# Calculates the two-norm/Euclidean length of a vector (problem 2a)
def euclideanNorm(vector):
    lengthSquared = 0
    for x in vector:
        lengthSquared += x * x
    return np.sqrt(lengthSquared)

# Calculates the one-norm length of a vector (problem 2b)
def oneNorm(vector):
    length = 0
    for x in vector:
        length += np.abs(x)
    return length
    
# Calculates the infinity norm/sup-norm length of a vector (problem 2c)
def supNorm(vector):
    length = 0
    for x in vector:
        if (np.abs(x) > length):
            length = np.abs(x)
    return length
    
# Calculates the two-norm/Euclidean distance between two vectors (problem 3a)
def euclideanDistance(vector1, vector2):
    # The two vectors must be the same length
    if len(vector1) != len(vector2):
        raise ValueError(f"Vectors are different lengths: {len(vector1)} and {len(vector2)}")
    
    distanceSquared = 0
    for i in range(len(vector1)):
        x = vector1[i] - vector2[i]
        distanceSquared += x * x
    return np.sqrt(distanceSquared)
    
# Calculates the one-norm distance between two vectors (problem 3b)
def oneNormDistance(vector1, vector2):
    # The two vectors must be the same length
    if len(vector1) != len(vector2):
        raise ValueError(f"Vectors are different lengths: {len(vector1)} and {len(vector2)}")
    
    distance = 0
    for i in range(len(vector1)):
        x = vector1[i] - vector2[i]
        distance += np.abs(x)
    return distance

# Calculates the infinity-nrom distance between two vectors (problem 3c)
def supNormDistance(vector1, vector2):
    # The two vectors must be the same length
    if len(vector1) != len(vector2):
        raise ValueError(f"Vectors are different lengths: {len(vector1)} and {len(vector2)}")
    
    distance = 0
    for i in range(len(vector1)):
        x = vector1[i] - vector2[i]
        if (np.abs(x) > distance):
            distance = np.abs(x)
    return distance