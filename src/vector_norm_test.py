# Test program for vector norms using randomly generated values
# (Not part of the assignment)

from vector_norms import *

# Functions used for generating test vectors
def randomVector(numDimensions):
    return [np.random.randint(-10, 10) for x in range(numDimensions)]

# Test norms
print("Testing norms:")

numDimensions = np.random.randint(1, 5)
vector = randomVector(numDimensions)
length = euclideanNorm(vector)
print(f"\tThe Euclidean norm of {vector} is {length}")

numDimensions = np.random.randint(1, 5)
vector = randomVector(numDimensions)
length = oneNorm(vector)
print(f"\tThe one-norm of {vector} is {length}")

numDimensions = np.random.randint(1, 5)
vector = randomVector(numDimensions)
length = supNorm(vector)
print(f"\tThe infinity norm of {vector} is {length}")

# Test distances
print("\nTesting distances:")

numDimensions = np.random.randint(1, 5)
vector1 = randomVector(numDimensions)
vector2 = randomVector(numDimensions)
distance = euclideanDistance(vector1, vector2)
print(f"\tThe Euclidean distance between {vector1} and {vector2} is {distance}")

numDimensions = np.random.randint(1, 5)
vector1 = randomVector(numDimensions)
vector2 = randomVector(numDimensions)
distance = oneNormDistance(vector1, vector2)
print(f"\tThe one-norm distance between {vector1} and {vector2} is {distance}")

numDimensions = np.random.randint(1, 5)
vector1 = randomVector(numDimensions)
vector2 = randomVector(numDimensions)
distance = supNormDistance(vector1, vector2)
print(f"\tThe infinity-norm distance between {vector1} and {vector2} is {distance}")