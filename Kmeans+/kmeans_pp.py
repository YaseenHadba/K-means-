import pandas as pd
import numpy as np
import mykmeanssp as kmp
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("k", type=int)
parser.add_argument("max_iter", type=int, nargs='?', default=300, const=0)
parser.add_argument("epsilon", type=float)
parser.add_argument("file_1", type=str)
parser.add_argument("file_2", type=str)
args = parser.parse_args()

def calc_norm(centroids, point):  
    closest_norm = 0
    closest_norm_value = diff_norm(centroids[0], point)
    for i in range(len(centroids)):
        x = diff_norm(centroids[i], point)
        if x < closest_norm_value:
            closest_norm = i
            closest_norm_value = x
    return closest_norm


def help2(vector1, vector2): #help func to for calculating the norm 
    arr = [0 for i in range(len(vector2))]
    for i in range(len(vector2)):
        arr[i] = vector1[i] - vector2[i]
    sum2 = 0
    for i in range(len(vector2)):
        sum2 += pow(arr[i], 2)
    return sum2


def sum(list1, list2): 
    result = []
    for i in range(len(list1)):
        result.append(list1[i] + list2[i])
    return result


def divide(vector, num):  
    for i in range(len(vector)):
        vector[i] = vector[i] / num
    return vector

def starting_centroids(k, vectors):
    N = len(vectors)
    i=0
    indexes = [0 for i in range(k)]
    np.random.seed(0)
    indexes[0] = np.random.choice(N)
    mu0 = vectors[indexes[0]]
    mu = np.zeros((k, len(vectors[0])))
    mu[0] = mu0
    while i < k:
        i+=1
        if(i==k) :break
        D = [0 for s in range(N)]
        for l in range(0,N):
            minimum = help2(vectors[l], mu[0])
            for j in range(0,i):
                calculation = help2(vectors[l], mu[j])
                if calculation < minimum:
                    minimum = calculation
            D[l] = minimum
        sum2 = 0
        for m in range(0,N):
            sum2 = sum2 + D[m]
        p = [0 for i in range(N)]
        for l in range(0,N):
            p[l] = D[l] / sum2
        indexes[i] = np.random.choice(N, p=p)
        mu[i] = vectors[indexes[i]]
    return mu, indexes


def merge(input_file1, input_file2):
    input1 = pd.read_csv(input_file1, header=None)
    input2 = pd.read_csv(input_file2, header=None)
    vectors = pd.merge(input1, input2, on=0)
    vectors.sort_values(by=[0], inplace=True, ascending=True)
    vectors = vectors.drop(columns=0)
    return vectors

def centroids_result(k, max_iter, epsilon, d,  vectors, starting):    
    starting = starting.tolist()
    centroids = kmp.fit(k, max_iter, epsilon, d, vectors, starting)   
    for i in range(len(centroids)):
        print(','.join([format(centroids[i][j], ".4f") for j in range(len(centroids[i]))]))
    print('\n')

def execute(k, max_iter, epsilon,vectors):
    starting = starting_centroids(k, vectors) 
    centroids = starting[0]
    d = len(centroids[0])
    indexes = starting[1]
    print(','.join([str(indexes[j]) for j in range(len(indexes))]))
    centroids_result(k, max_iter, epsilon, d, vectors, centroids)



def valid(input):
    if not (str(input[0]).isdigit()) or (not (str(input[1]).isdigit())):
        return 0
    if(input[0]<0 or input[1]<0): return 0
    return 1


input_args = [args.k, args.max_iter, args.epsilon, args.file_1, args.file_2]
file1 = merge(args.file_1, args.file_2)
vectors = file1.values.tolist()
if(len(vectors) < args.k) :
    print("Invalid Input!")
if(not valid(input_args)):
    print("Invalid Input!")
execute(input_args[0], input_args[1], input_args[2], vectors)
