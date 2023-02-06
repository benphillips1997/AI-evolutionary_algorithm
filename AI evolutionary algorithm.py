import random as rand
import copy
import matplotlib.pyplot as plt
import numpy as np
import csv


N = 20   # amount of genes for individual
P = 50   # amount of population of individuals
L = 10   # times to loop for accuracy
G = 50   # amount of generations
MUTRATE = 0.025  # mutation rate
MUTSTEP = 0.75  # amount to change gene by if mutation occurs

file = "datafile.csv"   # file to save data in
show_graph = True       # whether or not to show graph
gen_start = 0           # if graph is shown, which generation to start displaying the results

# func1
#GENE_MIN = -100
#GENE_MAX = 100

# func2
GENE_MIN = -5
GENE_MAX = 10


class Individual:
    def __init__(self):
        self.gene = [0]*N
        self.fitness = 0



# func 1
"""
def calculate_fitness(ind):
    fitness = 0
    for i in range(1, N-1):
        fitness = fitness + (100 * (ind.gene[i+1] - (ind.gene[i] ** 2)) ** 2) + ((1 - ind.gene[i]) ** 2)
    ind.fitness = fitness
    #print(str(fitness))
    return fitness
"""

# func 2
#"""
def calculate_fitness(ind):
    fitness = 0

    sum1 = 0
    for i in range(1, N):
        sum1 = sum1 + (ind.gene[i] ** 2)

    sum2 = 0
    for j in range(1, N):
        sum2 = sum2 + (0.5 * j * ind.gene[j])
    sum2 = sum2 ** 2

    sum3 = 0
    for k in range(1, N):
        sum3 = sum3 + (0.5 * k * ind.gene[k])
    sum3 = sum3 ** 4

    fitness = sum1 + sum2 + sum3
    ind.fitness = fitness
    #print(str(fitness))
    return fitness
#"""


def get_lowest_fitness(pop):
    lowest_fitness = pop[0].fitness
    for ind in pop:
        if ind.fitness < lowest_fitness:
            lowest_fitness = ind.fitness
    return lowest_fitness


def get_mean_fitness(pop):
    mean_fitness = 0
    for individual in pop:
        mean_fitness = mean_fitness + individual.fitness
    mean_fitness = mean_fitness / len(pop)
    return mean_fitness



total_lowest_fitnesses = np.zeros(G)
total_mean_fitnesses = np.zeros(G)

# loop algorithm for accuracy
for l in range(L):
    lowest_fitnesses = []
    mean_fitnesses = []

    individuals_array = []

    for x in range(P):
        # create array N long and populate with random floats between 1 and 0
        tempgene = []
        for y in range(N):
            tempgene.append(rand.uniform(GENE_MIN, GENE_MAX))
        individual = Individual()
        individual.gene = tempgene.copy()
        calculate_fitness(individual)
        individuals_array.append(individual)

    # run generations
    for gen in range(G):
        # get offspring from individuals array
        offspring_array = []
        for i in range(P):
            parent1 = rand.randint(0, P-1)
            off1 = copy.deepcopy(individuals_array[parent1])
            parent2 = rand.randint(0, P-1)
            off2 = copy.deepcopy(individuals_array[parent2])
            # add the better comparison to offspring
            if off1.fitness < off2.fitness:
                offspring_array.append(off1)
            else:
                offspring_array.append(off2)


        mutated_array = []
        # cross over
        toff1 = Individual()
        toff2 = Individual()
        toff1_copy = Individual()
        for i in range(0, P, 2):
            toff1 = copy.deepcopy(offspring_array[i])
            toff2 = copy.deepcopy(offspring_array[i+1])
            toff1_copy = copy.deepcopy(offspring_array[i])
            cross_over_point = rand.randint(1, N)
            # cross over genes starting from cross over point to the end
            for j in range(cross_over_point, N):
                toff1.gene[j] = toff2.gene[j]
                toff2.gene[j] = toff1_copy.gene[j]
            # add crossed over offspring to new array
            mutated_array.append(copy.deepcopy(toff1))
            mutated_array.append(copy.deepcopy(toff2))

        # mutation
        for i in range(P):
            new_individual = Individual()
            new_individual.gene = []
            for j in range(N):
                gene = mutated_array[i].gene[j]
                mutprob = rand.random()
                if mutprob < MUTRATE:
                    alter = rand.uniform(-MUTSTEP, MUTSTEP)
                    gene = gene + alter
                    if gene < GENE_MIN:
                        gene = GENE_MIN
                    if gene > GENE_MAX:
                        gene = GENE_MAX
                new_individual.gene.append(gene)
            calculate_fitness(new_individual)
            mutated_array[i] = copy.deepcopy(new_individual)

        # calculate lowest and mean fitness
        lowest_fitness = get_lowest_fitness(mutated_array)
        mean_fitness = get_mean_fitness(mutated_array)

        if gen > 0:
            if lowest_fitness > lowest_fitnesses[gen-1]:
                lowest_fitness = lowest_fitnesses[gen-1]

        # save lowest and mean fitness from this generation
        lowest_fitnesses.append(lowest_fitness)
        mean_fitnesses.append(mean_fitness)

        # copy mutated population into individuals array to begin new generation
        individuals_array = copy.deepcopy(mutated_array)

        # for following what point the algorithm is at, used for when it took a long time to run
        '''if (gen+1) % 100 == 0:
            print(f"Loop {l+1} gen {gen+1}")'''

    #add up total fitnesses
    for i in range(G):
        total_lowest_fitnesses[i] = total_lowest_fitnesses[i] + lowest_fitnesses[i]
        total_mean_fitnesses[i] = total_mean_fitnesses[i] + mean_fitnesses[i]


# divide by the amount of times looped for average
for i in range(G):
    total_lowest_fitnesses[i] = total_lowest_fitnesses[i] / L
    total_mean_fitnesses[i] = total_mean_fitnesses[i] / L


# print data per generation in terminal
'''for i in range(G):
    print("Gen " + str(i+1))
    print("lowest fitness = " + str('%.10E' % total_lowest_fitnesses[i]))
    print("mean fitness = " + str('%.10E' % total_mean_fitnesses[i]) + "\n")'''


# saves data per generation in file
with open(file, 'w') as f:
    writer = csv.writer(f)
    for i in range(G):
        writer.writerow([i+1, ' %.3E' % total_lowest_fitnesses[i], ' %.3E' % total_mean_fitnesses[i]])


if show_graph:
    plt.ylabel("fitness")
    plt.xlabel("generation")
    plt.title(f"N = {N}, P = {P}, G = {G}, MUTRATE = {MUTRATE}, MUTSTEP = {MUTSTEP}")
    plt.plot(np.arange(gen_start, G), total_lowest_fitnesses[gen_start:])
    plt.plot(np.arange(gen_start, G), total_mean_fitnesses[gen_start:])
    plt.show()