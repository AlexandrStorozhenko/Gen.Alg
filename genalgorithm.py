import math
from math import e
import random
from math import log10, floor

def round_sig(x, sig=8):
    if(x==0):
        return 0
    else:
        return round(x, sig-int(floor(log10(abs(x))))-1)

# System best vars
bestVal = 10**10
best_parameter_fit = [0,0]
best_model_S = []
best_model_I = []
best_model = [best_model_S,best_model_I]

population = []
newPopulation = []
count = 0

# Parameter Fit Data
y_obs = [[0.02518109693,1-0.02518109693],[0.07726802346,1-0.07726802346],[0.2590548465, 1-0.2590548465],[0.263884098,1-0.263884098], [0.4101414281, 1-0.4101414281], [1.,0.],[0.2404277337, 1-0.2404277337],
[0.2307692308, 1-0.2307692308],[0.1193515005,1-0.1193515005],[0.1193515005, 1-0.1193515005],[0.130044843,1-0.130044843],[0.02311141773,1-0.02311141773],[0.3204553294, 1-0.3204553294],[0.07933770266, 1-0.07933770266],
[0.1545360469, 1-0.1545360469],[0.0217316316, 1-0.0217316316]]

def logisticFunc(Max, distf, k,t):
    val = Max/(1+(e**(-(distf*(t-k)))))
    if(math.isinf(val)):
        print("wtf")
        print("break condition")

    return val

def gaussFunc(a,b,c,t):
    val = a*(e**((-((t-b)**2)/(2*(c**2)))))
    if(math.isinf(val)):
        print("wtf")
        print("break condition")

    return val

def EulerMethod(Si, Ii, alpha, Max, distf, k, a, b, c, dt):
    model_S = [0.97481890307]
    model_I = [0.02518109693]
    model_R = [0]
    temp_model_S = [0.97481890307]
    temp_model_I = [0.02518109693]
    temp_model_R = [0]
    inf_param = alpha
    distf_param = 0
    recov_param = 0
    currT = 0

    for i in range(0, 16*int(1/dt)):
        recov_param = gaussFunc(Max,distf,k,currT)
        distf_param = logisticFunc(a,b,c, currT)
        temp_model_S.append(round_sig(temp_model_S[i]-(alpha*temp_model_I[i]*temp_model_S[i]-recov_param*temp_model_R[i])*dt))
        temp_model_I.append(round_sig(temp_model_I[i]+(alpha*temp_model_I[i]*temp_model_S[i]-distf_param*temp_model_I[i])*dt))
        temp_model_R.append(round_sig(temp_model_R[i]+(distf_param*temp_model_I[i]-recov_param*temp_model_R[i])*dt))

        if i % int(1/dt) == 0:
            currT += 1
            model_S.append(round_sig(temp_model_S[i]-(alpha*temp_model_I[i]*temp_model_S[i]-recov_param*temp_model_R[i])*dt))
            model_I.append(round_sig(temp_model_I[i]+(alpha*temp_model_I[i]*temp_model_S[i]-distf_param*temp_model_I[i])*dt))
            model_R.append(round_sig(temp_model_R[i]+(distf_param*temp_model_I[i]-recov_param*temp_model_R[i])*dt))

    return [model_S, model_I]

def fitness(params):
    y_model = EulerMethod(y_obs[0][0], y_obs[0][1], params[0], params[1], params[2], params[3], params[4], params[5], params[6], 0.005)
                                        
    # RSS (Residual Sum of Squares)
    rssVal = 0

    for i in range(0,len(y_model[0])-1):
        rssVal += (y_obs[i][0] - y_model[1][i])**2

    global bestVal
    global best_model
    global best_parameter_fit
    if(rssVal < bestVal):
        bestVal = rssVal
        best_model = y_model
        best_parameter_fit = []
        for i in range(0, len(params)):
            best_parameter_fit.append(params[i])
        print("Best Current Parameter Fit Found")
        print(best_model[1])
        print("Parameters")
        print(best_parameter_fit)
        print("RSS value")
        print(bestVal)

    return rssVal

def selection():
    wheelSelection = []
    sum = 0
    onesum = 0
    for i in range(0,len(population)):
        fitn = fitness(population[i])
        wheelSelection.append(fitn)
        sum += fitn
    for i in range(0,len(wheelSelection)):
        wheelSelection[i] = (wheelSelection[i]/sum)
        onesum += (wheelSelection[i])
    return wheelSelection

def reproduction():
    global population
    global newPopulation
    global count
    select = selection()
    count+=1
    for a in range(0,len(population)-1):
        # 7 crossover populations
        if a<7:
            # roulette wheel seelction (2 random numbers prob1/2 are generated --> they are always less than 1, the probability of each population is then substracted, whichever one it stops on, the number of the population is assigned through i and j) 
            i = 0
            j = 0
            while i == j:
                i = 0
                j = 0
                prob1 = random.random()
                i=0
                while prob1 > 0:
                    prob1 -= select[i]
                    i += 1
                prob2 = random.random()
                j=0
                while prob2 > 0:
                    prob2 -= select[j]
                    j += 1
            i-=1
            j-=1
            n = 0
            k=0
            while n==k:
                n = random.randint(0,7)
                k = random.randint(0,7)
                if k<n:
                    temp = k
                    k=n
                    n=k
            newPopulation.append(population[(i)][0:n]+population[(j)][n:k]+population[i][k:7])
        # 3 new random populations
        else:
            for a in range(7, 10):
                paramset = []
                for i in range (0,7):
                    param_val = random.randint(1, 10000)/1000
                    paramset.append(param_val)
                newPopulation.append(paramset)

      
        mutation = random.randrange(0,2000)/1000
        index = random.randrange(0,6)
        signRand = random.random()
        if(signRand>=0.5):
            newPopulation[a][index] += mutation
        elif(newPopulation[a][index]>mutation):
            newPopulation[a][index] -= mutation

        population[a] = newPopulation[a]

    newPopulation = []

def geneticAlgorithm():
    for a in range(0, 10):
        paramset = []
        for i in range (0,7):
            param_val = random.randint(1, 10000)/1000
            paramset.append(param_val)
        
        population.append(paramset)

    while bestVal > 0.1:
        reproduction()

    return best_parameter_fit

solution = geneticAlgorithm()
print("Best Parameter Fit Found")
print(best_model[1])
print("Parameters")
print(solution)
print("RSS value")
print(bestVal)