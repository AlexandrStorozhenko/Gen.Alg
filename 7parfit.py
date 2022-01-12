import numpy as np
from math import e

# System best vars
bestVal = 10**82 #approx. number of atoms in the obs.univ.
best_parameter_fit = [0,0]
best_model_S = []
best_model_I = []
best_model = [best_model_S,best_model_I]

# Parameter Fit Data
y_obs = [[0.02518109693,1-0.02518109693],[0.07726802346,1-0.07726802346],[0.2590548465, 1-0.2590548465],[0.263884098,1-0.263884098], [0.4101414281, 1-0.4101414281], [1.,0.],[0.2404277337, 1-0.2404277337],
[0.2307692308, 1-0.2307692308],[0.1193515005,1-0.1193515005],[0.1193515005, 1-0.1193515005],[0.130044843,1-0.130044843],[0.02311141773,1-0.02311141773],[0.3204553294, 1-0.3204553294],[0.07933770266, 1-0.07933770266],
[0.1545360469, 1-0.1545360469],[0.0217316316, 1-0.0217316316]]

model_S = [0.97481890307]
model_I = [0.02518109693]
temp_model_S = [0.97481890307]
temp_model_I = [0.02518109693]

def EulerMethod(Si, Ii, alpha, k, distf, Max, a, b, c, dt):
    model_S = [0.97481890307]
    model_I = [0.02518109693]
    model_D = [distf]
    model_R = [0]
    model_w = [0]
    temp_model_S = [0.97481890307]
    temp_model_I = [0.02518109693]
    temp_model_D = [distf]
    temp_model_R = [0]
    temp_model_w = [0]
    
    for i in range(0, 16*int(1/dt)):
        temp_model_D.append(temp_model_D[i]+(k*temp_model_D[i]*(1-temp_model_D[i]/Max)*dt))
        temp_model_w.append(temp_model_w[i]+((a*(i-b)*e**(-((i-b)**2)/2*(c**2)))/(c**2))*dt)
        temp_model_S.append(temp_model_S[i]-(alpha*temp_model_I[i]*temp_model_S[i]-temp_model_w[i]*temp_model_R[i])*dt)
        temp_model_I.append(temp_model_I[i]+(alpha*temp_model_I[i]*temp_model_S[i]-temp_model_D[i]*temp_model_I[i])*dt)
        temp_model_R.append(temp_model_R[i]+(temp_model_D[i]*temp_model_I[i]-temp_model_w[i]*temp_model_R[i])*dt)

        if i % int(1/dt) == 0:
            model_w.append(temp_model_w[i]+((a*(i-b)*e**(-((i-b)**2)/2*(c**2)))/(c**2))*dt)
            model_D.append(temp_model_D[i]+(k*temp_model_D[i]*(1-temp_model_D[i]/Max)*dt))
            model_S.append(temp_model_S[i]-(alpha*temp_model_I[i]*temp_model_S[i]-temp_model_w[i]*temp_model_R[i])*dt)
            model_I.append(temp_model_I[i]+(alpha*temp_model_I[i]*temp_model_S[i]-temp_model_D[i]*temp_model_I[i])*dt)
            model_R.append(temp_model_R[i]+(temp_model_D[i]*temp_model_I[i]-temp_model_w[i]*temp_model_R[i])*dt)

    return [model_S, model_I]

# Numerical solution

#Implement Genetic Algorithm
for alpha in range(10, 20):
    alpha = alpha/10
    for k in range(1, 10):
        k = k/10
        for distf in range(1, 10):
            distf = distf/100
            for Max in range(1, 10):
                Max = Max/10
                for a in range(1, 10):
                    a = a/10
                    for b in range(1, 10):
                        b = b/10
                        for c in range(1, 10):
                            c = c/10
                            y_model = EulerMethod(model_S[0], model_I[0], alpha, k, distf, Max, a, b, c, 0.01)
                                        
                            # RSS (Residual Sum of Squares)
                            rssVal = 0

                            for i in range(0,len(y_model[0])-1):
                                rssVal += (y_obs[i][0] - y_model[1][i])**2

                            if(rssVal < bestVal):
                                bestVal = rssVal
                                best_model = y_model
                                best_parameter_fit = [alpha, k, distf, Max, a, b, c]
                                print(best_parameter_fit)
        
print("Best Parameter Fit Found")
print(best_model[1])
print("Parameters")
print(best_parameter_fit)
print("RSS value")
print(bestVal)