#Seperate a string into a list of characters
def split(word):
    return [char for char in word]

#
def convert(s): 
    # initialization of string to "" 
    str1 = "" 
  
    # using join function join the list s by 
    # separating words by str1 
    return(str1.join(s))

#Seperate the numbers in a string
def sep_num(s):
    head = s.rstrip('0123456789')
    tail = s[len(head):]
    return tail

def compute_max(data):
    from scipy.stats import gaussian_kde
    from scipy import optimize
    from scipy.signal import argrelextrema
    import numpy as np

    kde = gaussian_kde(data)
    opt = optimize.minimize_scalar(lambda x: -kde(x))
    samples = np.linspace(min(data), max(data), 50)
    probs = kde.evaluate(samples)
    maxima_index = probs.argmax()
    maxima = samples[maxima_index]
    
    return maxima

