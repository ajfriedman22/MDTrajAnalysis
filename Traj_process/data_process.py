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

