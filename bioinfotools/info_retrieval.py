from bioinfotools.utils import BASES

def count_bases(text, to_dict=False, alphabet=None):
    if alphabet is None:
        alphabet = BASES
    
    counter = [0 for _ in range(len(alphabet))]

    for char in text:
        for i, symbol in enumerate(alphabet):
            if char == symbol:
                counter[i] += 1
                break

    if to_dict:
        counter = {b: c for b, c in (alphabet, counter)}
    
    return counter

def gc_skews(text):
    skews = [0]
    cur_skew = 0

    for char in text:
        if char == "G":
            cur_skew += 1
        elif char == "C":
            cur_skew -= 1
        
        skews.append(cur_skew)
    
    return skews

def minimum_gc_skew(text):
    # Init value of 2 is to ensure we also process sequences with a positive minima. 
    # We could initialize it to a larger number, but given that each iteration can update the value of the skew by at most 1, this is a sufficiently large value.
    min_skew = 2
    min_loc = []
    cur_skew = 0

    for i, char in enumerate(text, start=1):
        if char == "G":
            cur_skew += 1
        elif char == "C":
            cur_skew -= 1
        
        if cur_skew == min_skew:
            min_loc.append(i)
        elif cur_skew < min_skew:
            min_skew = cur_skew
            min_loc = [i]
    
    return min_loc
