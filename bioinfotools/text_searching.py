from bioinfotools.utils import BASES
from bioinfotools.data_compression import compute_lps, compute_z

def naive_search(pattern, text):
    """
    Performs a naive pattern matching algorithm.

    This method checks every possible substring of `text` that matches 
    the length of `pattern` and compares it to `pattern`. It has a 
    worst-case time complexity of O(n * m), where n is the length of 
    `text` and m is the length of `pattern`.

    Parameters:
        pattern (str): The pattern to search for.
        text (str): The text in which to search for the pattern.

    Returns:
        list: A list of starting indices where `pattern` occurs in `text`.
    """
    k = len(pattern)
    l = len(text)
    result = []

    for i in range(l):
        substr = text[i:i+k]
        if pattern == substr:
            result.append(i)
    
    return result


def z_algorithm(pattern, text, separator='$'):
    """
    Performs pattern matching using the Z-Algorithm.

    The Z-Algorithm preprocesses the concatenated string (pattern + separator + text)
    and constructs a Z-array, which is then used to efficiently locate occurrences of 
    `pattern` in `text`. The Z-array stores the length of the longest substring 
    starting from each position that matches the prefix of the concatenated string.

    Parameters:
        pattern (str): The pattern to search for.
        text (str): The text in which to search for the pattern.
        separator (str): A unique separator character used to concatenate `pattern`
                         and `text`. Default is `'$'`.

    Returns:
        list: A list of starting indices where `pattern` occurs in `text`.
    """
    pattern_size = len(pattern)
    concat_string = pattern + separator + text
    z_array = compute_z(concat_string)
    result = []

    for i, z in enumerate(z_array[pattern_size+1:]):
        if z == pattern_size:
            result.append(i)

    return result


def kmp_search(pattern, text):
    """
    Performs the Knuth-Morris-Pratt (KMP) pattern matching algorithm.

    Parameters:
        text (str): The text to search within.
        pattern (str): The pattern to search for in the text.

    Returns:
        list: A list of starting indices where the pattern is found in the text.
    """
    m, n = len(pattern), len(text)
    lps = compute_lps(pattern)
    i, j = 0, 0
    result = []

    while i < n:
        if text[i] == pattern[j]:
            i += 1
            j += 1
        
        if j == m:
            result.append(i - j)
            j = lps[j - 1]
        
        elif i < n and text[i] != pattern[j]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1

    return result


def _init_profile(k, alphabet_size):
    profile = [[1 for _ in range(k)] for _ in range(alphabet_size)]
    return profile

def _update_profile(profile, kmer, alpha_map):
    for i, char in enumerate(kmer):
        profile[alpha_map[char]][i] += 1

    return profile

def _compute_val(kmer, profile, alpha_map):
    val = 1
    for i, char in enumerate(kmer):
        val *= profile[alpha_map[char]][i]

    return val

def _get_most_probable(string, profile, alpha_map):
    s_len = len(string)
    k = len(profile[0])
    max_val = 0
    max_k = ""
    for i in range(s_len-k+1):
        cur_k = string[i:i+k]
        if _compute_val(cur_k, profile, alpha_map) > max_val:
            max_val = _compute_val(cur_k, profile, alpha_map)
            max_k = cur_k
    return max_k

def score(motifs, alpha_map):
    score = 0
    k = len(motifs[0])
    t = len(motifs)
    for i in range(k):
        counter = [0 for _ in range(4)]
        for j in range(t):
            counter[alpha_map[motifs[j][i]]] += 1
        score += t - max(counter)
    return score

def greedy_motif_search(k, t, sequences, alphabet=None):
    if alphabet is None:
        alphabet = BASES

    alpha_map = {l: n for n, l in enumerate(alphabet)}
    best_motifs = [seq[:k] for seq in sequences]
    s1_len = len(sequences[0])
    for i in range(s1_len-k+1):
        profile = _init_profile(k, len(alphabet))
        motif_1 = sequences[0][i:i+k]
        profile = _update_profile(profile, motif_1, alpha_map)
        motifs = [motif_1]
        for j in range(1, t):
            motif_j = _get_most_probable(sequences[j], profile, alpha_map)
            profile = _update_profile(profile, motif_j, alpha_map)
            motifs.append(motif_j)
        if score(motifs, alpha_map) < score(best_motifs, alpha_map):
            best_motifs = motifs

    return best_motifs
