from bioinfotools.suffix_tree import SuffixTree

#@sartor-bot
def compute_z(text):
    """
    Computes the Z-array of a text.

    Parameters:
        text (str): Input string.

    Returns:
        list: An array where each index `i` contains the length of the longest
              substring starting from `text[i]` that is also a prefix of `text`.
    """
    n = len(text)
    z = [0] * n
    left, right = 0, 0

    for i in range(1, n):
        if i <= right:
            z[i] = min(right - i + 1, z[i - left])

        while i + z[i] < n and text[z[i]] == text[i + z[i]]:
            z[i] += 1

        if i + z[i] - 1 > right:
            left, right = i, i + z[i] - 1

    return z

def compute_lps(text):
    """
    Computes the LPS (Longest Prefix Suffix) array for a text.

    Parameters:
        text (str): The input string to compute the LPS array for.

    Returns:
        list: The LPS array for `text`.
    """
    n = len(text)
    lps = [0] * n
    length = 0
    i = 1

    while i < n:
        if text[i] == text[length]:
            length += 1
            lps[i] = length
            i += 1
        else:
            if length != 0:
                length = lps[length - 1]
            else:
                lps[i] = 0
                i += 1

    return lps


def encode_bwt(text, end_char='$'):
    """
    Encode a string into its Burrows-Wheeler Transform (BWT).

    Parameters:
        text (str): Text to encode.
        end_char (char): End character in the text. Default to `'$'`.
    """
    if end_char not in text:
        text += end_char

    rotations = [text[i:] + text[:i] for i in range(len(text))]
    sorted_rotations = sorted(rotations)

    bwt = ''.join(row[-1] for row in sorted_rotations)

    return bwt

def decode_bwt(bwt, end_char='$'):
    """
    Decodes a Burrows-Wheeler Transform (BWT) into its original text.

    Parameters:
        bwt (str): Transform to decode.
        end_char (char): End character in the text. Default to `'$'`.
    """
    if not len(bwt):
        raise ValueError("Invalid BWT string of length 0.")
    if end_char not in bwt:
        raise ValueError("End Character not in BWT string.")
    
    sorted_text = sorted(bwt, key=lambda c: (c != end_char, c))
    alphabet = set(sorted_text)
    position_chars = {}
    letters_count = {char: 0 for char in alphabet}
    inverse_map = []

    for i, char in enumerate(sorted_text):
        if char not in position_chars:
            position_chars[char] = i

    for c in bwt:
        inverse_map.append(position_chars[c] + letters_count[c])
        letters_count[c] += 1

    idx = 0
    res = [end_char]

    while True:
        if bwt[idx] == end_char:
            break
        res.append(bwt[idx])
        idx = inverse_map[idx]

    text = ''.join(reversed(res))

    return text

def build_suffix_tree(text: str) -> SuffixTree:
    tree = SuffixTree(text)
    return tree
