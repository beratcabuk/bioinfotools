import pytest
from bioinfotools.data_compression import compute_z

def test_compute_z_simple():
    """
    Test compute_z on a simple string with some repeated prefixes.
    """
    text = "aabcaabxaaaz"
    # Z-array for 'aabcaabxaaaz' is [0,1,0,0,3,1,0,0,2,2,1,0]
    expected = [0,1,0,0,3,1,0,0,2,2,1,0]
    assert compute_z(text) == expected

def test_compute_z_no_repeats():
    """
    Test compute_z on a string where all characters are unique.
    """
    text = "abcdefg"
    # Each position except for 0 should be 0, since no repeats.
    expected = [0,0,0,0,0,0,0]
    assert compute_z(text) == expected

def test_compute_z_all_same_char():
    """
    Test compute_z on a string with all the same characters.
    """
    text = "aaaaaa"
    # For text of length N of all same chars, Z[1] = N-1, Z[2] = N-2, etc.
    expected = [0,5,4,3,2,1]
    assert compute_z(text) == expected

def test_compute_z_empty_string():
    """
    Test compute_z on an empty string.
    """
    text = ""
    expected = []
    assert compute_z(text) == expected

def test_compute_z_single_char():
    """
    Test compute_z on a string with a single character.
    """
    text = "x"
    expected = [0]
    assert compute_z(text) == expected

def test_compute_z_partial_prefix_match():
    """
    Test compute_z when part of the string prefix matches itself offset.
    """
    text = "abcab"
    # prefix: "abcab"
    # Z[1] = 0, Z[2] = 0, Z[3]='a' == 'a', then 'b'=='b', so Z[3]=2, Z[4]='b'==a -> 0
    expected = [0,0,0,2,0]
    assert compute_z(text) == expected

def test_compute_z_prefix_multiple_occurrences():
    """
    Test compute_z where the prefix appears multiple times in various places.
    """
    text = "abcabcabc"
    # Z = [0,0,0,6,0,0,3,0,0]
    expected = [0,0,0,6,0,0,3,0,0]
    assert compute_z(text) == expected

def test_compute_z_nonalpha_characters():
    """
    Test compute_z on a string with non-alphabetic characters.
    """
    text = "ab#ab#ab#"
    # The Z-array is [0,0,0,0,5,0,0,0,1]
    expected = [0,0,0,0,5,0,0,0,1]
    assert compute_z(text) == expected

@pytest.mark.parametrize("text,expected", [
    ("abcabcd", [0,0,0,3,0,0,0]),
    ("xyzxyzx", [0,0,0,4,0,0,1]),
    ("aaaaa", [0,4,3,2,1]),
    ("abababab", [0,0,6,0,4,0,2,0]),
])
def test_compute_z_various_examples(text, expected):
    """
    Parametric test for compute_z with various examples for regression and coverage.
    """
    assert compute_z(text) == expected
