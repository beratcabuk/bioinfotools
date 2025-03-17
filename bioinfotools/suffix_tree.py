from bioinfotools.utils import BASES

class ArgBasedSingletonMeta(type):
    """
    Metaclass that ensures the same instance is returned for the same __init__ arguments.
    """
    _instances = {}

    def __call__(cls, *args, **kwargs):
        key = (cls, args, frozenset(kwargs.items()))
        if key not in cls._instances:
            cls._instances[key] = super().__call__(*args, **kwargs)
        return cls._instances[key]

class SuffixTreeNode:
    """
    Represents a node in the suffix tree.
    Each node contains child nodes corresponding to substrings.
    """
    def __init__(self):
        self.children = {}  # Dictionary mapping characters to child nodes
        self.indexes = []   # Stores indexes of suffixes starting from this node

class SuffixTree(metaclass=ArgBasedSingletonMeta):
    """
    Constructs and represents a suffix tree for a given text.
    """
    def __init__(self, text):
        self.root = SuffixTreeNode()
        self.text = text
        self.build_suffix_tree()

    def build_suffix_tree(self):
        """
        Constructs the suffix tree by inserting all suffixes of the given text.
        """
        for i in range(len(self.text)):
            self._insert_suffix(self.text[i:], i)

    def _insert_suffix(self, suffix, index):
        """
        Inserts a suffix into the suffix tree.

        Parameters:
            suffix (str): The suffix to insert.
            index (int): The starting index of the suffix in the original text.
        """
        node: SuffixTreeNode = self.root
        for char in suffix:
            if char not in node.children:
                node.children[char] = SuffixTreeNode()
            node = node.children[char]
            node.indexes.append(index)

    def search(self, pattern):
        """
        Searches for a pattern in the suffix tree.

        Parameters:
            pattern (str): The pattern to search for.

        Returns:
            list: List of starting indices where the pattern occurs in the text.
        """
        node = self.root
        for char in pattern:
            if char not in node.children:
                return []
            node = node.children[char]
        return node.indexes
