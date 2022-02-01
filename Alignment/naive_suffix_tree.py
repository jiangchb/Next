import os

def get_suffixes(string):
    """[summary]  
    list comprehension to get the suffix of the string

    Args:
        string ([string]): [description]

    Returns:
        [list]: [suffix]
    """
    return [string[i:] + "$" for i in range(len(string) + 1)] 

def shared_prefix(string_a, string_b):
    """[summary]
    common prefix: get the longest common prefix from the list
    """
    return os.path.commonprefix([string_a, string_b]) # common prefix: get the longest common prefix from the list


def insert_suffix(string, suffix_tree):

    #suffix_tree is the dictionary to store the information
    #string is the information

    if len(suffix_tree) == 0:
        suffix_tree[string] = []
        return suffix_tree

    found_match = False

    for key in list(suffix_tree):
        prefix = shared_prefix(string, key)
        n = len(prefix)
        if len(prefix) > 0:
            found_match = True
            key_suffix = key[n:]
            string_suffix = string[n:]
            del suffix_tree[key]
            suffix_tree[prefix] = [key_suffix, string_suffix]

    if not found_match:
        suffix_tree[string] = []
    return suffix_tree

def construct_suffix_tree(string):
    tree = {} # use disctionary to store the result
    for suffix in get_suffixes(string):
        tree = insert_suffix(suffix, tree)
    return tree

my_string = "banana"
print(construct_suffix_tree(my_string))