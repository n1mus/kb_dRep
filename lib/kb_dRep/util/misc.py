import re


def is_upa(o):
    return isinstance(o, str) and re.match(r'^\d+/\d+/\d+$', o)

def add_root(obj, root):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if is_upa(v):
                obj[k] = root + ';' + v
            else:
                add_root(v, root)
    elif isinstance(obj, list):
        for i, el in enumerate(obj):
            if is_upa(el):
                obj[i] = root + ';' + v
            else:
                add_root(el, root)
