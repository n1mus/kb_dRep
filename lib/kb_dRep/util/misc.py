import re


def is_ref(o):
    return isinstance(o, str) and re.match(r'^([1-9]\d*/[1-9]\d*/[1-9]\d*)+$', o)

def add_root(obj, root):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if is_ref(v):
                obj[k] = root + ';' + v
            else:
                add_root(v, root)
    elif isinstance(obj, list):
        for i, el in enumerate(obj):
            if is_ref(el):
                obj[i] = root + ';' + v
            else:
                add_root(el, root)
