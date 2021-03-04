from dotmap import DotMap
import pandas as pd

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 80)


config = dict(
    debug=True,
)

app = DotMap(config) # global

def reset_globals():
    app.clear()
    app.update(config)



def ref_leaf(ref):
    return ref.split(';')[-1]

def file_safe_ref(ref):
    return ref.replace('/', '.')

TRANSFORM_NAME_SEP = '__' # separate UPA, object names,bin name
