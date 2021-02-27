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

def reset(app: DotMap):
    app.clear()
    app.update(config)
