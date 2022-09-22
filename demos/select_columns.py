import pandas as pd

test = pd.DataFrame({'a': ['a', 'b'], 1: [1, 2], 2: [3, 4]})
print(test)
print(test[[1,2]])