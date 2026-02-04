import pandas as pd

df = pd.read_csv("Ms_res9U_U.csv", header=None);   # no column names

filtered = df[df[4] <= df[3]]; # column 6 >= column 5 (0-based so this is actually col7>=col6 if you meant literal numbering)
filtered.to_csv("Ms_res9U_N.csv", index=False, header=False);
