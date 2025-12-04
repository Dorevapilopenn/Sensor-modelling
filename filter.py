import pandas as pd

df = pd.read_csv("cons9U_sortedU.csv", header=None);   # no column names

filtered = df[df[4] >= df[5]]; # column 6 >= column 5 (0-based so this is actually col7>=col6 if you meant literal numbering)
filtered.to_csv("cons9U_sortedN.csv", index=False, header=False);
