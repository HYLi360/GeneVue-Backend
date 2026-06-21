from genevue.Taxonomy.SimpleTree import SimpleTree
import pandas as pd

st = SimpleTree()
st.add_paths(pd.read_csv("res.tsv", sep="\t"))
st.build()
print(st.to_newick())
print(st.draw_as_ascii())
