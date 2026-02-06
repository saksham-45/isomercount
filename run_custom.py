from cyclic_sequences import count_adjacency_burnside

n = 200
print(f"Calculating for n={n}...")
count = count_adjacency_burnside(n)
print(f"Result: {count}")
