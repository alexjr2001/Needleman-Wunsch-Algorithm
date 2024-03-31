import matplotlib.pyplot as plt

with open("output.txt", "r") as archivo:
    lineas = archivo.readlines()

seq1 = list(lineas[0].strip())
seq2 = list(lineas[1].strip())

# Assuming seq1 and seq2 are your two strings of the same size
# and that you have already loaded them from the file as you have in your code.

# Define smaller figure size
plt.figure(figsize=(8, 8))  # You can adjust these values according to your needs

# Iterate over each sequence position
for i in range(len(seq1)):
    for j in range(len(seq2)):
        # If the letters in positions i and j are equal, put a square
        if seq1[i] == seq2[j]:
            plt.plot(i, j, marker='s', color='black', markersize=5)  # ’s' represents a square

# Adjust the view
plt.xticks(range(len(seq1)), seq1)
plt.yticks(range(len(seq2)), seq2)
plt.xlabel('Seq1')
plt.ylabel('Seq2')
plt.title('Matriz de coincidencias')

# Show matrix
plt.grid(True)
plt.show()