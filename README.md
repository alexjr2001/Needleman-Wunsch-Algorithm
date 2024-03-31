# Needleman-Wunsch Algorithm for Sequence Alignment

## Overview
This project demonstrates the use of dynamic programming with the Needleman-Wunsch algorithm to align sequences in DNA or RNA. The algorithm provides both the best alignment result and all possible alignments.

## Table of Contents
- [Description](#description)
- [Usage](#usage)
- [Files](#files)
- [Contributing](#contributing)
- [License](#license)

## Description
In this piece of code, dynamic programming is employed to implement the Needleman-Wunsch algorithm for sequence alignment. The algorithm calculates the optimal alignment score and provides the alignment paths for given sequences.

In main function you can edit the sequences you want to input. After all the process and having used the "rebuiltPath()" function you can jump into "matrizDePuntos.py" and run it, thus get a points Matrix in MatplotLib.

## Usage
1. **Compile**: Compile the code using a C++ compiler.
2. **Input**: Ensure your sequence data is stored in a file named `Sequencias.txt`. Each sequence should be in a separate line.
3. **Execute**: Run the compiled executable.
4. **Output**: The program will output the optimal alignment or all possible alignments, depending on the function uncommented in the `main` function.
5. **Visualize Matrix of Points**: After running the main program and generating the `output.txt` file, execute the `matrizDePuntos.py` script to visualize a matrix of points representing the positions where the sequences match.

## Files
- `main.cpp`: Contains the C++ code implementing the Needleman-Wunsch algorithm.
- `Sequencias.txt`: Input file containing sequences for alignment.
- `output.txt`: Output file storing the alignment result.
- `matrizDePuntos.py`: Python script to visualize a matrix of points representing sequence matches.

## Contributing
Contributions are welcome! If you'd like to contribute to this project, feel free to open an issue or submit a pull request.

## License
This project is licensed under the [MIT License](LICENSE).
