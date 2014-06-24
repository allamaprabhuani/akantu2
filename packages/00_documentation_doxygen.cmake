option(AKANTU_DOCUMENTATION_DOXYGEN "Build source documentation using Doxygen." OFF)

set(AKANTU_DOCUMENTATION_DOXYGEN_DOCUMENTATION
"
This generates the Doxygen documantation of the source code.
It depends on:
\\begin{itemize}
\\item \\href{http://www.stack.nl/~dimitri/doxygen/}{Doxygen} an automated source code documentations system.
\\item Optional: \\href{http://www.graphviz.org/}{Graphviz} to generate the dependencies graph
\\end{itemize}

Under Ubuntu (14.04 LTS), the installation of the dependencies can be performed using the following command:
\\begin{command}
  > sudo apt-get install doxygen
  > sudo apt-get install graphviz
\\end{command}
")