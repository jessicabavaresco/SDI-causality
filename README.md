## Code to accompany: *[Semi-device-independent certification of indefinite causal order](https://arxiv.org/abs/1903.10526)*
#### Jessica Bavaresco, Mateus Araújo, Časlav Brukner, and Marco Túlio Quintino

This is a repository for the code used to calculate the numerical results summarized in Table II of the article "*Semi-device-independent certification of indefinite causal order*, Jessica Bavaresco, Mateus Araújo, Časlav Brukner, Marco Túlio Quintino, [*Quantum* **3**, 176 (2019)](https://quantum-journal.org/papers/q-2019-08-19-176/), [arXiv:1903.10526 [quant-ph]](https://arxiv.org/abs/1903.10526)".

All code is written in MATLAB and requires:
- [Yalmip](https://yalmip.github.io) - a free MATLAB toolbox for rapid prototyping of optimization problems
- [MOSEK](https://www.mosek.com) - a software package for solving mathematical optimization problems (under the free personal academic license)
- [QETLAB](http://www.qetlab.com/) - a free MATLAB toolbox for quantum entanglement theory

This repository consists in the following four SDPs plus a function that defines the quantum switch process and sets of instruments used in our work, and calculates the non-zero values of Table II.

- [randrobust_tripartiteW_TTT](https://github.com/jessicabavaresco/SDI-causality/blob/master/randrobust_tripartiteW_TTT.m):
Calculates the \eta^* of a tripartite process matrix in the TTT scenario.

- [randrobust_tripartiteW_TTU](https://github.com/jessicabavaresco/SDI-causality/blob/master/randrobust_tripartiteW_TTU.m):
Calculates a lower bound for \eta^* of a tripartite process matrix in the TTU scenario, using a given set of POVMs for Charlie.

- [randrobust_tripartiteW_TUU](https://github.com/jessicabavaresco/SDI-causality/blob/master/randrobust_tripartiteW_TUU.m):
Calculates a lower bound for \eta^* of a tripartite process matrix in the TUU scenario, using a given set of instruments for Bob and a set of POVMs for Charlie.

- [randrobust_tripartiteW_UTT](https://github.com/jessicabavaresco/SDI-causality/blob/master/randrobust_tripartiteW_UTT.m):
Calculates a lower bound for \eta^* of a tripartite process matrix in the UTT scenario, using a given set of instruments for Alice.

- [switch_randrobustness](https://github.com/jessicabavaresco/SDI-causality/blob/master/switch_randrobustness.m):
Defines the quantum switch process and the sets of instruments used in our work, provides them as input to the four previous SDPs and outputs the values of/bounds for \eta^*.


