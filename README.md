# MNK Landspace Problems written in Python

## Description

Implementation of the MNK Problems in Python.
This is meant to make our experiment runs faster when the solver is implemented in Python, as we no longer need to do C++ calls on every evaluation.

However, it tries to follow the structure of the C++ code as much as possible. This means that it may not be the most efficient code possible in Python

## TODOs: 
- Check the generation of mpList. It is a randint() from 0 to 0xFFFFFF, but is there a proper way to do this?
- Change the implementation to focus on efficiency rather than similarity to the C++ code
