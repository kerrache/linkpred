#!/bin/bash

# Java
cd Java
swig -o LinkPredJava.cpp -java -c++ LinkPredJava.i
javac *.java
jar cf LinkPredJava.jar *.class
cd ..

# Python
cd Python
swig -o LinkPredPython.cpp -python -c++ LinkPredPython.i
cd ..

