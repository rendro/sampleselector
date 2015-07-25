# sampleselector
Port of [Jasper Cooper's sampleselector](https://github.com/jaspercooper/sampleselector) from R to C++.

## Purpose
Given a set of geo locations find the biggest subset of locations where each location has no other location within a defined buffer zone.

The programm reads in a source file (csv) and writes out all found sets larger than a defined threshold.
