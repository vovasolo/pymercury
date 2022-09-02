# pymercury
Position reconstruction in scintillation detectors

## Installation
Check if the following dependencies are satisfied: 
* eigen3 linear algebra library (on Linux, it can be normally installed from repositories) with
  `sudo apt install libeigen3-dev` (Debian-based) or 
  `yum install eigen3-devel` (RPM-based)
* CERN ROOT package with Minuit2
* Python 3 development headers and pybind11 (`pip install pybind11`)

Build shared libraries `lrmodel.so` and `mercury.so`:
```bash
make pylrm
make pymercury
```

Then copy the shared libraries into the place where they will be picked up by your Python interpreter. 
You can get the list of such places by running the following script:
```python
import sys
print(sys.path)
```

To import, use
```python
import lrmodel
import mercury
```
