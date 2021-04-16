## LICH-TEST

LICH-TEST - Liquid Water, Interfacial, Cubic and Hexagonal Ice Classification through Eclipsed and Staggered Conformation Template Matching - is a tool for classifying different ice types in atomistic simulations.

![Alt text](/img/lich_toc_top_1280.png?raw=true "LICHTEST")

Algorithm is published in the Journal of Physical Chemistry B (2021)
[https://doi.org/10.1021/acs.jpcb.1c01926](https://doi.org/10.1021/acs.jpcb.1c01926)

## Requirements

Versions for Matlab and Python.
This module requires the following modules:

(For Python 3)
 * Numpy
 * SciPy
 ```bash
 conda activate base  # example for using the above libraries from Anaconda
```

## Installation

Download the code

## Usage

Matlab:

```bash
cd lich_test_matlab_public_*version*/
Main_Ice_Recog_func(filename,filetype)
```

Python

Set filename (the simulation cell including water as atom type
O/OW/mW) in script.py.

```bash
cd lich_test_python_public_07012021/
python script.py
```

Note that the Python version is missing some of the functionality of
the more recent Matlab version.

## Maintainers

Current maintainers:
 * Golnaz Roudsari (usrname)
 * Farshad G. Veshki (usrname)
 * Olli H. Pakarinen (usrname)


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License (maybe)
[MIT](https://choosealicense.com/licenses/mit/)
