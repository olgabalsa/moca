# mocha
A tool for creating mock RVs induced by exoplanets

<p align="center">
<img src="https://github.com/olgabalsa/mocca/assets/47603865/9920b787-6f1f-420d-be49-1f8b8713aa0d" width="450" />

## Warning
This documentation is under construction. We are working on provide more information on the code and its usage. If you have any question or suggestion, please contact me: obalsalobre@cab.inta-csic.es.
  
  
## Installation

Just download this folder or clone it to your computer.

## Usage

The only mandatory input you have to provide is the number of planets! Mocha can simulate all the parameters you do not want to fix.

Call it from the terminal as:
```python
python run_mocha.py -st M2 -np 2 -p 4.5 12.3 -mp 3.5 7 -star MyStar -nrv 10 -erv 3.5 -cad 12 -dir 'path/outputs/'
```

You can customize the setting by giving additional inputs. Do ```python
python run_mocha.py --help``` to see all the options.
