# moca
A tool for computing mock RVs induced by exoplanets

<p align="center">
<img src="https://github.com/olgabalsa/mocha/assets/47603865/8fb9ef02-b9c8-4089-a6ee-e36e2283f0f3" width="350" />

## Warning
This documentation is under construction. We are working on provide more information on the code and its usage. If you have any question or suggestion, please contact me: obalsalobre@cab.inta-csic.es.
  
  
## Installation

Just download this folder or clone it to your computer.

## Usage

The only mandatory input you have to provide is the number of planets! Moca can simulate all the parameters you do not want to fix.

Call it from the terminal as:
```python
python run_moca.py -st M2 -np 2 -p 4.5 12.3 -mp 3.5 7.1 -star MyStar -nrv 10 -erv 3.5 -cad 12 -dir 'path/outputs/'
```

You can customizes the setting by giving additional inputs. Do ```python run_moca.py --help``` to see all the options.
