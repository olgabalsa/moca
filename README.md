# mocha
A tool for creating mock RVs induced by exoplanets

<p align="center">
<img src="https://github.com/olgabalsa/mocha/assets/47603865/6acc7d7b-453a-4024-970d-59e9d2e74f6d" width="450" />

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
[TO BE DONE, file as input].
You can also give a file with previous RVs and Mocha can simulate new measurement with same uncertainties and noise. Additionally, if you do not provide addicional parameters, Mocha can infer the planetary orbital parameters from that previous dataset. Once again, you only have to worry about providing the number of planets. The format of your file must be '.ascii', '.txt', or '.csv', with the columns 'BJD', 'RV', 'ERV', and 'INSTRUMENT' separated by commas (,). The units of the RVs and their uncertainties must be m/s.

```python
python run_KOBEsim.py -file 'my_previous_RVs.txt' -np 3
```

You can customize the setting by giving additional inputs (to be specified here!).
