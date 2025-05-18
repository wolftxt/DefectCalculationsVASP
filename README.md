# DefectCalculationsVASP
A python program that downloads input files from the NOMAD database and uses them to run VASP calculations on crystal structures with defects.

The program validates and modifies the downloaded input settings in INCAR to fit the calculation.
The program runs multiple vasp calculations. Between runs it moves all the generated files do a different directory to prevent vasp crashes and allow the user to check the data manually afterward. After running all the runs it saves the result to a .csv file.
