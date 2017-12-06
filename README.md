# Execution:
1) python main_Preprocess.py
2) python main_Preprocess.py

Directiory preparation, charcode change, and intraceullar domain detection.
Please execute it twice.
Please remind the number of directories that you will make. 
{1..37} is default setting.
If the number of self.IDs in sub171111.py exceeds it,
please modify the number of directories: {1..37}.

3) python master_PBS.py

Simulation.
If you do not want to use PBS. Use:
python main_single.py

4) python Postprocess.py

It summarizes many of simulation results.

5) matlab -nodisplay -nojvm -r main_Postprocess2

Obtain total spine length, head length, and neck length. 

python main_Plot1.py

python main_Plot2.py



