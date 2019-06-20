"C:\Program Files\MATLAB\R2018a\bin\matlab.exe" -wait -nosplash -r "run('First_round.m');"
"C:\Program Files\MATLAB\R2018a\bin\matlab.exe" -wait -nosplash -r "run('localPatch.m');"
python.exe Fitting_1.py
"C:\Program Files\MATLAB\R2018a\bin\matlab.exe" -wait -nosplash -r "run('Second_round.m');"
python.exe Fitting_2.py
"C:\Program Files\MATLAB\R2018a\bin\matlab.exe" -wait -nosplash -r "run('ARAP.m');"
