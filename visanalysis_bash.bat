@ECHO OFF
ECHO Hello world
cd C:\Users\jcsimon\Documents\GitHub\visanalysis
call activate visanalysis
python process_data.py --experiment_file_directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\eyesss\JS015_x_JS251\fly_001 --rig Bruker --series_number 3
PAUSE