@ECHO OFF
ECHO Hello world
cd C:\Users\jcsimon\Documents\GitHub\visanalysis
call activate visanalysis
python process_data.py --experiment_file_directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\eyesss\JS139_x_JS252\fly_001 --rig Bruker --series_number 3 --run_gui True --attach_metadata True
