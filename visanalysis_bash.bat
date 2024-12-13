@ECHO OFF
ECHO Hello world
cd C:\Users\jcsimon\Documents\GitHub\visanalysis
call activate visanalysis
python attach_metadata.py --dataset_path 'C:\Users\jcsimon\Documents\Stanford\Data\Bruker\eyesss\JS015_x_JS251\fly_002' --experiment_name 'fly' --rig 'Bruker'
PAUSE