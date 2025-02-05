@ECHO OFF
ECHO Hello world
cd C:\Users\jcsimon\Documents\GitHub\visanalysis
call activate visanalysis
python attach_metadata.py --dataset_path C:\Users\jcsimon\Documents\Stanford\Data\Bruker\eyesss\JS015_x_JS251\fly_001 --experiment_name fly --rig Bruker
cd C:\Users\jcsimon\Documents\GitHub\visanalysis\gui
python .\DataGUI.py
PAUSE