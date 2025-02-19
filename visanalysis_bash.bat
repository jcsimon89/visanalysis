@ECHO OFF
ECHO Hello world
cd C:\Users\jcsimon\Documents\GitHub\visanalysis
call activate visanalysis
python process_data.py --experiment_file_directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\eyesss\JS139_x_JS251\fly_003 --rig Bruker --series_number 3 --run_gui True --attach_metadata True
PAUSE
python analyze_data.py --experiment_file_directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\eyesss\JS139_x_JS251\fly_003 --rig Bruker --show_figs False --save_figs True
PAUSE
select_final_rois.py --experiment_file_directory C:\Users\jcsimon\Documents\Stanford\Data\Bruker\eyesss\JS139_x_JS251\fly_003 --rig Bruker --save True
PAUSE