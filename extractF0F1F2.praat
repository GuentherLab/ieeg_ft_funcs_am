form Convert to pitch and formant files
    sentence fileDir File path
    sentence fileName File name
    positive formantCeiling hz
	positive pitchLower hz
	positive pitchUpper hz
endform

# fileDir$ = "Y:\DBS\groupanalyses\task-lombard\20210922-beh-lombard-effect-PLB\data\"
# fileName$ = "sub-DM1002_ses-intraop_task-lombard_run-03_recording-directionalmicaec"

Read from file: fileDir$ + fileName$ + ".wav"
selectObject: "Sound " + fileName$
To Formant (burg): 0.01, 5, formantCeiling, 0.025, 50
Save as short text file: fileDir$ + fileName$ + ".Formant"

selectObject: "Sound " + fileName$
To Pitch: 0.0, pitchLower, pitchUpper 
Save as short text file: fileDir$ + fileName$ + ".Pitch"