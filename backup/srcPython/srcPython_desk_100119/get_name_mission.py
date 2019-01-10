# This function is used to save results in a folder that has a particular name: name_mission
def get_name_mission(spacecraft_name_output):

    if ( 'cygnss' in spacecraft_name_output.lower() ):
        name_mission = 'CYGNSS' 
    elif ( 'cadre' in spacecraft_name_output.lower() ):
        name_mission = 'CADRE' 
    elif ( 'aerie' in spacecraft_name_output.lower() ):
        name_mission = 'AERIE' 
    elif ( 'scion' in spacecraft_name_output.lower() ):
        name_mission = 'SCION' 
    elif ( 'qb50' in spacecraft_name_output.lower() ):
        name_mission = 'QB50' 
    else:
        name_mission = 'other' 
        
    return name_mission
