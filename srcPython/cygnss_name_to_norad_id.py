# this script converts the FM01, FM02, ..., FM08 designation to the NORAD id of CYGNSS

def cygnss_name_to_norad_id(cygnss_name):
    cygnss_norad_array = ['41884', '41885', '41886', '41887', '41888', '41889', '41890', '41891']
    cygnss_name_array = ['FM05', 'FM04', 'FM02', 'FM01', 'FM08', 'FM06', 'FM07', 'FM03']
    cygnss_norad = cygnss_norad_array[cygnss_name_array.index(cygnss_name)]#cygnss_name_array[conversion_norad_to_cygnss_name[cygnss_norad.index(norad_id)]]
    return cygnss_norad

