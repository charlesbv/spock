# this cript converts the norad id of CYGNSS to its FM01, FM02, ..., FM08 designation

#def norad_id_to_cygnss_id(norad_id):
cygnss_norad = ['41884', '41885', '41886', '41887', '41888', '41889', '41890', '41891']
cygnss_id_array = ['FM01', 'FM02', 'FM03', 'FM04', 'FM05', 'FM06', 'FM07', 'FM08']
conversion_norad_to_cygnss_id = [4, 3, 1, 0, 7, 5, 6, 2]
for i in range(8):
    norad_id = cygnss_norad[i]
    cygnss_id = cygnss_id_array[conversion_norad_to_cygnss_id[cygnss_norad.index(norad_id)]]
    print cygnss_id
    #return cygnss_id
