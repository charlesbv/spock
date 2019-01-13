def find_in_read_input_order_variables(order_variables, var_to_find):
    
    index_of_var_to_find = (int)([s for s in order_variables if var_to_find in s][0].split('|')[1])

    return index_of_var_to_find
