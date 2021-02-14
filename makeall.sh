
# To set 
path_spock=$HOME # don't include the last '/' in the path
is_sift=0 # set this variable to 1 if you want to use SpOCK as part of SIFT. But if you want to use SpOCK independtly of SIFT, set it to 0
py_only=0 # if set to 1 and is_sift is 1 then doesnt compile SpOCK but just move the python scripts to ../ so that the demo can be run before installing SpOCK. cbv should st py_only to 1 before distributing SIFT so that users alreay ahve all these python scripts in ../ and can run the demo without running this makeall.sh
# end of to set

arg_compiler=mpicc # $1
if [ "$arg_compiler" = "gcc" ]; then
    compiler="gcc"
    path_compiler=gcc
else
    compiler="mpicc"
    path_compiler=mpicc #/opt/local/bin/mpicc-openmpi-gcc7  #$2 don't include the last '/' in the path
fi 


if [ $is_sift -eq 1 ];then
    path_spice=./cspice # when distributing: ./cspice
    path_gsl=./gsl_installation # when distributing: ./gsl_installation
else
    path_spice=${PWD}/cspice   #$3 don't include the last '/' in the path
    path_gsl=/usr/local    #$4 you can put ./ for current directory but don't put ../ for one level up. in this case, put absolute path. also, don't include the last '/' in the path
fi

# ot run SpOCK from python script with os.system(path_mpirun + " -np 4 spock in.txt") -> looks like if cygwin then the system won't recognize mpirun so need to run orterun. But if linux.unix/other than cygwin, won't regcognize orterun so need to put mpirun
if [ $OSTYPE = "cygwin" ]; then
    path_mpirun="orterun"
else
    path_mpirun="mpirun"
fi
#path_mpirun="$(which orterun)" #/opt/local/bin/mpirun-openmpi-gcc49   #$5


# if path to exectuable spock is not already in PATH then add it. Note: the user still needs to source the .bash_profile after running this makeall.sh file (not poassible to source it within this makeall.sh).
## create the varilable absolute path to executable spock
if [[ "${path_spock:0:1}" == "." ]]; then # if path_spock includes the current directory "." as first char then repalce "." with absolute path of current directory to create the absolute path to executable spock
    path_spock_after_pwd=$(cut -d "." -f 2 <<< "$path_spock") 
    path_spock_abso=$PWD$path_spock_after_pwd
elif [[ "${path_spock:0:1}" == "~" ]]; then # if path_spock includes the current directory "~" as first char then repalce "~" with absolute path of home directory to create the absolute path to executable spock
    path_spock_after_pwd=$(cut -d "." -f 0 <<< "$path_spock") # ignore warning with cut  that is likely to be printed when running this line
    path_spock_abso=$HOME$path_spock_after_pwd
else # otherwise that means that path_spock is already an absolute path
    path_spock_abso=$path_spock
fi
## now add path in .bash_profile if not already in PATH
if [[ ":$PATH:" == *"$path_spock_abso:"* ]]; then
    echo "already have path to SpOCK"
else
    echo "creating path to SpOCK"
    echo "" >> $HOME/.bash_profile
    echo "PATH=\$PATH:$path_spock_abso" >> $HOME/.bash_profile
fi

# add path to gsl lib and bin (suppsoed to prevent errors like "error while loading shared libraries"). Note: the user still needs to source the .bash_profile after running this makeall.sh file (not poassible to source it within this makeall.sh). !!! ALSO here for now only add ath to bin but mayeb need to add path to lib
## create the varilable absolute path to gsl
if [[ "${path_gsl:0:1}" == "." ]]; then # if path_gsl includes the current directory "." as first char then repalce "." with absolute path of current directory to create the absolute path to gsl
    path_gsl_after_pwd=$(cut -d "." -f 2 <<< "$path_gsl")
    path_gsl_abso=$PWD$path_gsl_after_pwd
else # otherwise that means that path_gsl is already an absolute path
    path_gsl_abso=$path_gsl
fi
## now add path in .bash_profile if not already in PATH
if [[ ":$PATH:" == *"$path_gsl_abso/bin:"* ]]; then
    echo "already have path to gsl lib and bin"
else
    echo "creating path to gsl lib and bin"
    echo "" >> $HOME/.bash_profile
    echo "PATH=\$PATH:$path_gsl_abso/bin" >> $HOME/.bash_profile
fi
# end end of add path to gsl lib and bin



## create the varilable absolute path to executable spice
if [[ "${path_spice:0:1}" == "." ]]; then # if path_spice includes the current directory "." as first char then repalce "." with absolute path of current directory to create the absolute path to executable spice
    path_spice_after_pwd=$(cut -d "." -f 2 <<< "$path_spice") 
    path_spice_abso=$PWD$path_spice_after_pwd
elif [[ "${path_spice:0:1}" == "~" ]]; then # if path_spice includes the current directory "~" as first char then repalce "~" with absolute path of home directory to create the absolute path to executable spice
    path_spice_after_pwd=$(cut -d "." -f 0 <<< "$path_spice") # ignore warning with cut  that is likely to be printed when running this line
    path_spice_abso=$HOME$path_spice_after_pwd
else # otherwise that means that path_spice is already an absolute path
    path_spice_abso=$path_spice   
fi
if [ $py_only -ne 1 ];then
    # cp SPICE files to SPICE directory if this has not already been done. Same for gravitational spherical harmonincs and specular points
    if [ ! -f "$path_spice_abso"/data/naif0012.tls ]; then
	cp ./src/naif0012.tls ./src/de432s.bsp ./src/earth_000101_210404_210111.bpc ./src/pck00010.tpc "$path_spice_abso"/data
	cp ./src/egm96_to360_not_norm.txt "$path_spice_abso"/data
	cp ./src/{antenna.bin,antenna.info,sigma0_table.bin} "$path_spice_abso"/data
    fi
    if [ $is_sift -eq 1 ];then # if compiling SpOCK to use it as part of SIFT then copy files for SIFT
	cp ./src/{cygnss_geometry_2016.txt,my_ground_stations.txt} ../../input_sift/
    fi

    # Compile SpOCK (and specular point codes)
    if [ $compiler = "gcc" ]; then
	# SpOCK
	make clean PATH_EXECUTABLE="$path_spock_abso"
	make all COMPILER=gcc PATH_COMPILER=$path_compiler PATH_SPICE="$path_spice_abso" PATH_GSL="$path_gsl_abso" PATH_EXECUTABLE="$path_spock_abso"
	# # Specular points binary files
	# $path_compiler -c -o ./src/mpi_fake/mpi.o ./src/mpi_fake/mpi.c 
	# ar r ./src/mpi_fake/libmpi.a ./src/mpi_fake/mpi.o 
	# ranlib ./src/mpi_fake/libmpi.a 
	# $path_compiler -c -o ./src/find_specular_points.o ./src/find_specular_points.c
	# $path_compiler -o "$path_spock_abso"/spec ./src/find_specular_points.o -L./src/mpi_fake -lmpi -lm  # change the path of the executable here
	# # Convert specular points binary files to txt files
	cd ./src/storm
	# make clean PATH_EXECUTABLE="$path_spock_abso"
	# make all COMPILER=gcc PATH_COMPILER=$path_compiler PATH_SPICE="$path_spice_abso" PATH_GSL="$path_gsl_abso" PATH_EXECUTABLE="$path_spock_abso"
	
    else
	# SpOCK
	make clean PATH_EXECUTABLE="$path_spock_abso"
	make all PATH_COMPILER=$path_compiler PATH_SPICE="$path_spice_abso" PATH_GSL="$path_gsl_abso" PATH_EXECUTABLE="$path_spock_abso"
	#echo sudo /usr/libexec/ApplicationFirewall/socketfilterfw -add "$path_spock_abso"/spock_dev
	#sudo /usr/libexec/ApplicationFirewall/socketfilterfw -add "$path_spock_abso"/spock_dev
	# # Specular points binary files
	#echo $path_compiler -o "$path_spock_abso"/spec ./src/find_specular_points_aspherical_dev.c -lm -w -I"$path_spice_abso" -I"$path_spice_abso"/include "$path_spice_abso"/lib/csupport.a "$path_spice_abso"/lib/cspice.a # change the path of the executable here
	echo $path_compiler -o "$path_spock_abso"/spec_gnss_rcg ./src/find_specular_points_aspherical_dev_gnss.c -lm -w -I"$path_spice_abso" -I"$path_spice_abso"/include "$path_spice_abso"/lib/csupport.a "$path_spice_abso"/lib/cspice.a # change the path of the executable here	# added on 2020-05-19
# 	sudo /usr/libexec/ApplicationFirewall/socketfilterfw -add "$path_spock_abso"/spec_gnss
# 	# Convert specular points binary files to txt files
	#cd ./src/storm
	#echo make all PATH_COMPILER=$path_compiler PATH_SPICE="$path_spice_abso" PATH_GSL="$path_gsl_abso" PATH_EXECUTABLE="$path_spock_abso"
	#make clean PATH_EXECUTABLE="$path_spock_abso"
	#make all PATH_COMPILER=$path_compiler PATH_SPICE="$path_spice_abso" PATH_GSL="$path_gsl_abso" PATH_EXECUTABLE="$path_spock_abso"
    fi
    #cd ../..

fi

exit

# Move python scripts to executable folders
if [ $is_sift -eq 1 ];then # if compiling SpOCK to run it within SIFT
    cp ./srcPython/{spock_cygnss_spec_parallel_for_sift.py,spock_main_input.py,cygnss_tle_for_sift.py,gps_tle_for_sift.py,report_coverage_ground_station_for_sift_parallel.py,nb_usable_proc.py,read_input_file.py,find_in_read_input_order_variables.py} ../ # moves these python scripts from srcPython (part of SpOCK) to the folder with the SIFT python scripts (the only diffence in spock_cygnss_spec_parallel_for_sift.py and spock_cygnss_spec_parallel.py is that spock_cygnss_spec_parallel_for_sift.py calls cygnss_tle.py snd gps_tle.py as python scripts not as executables
    i="spock_cygnss_spec_parallel_for_sift.py"
    echo "# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT">>  ../"$i"_python_path
    echo "path_mpirun = '$path_mpirun'">> ../"$i"_python_path
    echo "spice_path = '$path_spice_abso/data'">> ../"$i"_python_path
    echo "# end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT">> ../"$i"_python_path
    cat ../"$i" >> ../"$i"_python_path
    mv ../"$i"_python_path ../"$i"
    #chmod +x /"$i"

else # if compiling SpOCK indepently of SIFT
    declare -a py_script=("spock_cygnss_spec_parallel.py" "spock_cygnss_spec_parallel_for_sift_sftp.py"
                "cygnss_tle.py" "gps_tle.py" "report_coverage_ground_station_for_sift_parallel_sftp.py"
                "spock_main_input.py" "read_input_file.py"
		"norad_id_to_cygnss_name.py" "cygnss_read_spock_spec.py" "cygnss_read_spock_spec_no_prn.py" "cygnss_read_spock_spec_bin.py"
		"read_output_file.py" 
		 "report_coverage_ground_station.py"
		"spock_animation.py" "download_tle.py" "spock_run_tle.py"
		"kalman_9state.py" "noise.py" "noise_drag.py" "nb_usable_proc.py"
		)

    ## now loop through the above array
    if [ ! -d "$path_spock_abso" ]; then
	mkdir "$path_spock_abso"
    fi


    PYTHONPATH='#!/Users/cbv/Library/Enthought/Canopy_64bit/User/bin/python'
    #'#!'"$(which python)" # Python path
    for i in "${py_script[@]}"
    do
	cp -p ./srcPython/"$i" "$path_spock_abso"/"$i"
	# need to add this line at top of each python executable (path of Python)
	echo "${PYTHONPATH}" > "$path_spock_abso"/"$i"_python_path
	if [ $i = "spock_cygnss_spec_parallel.py" ]; then # for certain Python scripts, need to tell the path to the program that runs the executables (mpirun for instance)
	    echo "# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT">> "$path_spock_abso"/"$i"_python_path
	    echo "path_mpirun = '$path_mpirun'">> "$path_spock_abso"/"$i"_python_path
	    echo "spice_path = '$path_spice_abso/data'">> "$path_spock_abso"/"$i"_python_path
	    echo "# end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT">> "$path_spock_abso"/"$i"_python_path
	fi
	if [ $i = "spock_cygnss_spec_parallel_for_sift_sftp.py" ]; then # for certain Python scripts, need to tell the path to the program that runs the executables (mpirun for instance)
	    echo "# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT">> "$path_spock_abso"/"$i"_python_path
	    echo "path_mpirun = '$path_mpirun'">> "$path_spock_abso"/"$i"_python_path
	    echo "spice_path = '$path_spice_abso/data'">> "$path_spock_abso"/"$i"_python_path
	    echo "# end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT">> "$path_spock_abso"/"$i"_python_path
	fi
	if [ $i = "spock_run_tle.py" ]; then # for certain Python scripts, need to tell the path to the program that runs the executables (mpirun for instance)
	    echo "# PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT">> "$path_spock_abso"/"$i"_python_path
	    echo "path_mpirun = '$path_mpirun'">> "$path_spock_abso"/"$i"_python_path
	    echo "spice_path = '$path_spice_abso/data'">> "$path_spock_abso"/"$i"_python_path
	    echo "# end of PARAMETERS TO SET BEFORE RUNNING THIS SCRIPT">> "$path_spock_abso"/"$i"_python_path
	fi
	
	cat "$path_spock_abso"/"$i" >> "$path_spock_abso"/"$i"_python_path
	mv "$path_spock_abso"/"$i"_python_path "$path_spock_abso"/"$i"
	chmod +x "$path_spock_abso"/"$i"
    done

    #path_mpirun = 'mpirun-openmpi-gcc49'
fi
