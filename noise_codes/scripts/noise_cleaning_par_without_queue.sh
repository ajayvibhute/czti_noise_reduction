#noise_cleaning_without_queue.sh
#Ajay Vibhute,
#Mayuri Shinde
#19 Sept 2017
#Script takes level2 file, mkf file, bunch file, saaconfig, noise config and generates 
#level2 cleaned file
#Steps :
#1.cztdataselection- for removing saa
#2.cztnoisypixclean- for removing noisy pixels
#3.cztsuperbunchclean - for removing super bunches 
#4.cztheavybunchclean - for removing heavy bunches 
#5.cztflickpixclean - for removing flickering pixels
#6.cztdatasel
#7.cztevtclean
#8.cztbindata

if [ -z "$CZTNOISECLEAN" ] 
then
	echo "CZTNOISECLEAN variable not defined"
	exit -1;
fi

if [ -z "$CZTIMONITOR" ] 
then
	echo "CZTIMONITOR variable not defined"
	exit -1;
fi

if [ -z "$CZTIMONITOR_LOG" ] 
then
	echo "'CZTIMONITOR_LOG' variable not defined"
	exit -1;
fi
logmap=$CZTIMONITOR/queues/log.map
noise_cleaning_logfile=noise_cleaning.log
saaconfig=$CZTNOISECLEAN/config/saaThreshold.txt
caldb_path=$CALDB/data/as1/czti/bcf/
caldbbadpix=`ls $caldb_path | grep AS1cztbadpix201|tail -n1`
noiseconfigfile=$CZTNOISECLEAN/config/noiseReductionThreshold.txt
caldblld=`ls $caldb_path | grep AS1cztlld201|tail -n1`
sbcconfig=$CZTNOISECLEAN/config/superBunchCleanThreshold.txt
hbcconfig=$CZTNOISECLEAN/config/heavyBunchThreshold.txt
fcconfig=$CZTNOISECLEAN/config/flickPixThreshold.txt

l2dir=$1
#checking if there is new data arrived or not
if [ ! -z $l2dir ]; then
	error_flag=0
	substring=`echo $l2dir |grep "orbit"` #for merge processing
	if [ -z $substring ]
	then
		#logfile=$(basename $l2dir)".log"
		#logfile=`echo $logfile |sed 's/_level2.log//'`
		#logfile=$CZTIMONITOR_LOG/$logfile
		obs=$(basename $l2dir)
		export GLOG_log_dir="$l2dir/czti/"
		export logfile_env=${obs:0:30}
		echo "Process logfile: $noise_cleaning_logfile"
		echo "module logfile: $logfile_env"
		echo "module log dir: $GLOG_log_dir"
		echo -e "\npnoise_cleaningsh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: Assumes that data is for full observation">>$noise_cleaning_logfile
		prefix="czti/"
		orbno=""

	#for orbit wise processing
    	else
		echo -e "\n">>$noise_cleaning_logfile
		search=`echo $l2dir|rev|cut -d/ -f 4 |rev|sed 's/_level2//' |sed 's/_//'`
		echo "Search" $search
		orbno=$(basename $l2dir)
		echo "grep $search $logmap|grep $orbno"
		logfile=`grep $search $logmap|grep $orbno|head -1|cut -d'.' -f1`
		export GLOG_log_dir="$l2dir/"
		export logfile_env=`echo $logfile|sed 's/LEVL1/LEVL2/'`
		if [ -z $logfile_env ]
		then
			echo "Log file is empty... something went wrong... exiting..."
		    echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: Unable to create log file for :: $l2dir" 
			error_flag=1
		fi
		echo "Process logfile: $noise_cleaning_logfile"
		echo "module logfile: $logfile_env"
		echo "module log dir: $GLOG_log_dir"
		echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: Assumes that data is for a orbit">>$_logfile
		prefix=""

		filename=`echo $l2dir|rev|cut -d/ -f 4 |rev |sed 's/_level2//'` 
	fi

	
	mkffile=`echo $l2dir/$prefix*.mkf`
	#mkffile=`find  $l2dir/$prefix/  -name AS1*czt*level2.mkf`
	echo $mkffile

	if [[ -z $mkffile ]];then
		echo Mkf not found, exiting...
		echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: Mkf file not found... exiting.... ">>$noise_cleaning_logfile
		echo "find  $l2dir/$prefix  -name "AS1*czt*level2.mkf" "

		echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: Mkf file not found ::$l2dir"
	fi

	m0l2evt=`find  $l2dir/$prefix"modeM0"/  -name AS1*cztM0*level2.evt`
	if [[ -z $m0l2evt ]];then
		echo  "Level 2 evt not found, exiting..."
		echo "find  $l2dir/$prefix"modeM0/"  -name AS1*cztM0*level2.evt"
		echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: Level2 evt file not found... exiting....">>$noise_cleaning_logfile

		echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: Level2 evt file not found :: $l2dir"
		#exit -1;
	fi
	
	bunchfile=`echo $m0l2evt | sed 's/.evt/_bunch.fits/'`
	if [[ -z $bunchfile ]];then
		echo  "Level 2 bunch file not found, exiting..."
		echo "find  $l2dir/$prefix"modeM0/"  -name AS1*cztM0*bunch.fits"
		echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: Level2 bunch file not found... exiting....">>$noise_cleaning_logfile

		echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: Level2 bunch file not found :: $l2dir"
		#exit -1;
	fi
	
	#temporary
	#czt bunch clean
	bunch_clean_evt=`echo $(basename $m0l2evt)| sed 's/.evt/_bc.evt/'`
	bunch_clean_livetime=`echo $(basename $m0l2evt)| sed 's/.evt/_bc_livetime.fits/'`
	if [ $error_flag -eq 0 ];then	
		echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: Executing cztbunchclean">>$noise_cleaning_logfile
		echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztbunchclean par_infile=$dsoutput par_outfile=$bunch_clean_evt par_bunchdeftime=30 par_bunch_length_thresh=20 par_skipT1=0 par_skipT2=0 par_skipT3=0 clobber=y history=y">>$noise_cleaning_logfile
		cztbunchclean par_infile=$m0l2evt par_bunchfile=$bunchfile par_outfile=$bunch_clean_evt par_livetimefile=$bunch_clean_livetime par_bunchdeftime=30 par_bunch_length_thresh=20 par_livetime_binsize=1 par_skipT1=0 par_skipT2=0 par_skipT3=0 clobber=y history=y>>$noise_cleaning_logfile

		p_status=$?
		if [ $p_status -ne 0 ];then
			error_flag=1
			echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztbunchclean Execution: `perror $p_status`  "
			echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztbunchclean Execution: `perror $p_status` ">>$noise_cleaning_logfile
		else
			echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztbunchclean execution completed"
			echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztbunchclean execution completed">>$noise_cleaning_logfile

		fi	
	fi


	#dsoutput=$(basename $m0l2evt)
	dsoutput=$(basename $bunch_clean_evt)
	dsoutput=`echo $dsoutput| sed 's/.evt//'`
	if [ $error_flag -eq 0 ];then
			echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztdataselection $mkffile $saaconfig $m0l2evt $dsoutput">>$noise_cleaning_logfile
			echo "Executing data selection, will create multiple files";
			#char *infile, char *mkffile, char *thresholdfile,char *outfile,int *clobber,int *history)
			cztdataselection infile=$bunch_clean_evt mkffile=$mkffile thresholdfile=$saaconfig outfile=$dsoutput clobber="y" history="y" 
			p_status=$?
			if [ $p_status -ne 0 ];then
				error_flag=1
				echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztdataselection Execution: `perror $p_status`  "
				echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztdataselection Execution: `perror $p_status` ">>$noise_cleaning_logfile
			else
				echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztdataselection execution completed"
				echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztdataselection execution completed">>$noise_cleaning_logfile
			fi	
	fi

	temp=`echo $m0l2evt| sed 's/.evt//'`
	lsstr="$dsoutput""_""*[0-9].evt"
	#lsstr="$dsoutput_[0-9]*.evt"
	if [ $error_flag -eq 0 ];then
		echo "Executing cztnoisypix clean"
		for f in `ls $lsstr`;
		do
			error_flag=0
			if [ $error_flag -eq 0 ];then
				echo "Running cztnoisypix clean for $f"
				noisecleanedfile=`echo $f|sed 's/.evt/_nc.evt/'`
				noisebadpixfile=`echo $f|sed 's/.evt/_badpix_nc.fits/'`
				noiselivetime=`echo $f|sed 's/.evt/_nc_livetime.fits/'`
				orbitno=`echo $noisecleanedfile |cut -d'_' -f 4|sed 's/cztM0//'`

				echo "cztnoisypixclean infile=$f caldb_badpix_file=$caldb_path$caldbbadpix thresholdfile=$noiseconfigfile outbadpixfile=$noisebadpixfile outfile=$noisecleanedfile clobber="y" history="y""

				echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztnoisypixclean infile=$f caldb_badpix_file=$caldb_path$caldbbadpix thresholdfile=$noiseconfigfile outbadpixfile=$noisebadpixfile outfile=$noisecleanedfile clobber="y" history="y"">>$noise_cleaning_logfile

			cztnoisypixclean infile=$f caldb_badpix_file=$caldb_path$caldbbadpix thresholdfile=$noiseconfigfile outbadpixfile=$noisebadpixfile outfile=$noisecleanedfile clobber="y" history="y" 
			
				#p_status=$?
				#if [ $p_status -ne 0 ];then
				#	error_flag=1
				#	echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztnoisypixclean Execution: `perror $p_status`  "
				#	echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztnoisypixclean Execution: `perror $p_status` ">>$noise_cleaning_logfile
				#else
				#	echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztnoisypixclean execution completed"
				#	echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztnoisypixclean execution completed">>$noise_cleaning_logfile
				#fi	
			fi



			if [ $error_flag -eq 0 ];then
				echo "Running cztsuperbunchclean for $noisecleanedfile"
				
				superbunchcleanevt=`echo $f|sed 's/.evt/_nc_sbc.evt/'`
				superbunchlivetime=`echo $f|sed 's/.evt/_nc_sbc_livetime.fits/'`
				superbunchbti=`echo $f|sed 's/.evt/_nc_sbc_bti.fits/'`

				echo "cztsuperbunchclean infile=$noisecleanedfile inbunchfile=$bunchfile thresholdfile=$sbcconfig inlivetimefile=$bunch_clean_livetime outlivetimefile=$superbunchlivetime outbtifile=$superbunchbti outfile=$superbunchcleanevt clobber="y" history="y""
				echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztsuperbunchclean infile=$noisecleanedfile inbunchfile=$bunchfile thresholdfile=$sbcconfig inlivetimefile=$bunch_clean_livetime outlivetimefile=$superbunchlivetime outbtifile=$superbunchbti outfile=$superbunchcleanevt clobber="y" history="y"">>$noise_cleaning_logfile
				
			cztsuperbunchclean infile=$noisecleanedfile inbunchfile=$bunchfile thresholdfile=$sbcconfig inlivetimefile=$bunch_clean_livetime outlivetimefile=$superbunchlivetime outbtifile=$superbunchbti outfile=$superbunchcleanevt clobber="y" history="y"		
				p_status=$?
				if [ $p_status -ne 0 ];then
					error_flag=1
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztsuperbunchclean Execution: `perror $p_status`  "
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztsuperbunchclean Execution: `perror $p_status` ">>$noise_cleaning_logfile
				else
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztsuperbunchclean execution completed"
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztsuperbunchclean execution completed">>$noise_cleaning_logfile
				fi	
			fi


			if [ $error_flag -eq 0 ];then
				echo "Running cztheavybunchclean for $superbunchcleanevt"
				heavybunchcleanevt=`echo $f|sed 's/.evt/_nc_sbc_hbc.evt/'`

				echo "cztheavybunchclean infile=$superbunchcleanevt inbunchfile=$bunchfile caldb_lld_file=$caldb_path$caldblld thresholdfile=$hbcconfig outfile=$heavybunchcleanevt clobber="y" history="y""
				echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztheavybunchclean infile=$superbunchcleanevt inbunchfile=$bunchfile caldb_lld_file=$caldb_path$caldblld thresholdfile=$hbcconfig outfile=$heavybunchcleanevt clobber="y" history="y"">>$noise_cleaning_logfile
				
			cztheavybunchclean infile=$superbunchcleanevt inbunchfile=$bunchfile caldb_lld_file=$caldb_path$caldblld thresholdfile=$hbcconfig outfile=$heavybunchcleanevt clobber="y" history="y"

				p_status=$?
				if [ $p_status -ne 0 ];then
					error_flag=1
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztheavybunchclean Execution: `perror $p_status`  "
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztheavybunchclean Execution: `perror $p_status` ">>$noise_cleaning_logfile
				else
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztheavybunchclean execution completed"
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztheavybunchclean execution completed">>$noise_cleaning_logfile
				fi	
			fi

			if [ $error_flag -eq 0 ];then
				echo "Running cztflickpixclean for $heavybunchcleanevt"
				flickcleanevt=`echo $f|sed 's/.evt/_nc_sbc_hbc_fc.evt/'`
				flickbadpixfile=`echo $f|sed 's/.evt/_badpix_nc_fc.fits/'`
				flickcleanbti=`echo $f|sed 's/.evt/_nc_sbc_hbc_fc.bti/'`

				echo "cztflickpixclean infile=$heavybunchcleanevt inbadpixfile=$noisebadpixfile thresholdfile=$fcconfig outfile=$flickcleanevt outbadpixfile=$flickbadpixfile outbtifile=$flickcleanbti clobber="y" history="y""
				echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S")::cztflickpixclean infile=$heavybunchcleanevt inbadpixfile=$noisebadpixfile thresholdfile=$fcconfig outfile=$flickcleanevt outbadpixfile=$flickbadpixfile outbtifile=$flickcleanbti clobber="y" history="y""
	
		cztflickpixclean infile=$heavybunchcleanevt inbadpixfile=$noisebadpixfile thresholdfile=$fcconfig outfile=$flickcleanevt outbadpixfile=$flickbadpixfile outbtifile=$flickcleanbti clobber="y" history="y"

				p_status=$?
				if [ $p_status -ne 0 ];then
					error_flag=1
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztheavybunchclean Execution: `perror $p_status`  "
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztheavybunchclean Execution: `perror $p_status` ">>$noise_cleaning_logfile
				else
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztheavybunchclean execution completed"
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztheavybunchclean execution completed">>$noise_cleaning_logfile
				fi	
			fi

			if [ $error_flag -eq 0 ];then
				echo "Running czteventseperation for $flickcleanevt"
				evtsep_single=`echo $f|sed 's/.evt/_single.evt/'`
				evtsep_double=`echo $f|sed 's/.evt/_double.evt/'`
				dscleanevt=`echo $f|sed 's/.evt/_single_ds.evt/'`
				echo "czteventsep infile=$flickcleanevt outfile_single=$evtsep_single outfile_double=$evtsep_double clobber="y" history="y""
				echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: czteventsep infile=$flickcleanevt outfile_single=$evtsep_single outfile_double=$evtsep_double clobber="y" history="y"">>$noise_cleaning_logfile

		czteventsep infile=$flickcleanevt outfile_single=$evtsep_single outfile_double=$evtsep_double clobber="y" history="y"					
				p_status=$?
				if [ $p_status -ne 0 ];then
					error_flag=1
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in czteventseperation Execution: `perror $p_status`  "
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in czteventseperation Execution: `perror $p_status` ">>$noise_cleaning_logfile
				else
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: czteventseperation execution completed"
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: czteventseperation execution completed">>$noise_cleaning_logfile
				fi	
			fi
	
			if [ $error_flag -eq 0 ];then
				echo "Running cztdatasel"
				
				echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: Executing cztdatasel with GTI TYPE QUAD ">>$noise_cleaning_logfile
				echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S")::cztdatasel infile=$evtsep_single gtifile=$evtsep_single gtitype="quad" outfile=$dscleanevt clobber="y" history="y" ">>$noise_cleaning_logfile
			cztdatasel infile=$evtsep_single gtifile=$evtsep_single gtitype="quad" outfile=$dscleanevt clobber="y" history="y"
				p_status=$?
				if [ $p_status -ne 0 ];then
					error_flag=1
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztdatasel Execution: `perror $p_status`  "
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztdatasel Execution: `perror $p_status` ">>$noise_cleaning_logfile
				else
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztdatasel execution completed"
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztdatasel execution completed">>$noise_cleaning_logfile
				fi	
			fi

			cleanevt=`echo $f|sed 's/.evt/_clean.evt/'`
			if [ $error_flag -eq 0 ];then
				echo "Running cztevt clean"
				echo "cztevtclean infile=$dscleanevt outfile=$cleanevt alphaval="0" vetorange="0-0" clobber="y" history="y" "

				echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: Executing cztevtclean on quad and keeping only alpha=0 and veto=0 ">>$noise_cleaning_logfile
				echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztevtclean infile=$quad_pc_evt outfile=$quad_clean_evt alphaval="0" vetorange="0" clobber=y isdoubleEvent=n history=y  ">>$noise_cleaning_logfile
			cztevtclean infile=$dscleanevt outfile=$cleanevt alphaval="0" vetorange="0-0" clobber="y" history="y" 	
				
				bindataout=`echo $f|sed 's/.evt/_/'`
				wtevt=`echo $f|sed 's/.evt/.wevt/'`
		
				p_status=$?
				if [ $p_status -ne 0 ];then
					error_flag=1
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztevtclean  Execution: `perror $p_status`  "
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztevtclean Execution: `perror $p_status` ">>$noise_cleaning_logfile
				else
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztevtclean execution completed"
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztevtclean execution completed">>$noise_cleaning_logfile
				fi	
			fi

			
		#	if [ $error_flag -eq 0 ];then
		#		rapnt=`fkeyprint "$m0l2evt" RA_PNT |grep RA_PNT |tail -1 |cut -d "=" -f 2|cut -d'/'  -f 1|sed "s/^[ \t]*//"`
		#		p_status=$?
		#		if [ $p_status -ne 0 ];then
		#			error_flag=1
		#			echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in fkeyprint, Please execute heainit then try: `perror $p_status`  "
		#			echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in fkeyprint, Please execute heainit then try: `perror $p_status` ">>$noise_cleaning_logfile
		#		else
		#			echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: fkeyprint execution completed"
		#			echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: fkeyprint execution completed">>$noise_cleaning_logfile
		#		fi	
		#		decpnt=`fkeyprint "$m0l2evt" DEC_PNT|grep DEC_PNT |tail -1 |cut -d "=" -f 2|cut -d'/'  -f 1|sed "s/^[ \t]*//"`
		#	fi
		#	if [ -z "$rapnt" ] 
		#	then
		#		echo "Unable to read RA PNT value"
		#		exit -1;
		#	fi

		#	if [ -z "$decpnt" ] 
		#	then
		#		echo "Unable to read DEC PNT value"
		#		exit -1;
		#	fi
		
			rapnt=1
			decpnt=1
			if [ $error_flag -eq 0 ];then
				echo "Running cztbindata"
				echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztbindata inevtfile=$cleanevt mkffile=$mkffile  badpixfile=$flickbadpixfile livetimefile=$superbunchlivetime outfile=$bindataout outevtfile=$wtevt maskWeight="no" rasrc=$rapnt decsrc=$decpnt badpixThreshold=0 outputtype="lc" timebinsize="10" energyrange="-" clobber="y"">>$noise_cleaning_logfile

				cztbindata inevtfile=$cleanevt mkffile=$mkffile  badpixfile=$flickbadpixfile livetimefile=$superbunchlivetime outfile=$bindataout outevtfile=$wtevt maskWeight="no" rasrc=$rapnt decsrc=$decpnt badpixThreshold=0 outputtype="lc" timebinsize="10" energyrange="-" clobber="y"
				p_status=$?
				if [ $p_status -ne 0 ];then
					error_flag=1
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztbindata Execution: `perror $p_status`  "
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: ERROR in cztbindata Execution: `perror $p_status` ">>$noise_cleaning_logfile
				else
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztbindata execution completed"
					echo "noise_cleaning.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: cztbindata execution completed">>$noise_cleaning_logfile
				fi	
			fi
			
		done
	fi

	echo "Script Execution completed"
	echo "noise_cleaning_without_queue.sh: $USER $(date "+%Y-%m-%d %H:%M:%S"):: Script Execution completed">>$noise_cleaning_logfile
	read	
else
	echo "Please give the level2 path of orbit or observation ID"
fi
