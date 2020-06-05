#clean_data.sh
#Ajay Vibhute,
#Mayuri Shinde
#19 Sept 2017 
# Updated by Ajay Ratheesh on 10 Nov 2017

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


#datadir="/home/cztipoc/Mayuri/cztipoc/noise_reduction/data_for_testing/grbtest/orbit_04229"
#mkffile=$datadir/AS1G05_240T01_9000000536_04229czt_level2.mkf
#saaconfig=$datadir/saaThreshold
#eventfile=$datadir/AS1G05_240T01_9000000536_04229cztM0_level2.evt
#caldbbadpix=$datadir/AS1cztbadpix20160809v01.fits
#bunchfile=$datadir/AS1G05_240T01_9000000536_04229cztM0_level2_bunch.fits
#noiseconfigfile=$datadir/noiseReductionThreshold
#caldblld=$datadir/AS1cztlld20160517v01.fits

mkffile=$1
eventfile=$2
bunchfile=$3
mkfthresh=$4
usergti=$5
caldb_path=$CALDB/data/as1/czti/bcf
caldbbadpix=`ls $caldb_path | grep AS1cztbadpix201|tail -n1`
caldbbadpix=$caldb_path/$caldbbadpix
caldblld=`ls $caldb_path | grep AS1cztlld201|tail -n1`
caldblld=$caldb_path/$caldblld
saaconfig=$CZTNOISECLEAN/config/saaThreshold.txt
noiseconfigfile=$CZTNOISECLEAN/config/noiseReductionThreshold.txt
sbcconfig=$CZTNOISECLEAN/config/superBunchCleanThreshold.txt
hbcconfig=$CZTNOISECLEAN/config/heavyBunchThreshold.txt
fcconfig=$CZTNOISECLEAN/config/flickPixThreshold.txt

bunchcleanfile=`echo $(basename $eventfile)|sed 's/.evt/_bc.evt/'`
bunchlivetime=`echo $(basename $eventfile)|sed 's/.evt/_bc_livetime.fits/'`

#CZTBUNCHCLEAN
cztbunchclean par_infile=$eventfile par_bunchfile=$bunchfile par_outfile=$bunchcleanfile par_livetimefile=$bunchlivetime par_skipT1=0.0 par_skipT2=0.0 par_skipT3=0.0 par_bunchdeftime=20 par_bunch_length_thresh=20 par_livetime_binsize=1.0 clobber=yes history=yes


#CZTDATASELECTION
dsoutput=$(basename $bunchcleanfile)
dsoutput=`echo $dsoutput| sed 's/.evt//'`
#echo "Executing data selection, will create multiple files";
#echo "cztdataselection infile=$bunchcleanfile mkffile=$mkffile thresholdfile=$saaconfig outfile=$dsoutput clobber="y" history="y""
#cztdataselection infile=$bunchcleanfile mkffile=$mkffile thresholdfile=$saaconfig outfile=$dsoutput clobber="y" history="y"

#temp=`echo $bunchcleanfile| sed 's/.evt//'`
#lsstr="$temp""_""*[0-9].evt"
f=`echo $bunchcleanfile`
echo "Executing cztnoisypix clean"

#CZTNOISYPIXCLEAN
echo "Running cztnoisypix clean for $f"
noisecleanedfile=`echo $f|sed 's/.evt/_nc.evt/'`
noisebadpixfile=`echo $f|sed 's/.evt/_badpix_nc.fits/'`
echo "cztnoisypixclean infile=$f caldb_lld_file=$caldbbadpix thresholdfile=$noiseconfigfile outbadpixfile=$noisebadpixfile outfile=$noisecleanedfile clobber="y" history="y""
cztnoisypixclean infile=$f caldb_badpix_file=$caldbbadpix thresholdfile=$noiseconfigfile outbadpixfile=$noisebadpixfile outfile=$noisecleanedfile clobber="y" history="y"
#cztnoisypixclean $f $caldbbadpix $noiseconfigfile
#bunchlivetime=`echo $f|sed 's/.evt/_nc_livetime.fits/'`
	
#CZTSUPERBUNCHCLEAN
#orbitno=`echo $noisecleanedfile |cut -d'_' -f 4|sed 's/cztM0//'`
superbunchcleanevt=`echo $f|sed 's/.evt/_nc_sbc.evt/'`
superbunchlivetime=`echo $f|sed 's/.evt/_nc_sbc_livetime.fits/'`
superbunchbti=`echo $f|sed 's/.evt/_nc_sbc_bti.fits/'`
echo "Running cztsuperbunchclean for $noisecleanedfile"
echo "cztsuperbunchclean infile=$noisecleanedfile inbunchfile=$bunchfile thresholdfile=$sbcconfig inlivetimefile=$bunchlivetime outlivetimefile=$superbunchlivetime outbtifile=$superbunchbti outfile=$superbunchcleanevt clobber="y" history="y""
cztsuperbunchclean infile=$noisecleanedfile inbunchfile=$bunchfile thresholdfile=$sbcconfig inlivetimefile=$bunchlivetime outlivetimefile=$superbunchlivetime outbtifile=$superbunchbti outfile=$superbunchcleanevt clobber="y" history="y"	
#cztsuperbunchclean $noisecleanedfile $bunchfile $bunchlivetime $orbitno 

#CZTHEAVYBUNCHCLEAN
heavybunchcleanevt=`echo $f|sed 's/.evt/_nc_sbc_hbc.evt/'`
echo "Running cztheavybunchclean for $superbunchcleanevt"
echo "cztheavybunchclean infile=$superbunchcleanevt inbunchfile=$bunchfile caldb_lld_file=$caldblld thresholdfile=$hbcconfig outfile=$heavybunchcleanevt clobber="y" history="y""
cztheavybunchclean infile=$superbunchcleanevt inbunchfile=$bunchfile caldb_lld_file=$caldblld thresholdfile=$hbcconfig outfile=$heavybunchcleanevt clobber="y" history="y"
#cztheavybunchclean $superbunchcleanevt $bunchfile $caldblld

#CZTFLICKPIXCLEAN
flickcleanevt=`echo $f|sed 's/.evt/_nc_sbc_hbc_fc.evt/'`
flickbadpixfile=`echo $f|sed 's/.evt/_badpix_nc_fc.fits/'`
flickcleanbti=`echo $f|sed 's/.evt/_nc_sbc_hbc_fc.bti/'`
echo "Running cztflickpixclean for $heavybunchcleanevt"
cztflickpixclean infile=$heavybunchcleanevt inbadpixfile=$noisebadpixfile thresholdfile=$fcconfig outfile=$flickcleanevt outbadpixfile=$flickbadpixfile outbtifile=$flickcleanbti clobber="y" history="y"

#cztflickpixclean $heavybunchcleanevt $noisebadpixfile	 

#CZTEVENTSEP
evtsep_single=`echo $f|sed 's/.evt/_nc_sbc_hbc_fc_single.evt/'`
evtsep_double=`echo $f|sed 's/.evt/_nc_sbc_hbc_fc_double.evt/'`
dscleanevt=`echo $f|sed 's/.evt/_nc_sbc_hbc_fc_single_ds.evt/'`
dscleandoubleevt=`echo $f|sed 's/.evt/_nc_sbc_hbc_fc_double_ds.evt/'`
echo "Running czteventseperation for $flickcleanevt"
czteventsep infile=$flickcleanevt outfile_single=$evtsep_single outfile_double=$evtsep_double clobber="y" history="y"
#czteventsep $flickcleanevt
     
#CZTGTIGEN	
echo "Running cztgtigen"	
gtievt=`echo $f|sed 's/.evt/_nc_sbc_hbc_fc_single.gti/'`
cztgtigen $evtsep_single mkffile=$mkffile outfile=$gtievt $mkfthresh usergtifile=$usergti
#cztgtigen $evtsep_double mkffile=$mkffile outfile=$gtievt $mkfthresh usergtifile=$usergti
	
#CZTDATASEL
echo "Running cztdatasel"
cztdatasel infile=$evtsep_single gtifile=$gtievt gtitype="quad" outfile=$dscleanevt clobber="y" history="y"
cztdatasel infile=$evtsep_double gtifile=$gtievt gtitype="quad" outfile=$dscleandoubleevt clobber="y" history="y"


#CZTEVTCLEAN
cleanevt=`echo $f|sed 's/.evt/_quad.evt/'`
cleandblevt=`echo $f|sed 's/.evt/_quad.dblevt/'`
echo "Running cztevt clean"
cztevtclean infile=$dscleanevt outfile=$cleanevt alphaval="0" vetorange="0-0" clobber="y" history="y" 
cztevtclean infile=$dscleandoubleevt outfile=$cleandblevt IsdoubleEvent="y" alphaval="0" vetorange="0-0" clobber="y" history="y" 
	
python $CZTNOISECLEAN/scripts/modify_cztcleanevtfile_cztntick.py --indir ./ --outdir ./ --Eventfile $cleanevt
python $CZTNOISECLEAN/scripts/modify_cztcleanevtfile_cztntick.py --indir ./ --outdir ./ --Eventfile $cleandblevt
	
	
	
#CZTBINDATA
#~ echo "Running cztbindata"

#~ bindataout=`echo $f|sed 's/.evt/_/'`
#~ bindataout=$bindataout"_quad_clean"
#~ wtevt=`echo $f|sed 's/.evt/.wevt/'`
#~ echo $bindataout $wtevt 

#cztbindata inevtfile=$cleanevt mkffile=$mkffile  badpixfile=$flickbadpixfile livetimefile=$superbunchlivetime outfile=$bindataout outevtfile=$wtevt maskWeight="yes" rasrc="182.6357" decsrc="39.4058" badpixThreshold=0 outputtype="both" timebinsize="1" energyrange=$i clobber="y"
#cztbindata inevtfile=$cleanevt mkffile=$mkffile  badpixfile=$flickbadpixfile livetimefile=$superbunchlivetime outfile=$bindataout outevtfile=$wtevt maskWeight="no" rasrc="1" decsrc="1" badpixThreshold=0 outputtype="lc" timebinsize="1" energyrange='-' clobber="y"

		




