
#run_noiseclean.sh 

#Modified version of clean_data.sh by Ajay Vibhute, Mayuri Shinde, Ajay Ratheesh 

#Mithun NPS (Updated on 9 Oct 2020)

if [ $# -ne 2 ]; then
    echo "Incorrect number of arguments"
    echo "Usage: run_noiseclean.sh indir outdir"
    echo "indir: Directory with input files (bc.evt, bc_bunch.fits, bc_livetime.fits, .gti, .mkf)"
    echo "outdir: Directory where output files will be saved"
    exit 1    
fi

datadir=$1
outdir=$2

if [ ! -d "$datadir" ]; then
    echo "indir: $datadir does not exist"
    exit 1
fi    

if [ ! -d "$outdir" ]; then
    echo "outdir does not exist, creating $outdir"
    mkdir -p $outdir
fi

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

## Get input files required

bcevt=`ls $datadir/*_level2_bc.evt`
bclivetime=`ls $datadir/*_level2_bc_livetime.fits`
gtifile=`ls $datadir/*_level2.gti`
bunchfile=`ls $datadir/*_level2_bc_bunch.fits`
mkffile=`ls $datadir/*_level2.mkf`

bcevtbase=`basename $bcevt`

echo "Executing cztnoisypixclean"

#CZTNOISYPIXCLEAN

echo "Running cztnoisypix clean for $bcevt"
ncevt=$outdir/`echo $bcevtbase|sed 's/.evt/_nc.evt/'`
ncbadpixfile=$outdir/`echo $bcevtbase|sed 's/_bc.evt/_badpix_nc.fits/'`

cztnoisypixclean infile=$bcevt caldb_badpix_file=$caldbbadpix thresholdfile=$noiseconfigfile outbadpixfile=$ncbadpixfile outfile=$ncevt clobber="y" history="y"

##CZTDATASEL

echo "Running cztdatasel for $ncevt"
dsevt=`echo $ncevt|sed 's/.evt/_ds.evt/'`

cztdatasel infile=$ncevt gtifile=$gtifile gtitype='QUAD' outfile=$dsevt clobber=yes history=yes 
	
#CZTSUPERBUNCHCLEAN

sbcevt=`echo $dsevt|sed 's/.evt/_sbc.evt/'`
sbclivetime=$outdir/`echo $bcevtbase|sed 's/_bc.evt/_bc_sbc_livetime.fits/'`
sbcbti=`echo $dsevt|sed 's/.evt/_sbc_bti.fits/'`

echo "Running cztsuperbunchclean for $dsevt"

cztsuperbunchclean infile=$dsevt inbunchfile=$bunchfile thresholdfile=$sbcconfig inlivetimefile=$bclivetime outlivetimefile=$sbclivetime outbtifile=$sbcbti outfile=$sbcevt clobber="y" history="y"	

#CZTHEAVYBUNCHCLEAN

hbcevt=`echo $sbcevt|sed 's/.evt/_hbc.evt/'`

echo "Running cztheavybunchclean for $sbcevt"

cztheavybunchclean infile=$sbcevt inbunchfile=$bunchfile caldb_lld_file=$caldblld thresholdfile=$hbcconfig outfile=$hbcevt clobber="y" history="y"

#CZTFLICKPIXCLEAN

fcevt=`echo $hbcevt|sed 's/.evt/_fc.evt/'`
fcbadpixfile=`echo $ncbadpixfile | sed 's/_badpix_nc.fits/_badpix_nc_fc.fits/'`
fcbti=`echo $hbcevt|sed 's/.evt/_fc_bti.fits/'`

echo "Running cztflickpixclean for $hbcevt"

cztflickpixclean infile=$hbcevt inbadpixfile=$ncbadpixfile thresholdfile=$fcconfig outfile=$fcevt outbadpixfile=$fcbadpixfile outbtifile=$fcbti clobber="y" history="y"

#CZTEVENTSEP

evtsing=`echo $hbcevt|sed 's/.evt/_single.evt/'`
evtdoub=`echo $hbcevt|sed 's/.evt/_double.evt/'`

echo "Running czteventseperation for $fcevt"

czteventsep infile=$fcevt outfile_single=$evtsing outfile_double=$evtdoub clobber="y" history="y"


#CZTEVTCLEAN

cleanevt=$outdir/`echo $bcevtbase|sed 's/_bc.evt/_quad_clean.evt/'`
cleandblevt=$outdir/`echo $bcevtbase|sed 's/_bc.evt/_quad_clean.dblevt/'`

echo "Running cztevtclean on $evtsing"
cztevtclean infile=$evtsing outfile=$cleanevt alphaval="0" vetorange="0-0" clobber="y" history="y" IsdoubleEvent="n"

echo "Running cztevtclean on $evtdoub"
cztevtclean infile=$evtdoub outfile=$cleandblevt alphaval="0" vetorange="0-0" clobber="y" history="y" IsdoubleEvent="y"

	
#CZTBINDATA
echo "Running cztbindata"

bindataout=`echo $cleanevt|sed 's/.evt//'`
wtevt=`echo $cleanevt|sed 's/.evt/.wtevt/'`

cztbindata inevtfile=$cleanevt mkffile=$mkffile  badpixfile=$fcbadpixfile livetimefile=$sbclivetime outfile=$bindataout outevtfile=$wtevt maskWeight="no" rasrc="1" decsrc="1" badpixThreshold=0 outputtype="lc" timebinsize="1" energyrange='-' clobber="y"



