#!/bin/bash
# LMU Munich - AG Enard / Karolinska Institute - Sandberg Lab
# Pipeline to run UMI-seq analysis from fastq to read count tables.
# Authors: Swati Parekh, Christoph Ziegenhain, Beate Vieth & Ines Hellmann
# Contact: sparekh@age.mpg.de or christoph.ziegenhain@ki.se
vers=2.9.7
currentv=$(curl -s https://raw.githubusercontent.com/sdparekh/zUMIs/main/zUMIs.sh | grep '^vers=' | cut -f2 -d "=")
if [ "$currentv" != "$vers" ] ; then
    echo -e "------------- \n\n Good news! A newer version of zUMIs is available at https://github.com/sdparekh/zUMIs \n\n-------------";
fi

function check_opts() {
    value=$1
    name=$2
    flag=$3

    if [[ -z "${value}" ]] ; then
        failure "No ${name}!! One can not run this pipeline without ${flag} option."
    fi
}

function failure() {
	echo -e "\n There seems to be a problem. Please check the usage: \n $1 \n\n"
	usage
	exit 1
}

zumis=$0

function usage () {
    cat >&2 <<EOF

  USAGE: ${zumis} [options]
	-h  Print the usage info.

## Required parameters ##

	-y  <YAML config file> : Path to the YAML config file. Required.

## Program path ##
	-d  <zUMIs-dir>   	 : Directory containing zUMIs scripts.  Default: path to this script.

## Miniconda environment

  -c : Use zUMIs dependencies in the preinstalled conda enviroment.

zUMIs version ${vers}

EOF
}

# Define the default variables #
zumisdir=$(dirname $(readlink -f $0))


while getopts ":y:d:ch" options; do #Putting <:> between keys implies that they can not be called without an argument.
  case ${options} in
  y ) yaml=${OPTARG};;
  d ) zumisdir=${OPTARG};;
  c ) conda=true;;
  h ) usage
          exit 1;;
  \? ) echo -e "\n This key is not available! Please check the usage again: -${OPTARG}"
  	usage
  	exit 1;;
  esac
done

if [[ ${OPTIND} -eq 1 ]] ; then
    usage
    exit 1
fi

check_opts "${yaml}" "YAML" "-y"

# create temporary YAML file for corrected options
yaml_orig=${yaml}
yaml=$(dirname ${yaml})/$(basename ${yaml} .yaml).run.yaml
cp ${yaml_orig} ${yaml}

#now get some variables from YAML
num_threads=$(grep 'num_threads' ${yaml} | awk '{print $2}')
project=$(grep 'project:' ${yaml} | awk '{print $2}')
whichStage=$(grep 'which_Stage:' ${yaml} | awk '{print $2}')
outdir=$(grep 'out_dir' ${yaml} | awk '{print $2}')
genomedir=$(grep 'STAR_index:' ${yaml} | awk '{print $2}')
mem_limit=$(grep 'mem_limit:' ${yaml} | awk '{print $2}')
isstats=$(grep 'make_stats:' ${yaml} | awk '{print $2}')
fqfiles=$(grep 'name:' ${yaml} | awk '{print $2}')
velo=$(grep 'velocyto:' ${yaml} | awk '{print $2}')
staridxdir=$(grep 'STAR_index:' ${yaml} | awk '{print $2}')


if grep -q 'samtools_exec:' ${yaml} ; then
    samtoolsexc=$(grep 'samtools_exec' ${yaml} | awk '{print $2}')
else
    samtoolsexc=samtools
    echo "samtools_exec: ${samtoolsexc}" >> ${yaml}
fi

if grep -q 'pigz_exec:' ${yaml} ; then
    pigzexc=$(grep 'pigz_exec' ${yaml} | awk '{print $2}')
else
    echo "Warning: YAML file doesn't include 'pigz_exec' option; setting to 'pigz'"
    pigzexc=pigz
    echo "pigz_exec: ${pigzexc}" >> ${yaml}
fi

if grep -q 'STAR_exec:' ${yaml} ; then
    starexc=$(grep 'STAR_exec' ${yaml} | awk '{print $2}')
else
    echo "Warning: YAML file doesn't include 'STAR_exec' option; setting to 'STAR'"
    starexc=STAR
    echo "STAR_exec: ${starexc}" >> ${yaml}
fi

if grep -q 'Rscript_exec:' ${yaml} ; then
    Rexc=$(grep 'Rscript_exec' ${yaml} | awk '{print $2}')
else
    echo "Warning: YAML file doesn't include 'Rscript_exec' option; setting to 'Rscript'"
    Rexc=Rscript
    echo "Rscript_exec: ${Rexc}" >> ${yaml}
fi

#check for conda usage!
if [[ ${conda} = true ]] ; then
  echo "Using miniconda environment for zUMIs!"
  echo " note: internal executables will be used instead of those specified in the YAML file!"
  samtoolsexc=samtools
  if grep -q 'samtools_exec:' ${yaml} ; then
      sed -i '/samtools_exec:/d' ${yaml}
  fi
  echo "samtools_exec: ${samtoolsexc}" >> ${yaml}
  pigzexc=pigz
  if grep -q 'pigz_exec:' ${yaml} ; then
      sed -i '/pigz_exec:/d' ${yaml}
  fi
  echo "pigz_exec: ${pigzexc}" >> ${yaml}
  starexc=STAR
  if grep -q 'STAR_exec:' ${yaml} ; then
      sed -i '/STAR_exec:/d' ${yaml}
  fi
  echo "STAR_exec: ${starexc}" >> ${yaml}
  Rexc=Rscript
  if grep -q 'Rscript_exec:' ${yaml} ; then
      sed -i '/Rscript_exec:/d' ${yaml}
  fi
  echo "Rscript_exec: ${Rexc}" >> ${yaml}

  zumisenv=${zumisdir}/zUMIs-env
  miniconda=${zumisdir}/zUMIs-miniconda.tar.bz2
  #check if zUMIs environment has been unpacked from tar
  if [[ ! -d ${zumisenv} ]] || [[ ${zumisdir}/zUMIs-miniconda.partaa -nt ${zumisenv} ]] ; then
    [ -d ${zumisenv} ] || mkdir -p ${zumisenv}
    cat ${zumisdir}/zUMIs-miniconda.parta* > ${miniconda}
    tar -xj --overwrite -f ${miniconda} -C ${zumisenv}
  fi
  #activate zUMIs environment!
  unset PYTHONPATH
  unset PYTHONHOME
  source ${zumisenv}/bin/activate
  conda-unpack
fi

if grep -q 'zUMIs_directory:' ${yaml} ; then
    sed -i "s|zUMIs_directory:.*|zUMIs_directory: ${zumisdir}|" ${yaml}
else
    echo "zUMIs_directory: ${zumisdir}" >> ${yaml}
fi

${Rexc} ${zumisdir}/checkyaml.R ${yaml} > ${project}.zUMIs_YAMLerror.log
iserror=$(tail ${project}.zUMIs_YAMLerror.log -n1 | awk '{print $2}')

if [[ ${iserror} -eq 1 ]] ; then
    echo "YAML file has an error. Look at the zUMIs_YAMLerror.log or contact developers."
    exit 1
fi

#create main output folder if it didn't exist
if [[ ! -d ${outdir} ]] ; then
  mkdir -p ${outdir}
  if [ $? -ne 0 ] ; then
      echo "Please provide a valid output directory path."
      exit 1
  fi
fi

echo -e "\n\n You provided these parameters:
 YAML file:	${yaml_orig}
 zUMIs directory:		${zumisdir}
 STAR executable		${starexc}
 samtools executable		${samtoolsexc}
 pigz executable		${pigzexc}
 Rscript executable		${Rexc}
 RAM limit:   ${mem_limit}
 zUMIs version ${vers} \n\n" | tee "${outdir}/${project}.zUMIs_runlog.txt"
date

#check for executables
sam_exc_check=$(which ${samtoolsexc})
pigz_exc_check=$(which ${pigzexc})
r_exc_check=$(which ${Rexc})
star_exc_check=$(which ${starexc})

if [[ -z "${sam_exc_check}" ]] ||
   [[ -z "${pigz_exc_check}" ]] ||
   [[ -z "${r_exc_check}" ]] ||
   [[ -z "${star_exc_check}" ]] ; then
    echo "One or more of your executables were not found. Please check back."
    exit 1
fi

# Check if the STAR version used for mapping and the one in the provided STAR index are the same
starver=$(${starexc} --version | sed 's/STAR_//g' | sed 's/\s+//g')
staridxver=$(grep "versionGenome" ${staridxdir}/genomeParameters.txt | awk '{print $2}' | sed 's/\s+//g')

if [[ "${starver}" != "${staridxver}" ]] ; then
  echo "WARNING: The STAR version used for mapping is ${starver} and the STAR index was created using the version ${staridxver}. This may lead to an error while mapping. If you encounter any errors at the mapping stage, please make sure to create the STAR index using STAR ${starver}."
  #exit 1
fi

#create output folders
outdir=$(grep 'out_dir' ${yaml} | awk '{print $2}')
#[ -d ${outdir} ] || mkdir ${outdir}
[ -d ${outdir}/zUMIs_output/ ] || mkdir -p ${outdir}/zUMIs_output/
[ -d ${outdir}/zUMIs_output/expression ] || mkdir -p ${outdir}/zUMIs_output/expression
[ -d ${outdir}/zUMIs_output/stats ] || mkdir -p ${outdir}/zUMIs_output/stats
[ -d ${outdir}/zUMIs_output/.tmpMerge ] || mkdir -p ${outdir}/zUMIs_output/.tmpMerge


if [[ "${whichStage}" == "Filtering" ]] ; then
  echo "Filtering..."

  f=$(cut -d' ' -f1 <(echo ${fqfiles})) # the first fastq file to determine gzip status
  fullsize=$(stat -L --printf="%s" ${f})

  tmpMerge=${outdir}/zUMIs_output/.tmpMerge/

  if [[ ${f} =~ \.gz$ ]] ; then
      ${pigzexc} -dc ${f} | head -n 4000000 | ${pigzexc} > ${tmpMerge}/${project}.1mio.check.fq.gz
      smallsize=$(stat --printf="%s" ${tmpMerge}/${project}.1mio.check.fq.gz)
      rm ${tmpMerge}/${project}.1mio.check.fq.gz
      nreads=$(expr ${fullsize} \* 1000000 / ${smallsize})

      for i in ${fqfiles} ; do bash ${zumisdir}/splitfq.sh ${i} ${pigzexc} ${num_threads} ${tmpMerge} splitfqgz ${project} ${nreads} & done
      wait
      pref=$(basename ${f} .gz)
      l=$(ls ${tmpMerge}${pref}* | sed "s|${tmpMerge}${pref}||" | sed 's/.gz//')
  else
      head -n 4000000 ${f} > ${tmpMerge}/${project}.1mio.check.fq
      smallsize=$(stat --printf="%s" ${tmpMerge}/${project}.1mio.check.fq)
      rm ${tmpMerge}/${project}.1mio.check.fq
      nreads=$(expr ${fullsize} \* 1000000 / ${smallsize})

      for i in ${fqfiles} ; do bash ${zumisdir}/splitfq.sh ${i} ${pigzexc} ${num_threads} ${tmpMerge} splitfq ${project} ${nreads} & done
      wait
      pref=$(basename ${f})
      l=$(ls ${tmpMerge}${pref}* | sed "s|${tmpMerge}${pref}||")
  fi

  for x in ${l} ; do perl ${zumisdir}/fqfilter_v2.pl ${yaml} ${samtoolsexc} ${Rexc} ${pigzexc} ${zumisdir} ${x} & done
  wait
  bash ${zumisdir}/mergeBAM.sh ${zumisdir} ${tmpMerge} ${num_threads} ${project} ${outdir} ${yaml} ${samtoolsexc}
  for i in ${fqfiles} ; do
      pref=$(basename ${i} | sed 's/.fastq.gz//' | sed 's/.fq.gz//')
      rm ${tmpMerge}${pref}*gz
  done
  date

  #run barcode detection
  ${Rexc} ${zumisdir}/zUMIs-BCdetection.R ${yaml}

  #check if BC correction should be performed!
  BCbinTable=${outdir}/zUMIs_output/"${project}".BCbinning.txt
  if [[ -f "${BCbinTable}" ]] ; then
      for x in ${l} ; do
        rawbam="${tmpMerge}/${project}.${x}.raw.tagged.bam"
        fixedbam="${tmpMerge}/${project}.${x}.filtered.tagged.bam"
        mv ${fixedbam} ${rawbam}
        perl ${zumisdir}/correct_BCtag.pl ${rawbam} ${fixedbam} ${BCbinTable} ${samtoolsexc} &
      done
      wait
  fi

fi

if
[[ "${whichStage}" == "Filtering" ]] ||
[[ "${whichStage}" == "Mapping" ]] ; then
  echo "Mapping..."
    ${Rexc} ${zumisdir}/zUMIs-mapping.R ${yaml}
  date
fi

if
[[ "${whichStage}" == "Filtering" ]] ||
[[ "${whichStage}" == "Mapping" ]] ||
[[ "${whichStage}" == "Counting" ]] ; then
  echo "Counting..."
  ${Rexc} ${zumisdir}/zUMIs-dge2.R ${yaml}
  date
  ${Rexc} ${zumisdir}/misc/rds2loom.R ${yaml}
  date
  if [[ "${velo}" == "yes" ]] ; then
    ${Rexc} ${zumisdir}/runVelocyto.R ${yaml}
  fi
fi

if
[[ "${whichStage}" == "Filtering" ]] ||
[[ "${whichStage}" == "Mapping" ]] ||
[[ "${whichStage}" == "Counting" ]] ||
[[ "${whichStage}" == "Summarising" ]] ; then
  if [[ "${isstats}" == "yes" ]] ; then
    echo "Descriptive statistics..."
      ${Rexc} ${zumisdir}/zUMIs-stats2.R ${yaml}
  fi
  date
fi

#close conda enviroment if necessary
if [[ ${conda} = true ]] ; then
  source ${zumisenv}/bin/deactivate
fi
