#!/bin/sh

if [ -z "${DistFeatureHome}" ]; then
        echo "ERROR: Please set environmental variable DistFeatureHome to the installation folder of BuildFeatures"
        exit 1
fi

out_root=`pwd`

function Usage
{
        echo $0 "[-o out_root ]  input_A3M "
        echo "  This script generates input features for distance/orientation prediction from an MSA file in a3m format"
        echo "  input_A3M: the MSA file in a3m format, in which the first sequence shall be query and not contain any gaps."
        echo "  out_root: the output folder for result saving, default current work directory"
}

#-> parse arguments
while getopts ":o:" opt; do
        case $opt in
        o)
                out_root=$OPTARG
                ;;
        #-> help
        \?)
                echo "Invalid option: -$OPTARG" >&2
                exit 1
                ;;
        :)
                echo "Option -$OPTARG requires an argument." >&2
                exit 1
                ;;
        esac
done
shift $((OPTIND -1))
                                                           
if [ $# -ne 1 ]; then
        Usage
        exit 1
fi

input_file=$1
if [ ! -f "$input_file" ]; then
        echo "ERROR: invalid input MSA file $input_file " >&2
        exit 1
fi

## get target name
fulnam=`basename $input_file`
relnam=${fulnam%.*}

if [ ! -d $out_root ]; then
	mkdir -p $out_root
fi

# ---- generate A2M file from A3M file ---------#
a2m_file=$relnam.a2m
$DistFeatureHome/util/A3M_To_PSI $input_file $a2m_file.tmp
if [ $? -ne 0 ]; then
        echo "ERROR: failed to run $DistFeatureHome/util/A3M_To_PSI $input_file $a2m_file.tmp"
        exit 1
fi
grep -v "ss_pred\|ss_conf" $a2m_file.tmp | awk '{print substr($0,34,length($0)-32) }' > $out_root/$a2m_file
rm -f $a2m_file.tmp
