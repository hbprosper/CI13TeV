export CIPATH="$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PATH=$CIPATH/bin:$PATH
export LD_LIBRARY_PATH=$CIPATH/CI/lib:$CIPATH/CI/src:$LD_LIBRARY_PATH
export LHAPDF_DATA_PATH=$HOME/external/share/LHAPDF/PDFsets
echo "LHAPDF_DATA_PATH=$LHAPDF_DATA_PATH"
echo "CIPATH=$CIPATH"

