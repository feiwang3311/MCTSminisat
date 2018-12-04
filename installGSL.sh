target="$(pwd)/GSL"
mkdir $target
echo "building GSL at $target"
wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz
tar -zxvf gsl-2.4.tar.gz
cd gsl-2.4
./configure --prefix=$target
make
make check
make install
