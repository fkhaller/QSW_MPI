#python -m numpy.f2py -h --overwrite-signature fMPI.pyf -m fMPI fMPI.f90
#python -m numpy.f2py --f90exec="mpifort" --f90flags="-fopenmp -fbounds-check" -lgomp -c fMPI.pyf ../sparse.f90 ../one_norms.f90 ../expm.f90 ../operators.f90 -m fMPI fMPI.f90
#python -m numpy.f2py -h --overwrite-signature foperators.pyf -m foperators foperators.f90
##f2py --f90exec="mpifort" --f90flags="-fopenmp -fbounds-check" -lgomp -c ../sparse.f90 ../operators.f90  -m foperators foperators.f90
#python -m numpy.f2py --f90exec="mpifort" --f90flags="-fopenmp -fbounds-check" -lgomp -c foperators.pyf ../sparse.f90 ../operators.f90 -m foperators foperators.f90
#rm ../../freeqsw/*.so
#cp *.so ../../freeqsw/
#python -m numpy.f2py -h --overwrite-signature fMPI.pyf -m fMPI fMPI.f90
python -m numpy.f2py --f90exec="mpifort" --f90flags="-g -fbounds-check" -lgomp -c fMPI.pyf ../sparse.f90 ../one_norms.f90 ../expm.f90 ../operators.f90 -m fMPI fMPI.f90
#python -m numpy.f2py -h --overwrite-signature foperators.pyf -m foperators foperators.f90
#f2py --f90exec="mpifort" --f90flags="-fopenmp -fbounds-check" -lgomp -c ../sparse.f90 ../operators.f90  -m foperators foperators.f90
python -m numpy.f2py --f90exec="mpifort" --f90flags="-g -fbounds-check" -lgomp -c foperators.pyf ../sparse.f90 ../operators.f90 -m foperators foperators.f90
rm ../../freeqsw/*.so
cp *.so ../../freeqsw/