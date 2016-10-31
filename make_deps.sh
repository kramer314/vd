#!/bin/bash

devnull="/dev/null"

depsdir="./deps/"
flibdir=${depsdir}"fortran-lib/"

threads=1

function usage {
    # Displays usage and arguments
    echo ""
    echo "make_deps [[-c | --clean] | [-h | --help]] | [-j | --jobs [N] ]]"
    echo ""
    echo "--clean: removes dependency source and build subdirectories."
    echo "--help: print this usage statement"
    echo ""
    echo "--jobs: number of build threads"
    echo ""
    echo "If no arguments are specified, dependencies will be updated and"\
	 "built in serial."
    echo ""
}


function clean {
    # Removes dependency source and builds
    rm -rf $depsdir
    echo ""
    echo "All subdirectories associated with dependencies removed."
    echo ""
}


function build {
    # Gets and builds up-to-date dependencies
    echo "=================================="
    echo "Updating and building dependencies"
    echo "=================================="
    echo ""

    mkdir -p $depsdir

    # Pull and build generic Fortran helper library
    echo "Updating generic Fortran helper library source"
    echo "----------------------------------------------"
    git -C $flibdir pull origin master || git clone https://github.com/kramer314/fortran-lib.git $flibdir
    echo ""

    echo "Building generic Fortran helper library"
    echo "---------------------------------------"
    pushd $flibdir > $devnull
    scons -j ${threads}
    popd > $devnull
    echo ""

    echo "Dependency build complete; check log above for any build errors."
    echo ""
}

# Parse arguments
if [ -z "$1" ]; then
    build
    exit
fi

while [ "$1" != "" ]; do
    case $1 in
	-c | --clean)
	    clean
	    exit
	    ;;
	-h | --help)
	    usage
	    exit
	    ;;
	-j )
	    shift
	    threads=$1
	    build
	    exit
	    ;;
	* )
	    usage
	    exit
    esac
    shift
done
