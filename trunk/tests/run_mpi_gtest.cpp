/*  Copyright 2014 IST Austria

Contributed by: Jan Reininghaus

This file is part of DIPHA.

DIPHA is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DIPHA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with DIPHA.  If not, see <http://www.gnu.org/licenses/>. */

#include <cstdlib>
#include <string>
#include <iostream>

void print_help( )
{
    std::cerr << "Usage: " << "run_mpi_gtest " << "mpiexec_filename executable_filename" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << std::endl;
    std::cerr << "--help    --  prints this screen" << std::endl;
}

void print_help_and_exit( )
{
    print_help( );
    std::exit( 1 );
}

void parse_command_line( int argc, char** argv, std::string& mpiexec_filename, std::string& executable_filename )
{

    if( argc != 3 )
        print_help_and_exit( );

    executable_filename = argv[ argc - 1 ];
    mpiexec_filename = argv[ argc - 2 ];
}

int main( int argc, char **argv )
{
    std::string executable_filename; // executable
    std::string mpiexec_filename; // mpiexec
    parse_command_line( argc, argv, mpiexec_filename, executable_filename );
    std::string command = std::string( "\"" ) + mpiexec_filename + std::string( "\"" );
    command += " -n 1 ";
    command += executable_filename;
    std::cout << std::endl << command << std::endl;
    return std::system( command.c_str() );
}

//
//add_test( NAME unit_tests
//          COMMAND ${ MPIEXEC } -n ${ NUM_PROCESSES } $<TARGET_FILE:unit_tests> --gtest_shuffle
//          --gtest_repeat = ${ NUM_TEST_RUNS } --gtest_break_on_failure --gtest_random_seed = ${ RANDOM_SEED }
//WORKING_DIRECTORY ${ CMAKE_HOME_DIRECTORY } / test_data
//)


