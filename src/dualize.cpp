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

#include <dipha/includes.h>

void print_help()
{
  std::cerr << "Usage: " << "dualize " << "input_filename output_filename" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << std::endl;
  std::cerr << "--help    --  prints this screen" << std::endl;
  std::cerr << "--benchmark --  prints timing info" << std::endl;
}

void print_help_and_exit()
{
  print_help();
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

void parse_command_line(int argc, char** argv, bool& benchmark, std::string& input_filename, std::string& output_filename)
{
  if (argc < 3)
    print_help_and_exit();

  input_filename = argv[argc - 2];
  output_filename = argv[argc - 1];

  for (int idx = 1; idx < argc - 2; idx++)
  {
    const std::string option = argv[idx];
    if (option == "--benchmark") benchmark = true;
    else print_help_and_exit();
  }
}

int main(int argc, char** argv)
{
  // mandatory MPI initilization call
  MPI_Init(&argc, &argv);

  // take time at beggining of execution
  double time_at_start = MPI_Wtime();

  std::string input_filename; // name of file that contains the weighted cell complex
  std::string output_filename; // name of file that will contain the persistence diagram
  bool benchmark = false; // print timings / info
  parse_command_line(argc, argv, benchmark, input_filename, output_filename);

  if (dipha::file_types::get_file_type(input_filename) == dipha::file_types::WEIGHTED_BOUNDARY_MATRIX)
  {
    dipha::algorithms::dualize_explicit_complex(input_filename, output_filename);
    if (benchmark)
    {
      dipha::mpi_utils::cout_if_root() << std::endl << "Overall running time in seconds: " << std::endl;
      dipha::mpi_utils::cout_if_root() << std::setprecision(1) << MPI_Wtime() - time_at_start << std::endl;
    }
  }
  else
  {
    dipha::mpi_utils::error_printer_if_root() << "Input file " << input_filename << " is not of type WEIGHTED_BOUNDARY_MATRIX" 
                                              << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  MPI_Finalize();
}
