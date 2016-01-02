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

void print_help_and_exit()
{
  std::cerr << "Usage: " << "create_phat_filtration [options] input_filename output_filename" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << std::endl;
  std::cerr << "--help    --  prints this screen" << std::endl;
  std::cerr << "--upper_dim N   --  maximal dimension to compute" << std::endl;
  std::cerr << "--dual    --  saves the dualized filtation" << std::endl;
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

void parse_command_line(int argc, char** argv, bool& dualize, int64_t& upper_dim, std::string& input_filename, 
                        std::string& output_filename)
{
  if (argc < 3)
    print_help_and_exit();

  input_filename = argv[argc - 2];
  output_filename = argv[argc - 1];

  for (int idx = 1; idx < argc - 2; idx++)
  {
    const std::string option = argv[idx];
    if (option == "--help")
    {
      print_help_and_exit();
    }
    else if (option == "--dual")
    {
      dualize = true;
    }
    else if (option == "--upper_dim")
    {
      idx++;
      if (idx >= argc - 2)
        print_help_and_exit();
      std::string parameter = std::string(argv[idx]);
      size_t pos_last_digit;
      upper_dim = std::stoll(parameter, &pos_last_digit);
      if (pos_last_digit != parameter.size())
        print_help_and_exit();
    }
    else print_help_and_exit();
  }
}

template< typename Complex >
void create_phat_filtration(const std::string& input_filename, bool dualize, int64_t upper_dim, const std::string& output_filename)
{
  Complex complex;
  complex.load_binary(input_filename, upper_dim);
  dipha::data_structures::distributed_vector< int64_t > filtration_to_cell_map;
  dipha::algorithms::get_filtration_to_cell_map(complex, dualize, filtration_to_cell_map);
  dipha::data_structures::distributed_vector< int64_t > cell_to_filtration_map;
  dipha::algorithms::get_cell_to_filtration_map(complex.get_num_cells(), filtration_to_cell_map, cell_to_filtration_map);

  const int64_t nr_columns = complex.get_num_cells();
  std::vector< std::vector< int64_t > > matrix(nr_columns);
  std::vector< int64_t > dims(nr_columns, -1);

  for (int64_t cur_dim = 0; cur_dim <= complex.get_max_dim(); cur_dim++)
  {
    dipha::data_structures::flat_column_stack unreduced_columns;
    dipha::algorithms::generate_unreduced_columns(complex, filtration_to_cell_map, cell_to_filtration_map, cur_dim, dualize, 
                                                  unreduced_columns);
    dipha::data_structures::heap_column col;
    while (!unreduced_columns.empty())
    {
      int64_t index;
      unreduced_columns.pop(index, col);
      std::sort(col.begin(), col.end());
      matrix[index] = col;
      dims[index] = dualize ? complex.get_max_dim() - cur_dim : cur_dim;
    }
  }

  std::ofstream output_stream(output_filename.c_str(), std::ios_base::binary | std::ios_base::out);
  output_stream.write((char*)&nr_columns, sizeof(int64_t));
  for (int64_t cur_col = 0; cur_col < nr_columns; cur_col++)
  {
    int64_t cur_dim = dims[cur_col];
    output_stream.write((char*)&cur_dim, sizeof(int64_t));
    int64_t cur_nr_rows = matrix[cur_col].size();
    output_stream.write((char*)&cur_nr_rows, sizeof(int64_t));
    for (int64_t cur_row_idx = 0; cur_row_idx < cur_nr_rows; cur_row_idx++)
    {
      int64_t cur_row = matrix[cur_col][cur_row_idx];
      output_stream.write((char*)&cur_row, sizeof(int64_t));
    }
  }

  output_stream.close();
}

int main(int argc, char** argv)
{
  // mandatory MPI initilization call
  MPI_Init(&argc, &argv);

  if (dipha::mpi_utils::get_num_processes() != 1)
  {
    dipha::mpi_utils::error_printer_if_root() << "Error: supports only one process." << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  std::string input_filename; // name of file that contains the weighted cell complex
  std::string output_filename; // name of file that will contain the PHAT filtration
  bool dualize = false; // primal / dual computation toggle
  int64_t upper_dim = std::numeric_limits< int64_t >::max();
  parse_command_line(argc, argv, dualize, upper_dim, input_filename, output_filename);

  switch (dipha::file_types::get_file_type(input_filename))
  {
  case dipha::file_types::IMAGE_DATA:
    if (upper_dim != std::numeric_limits< int64_t >::max())
    {
      dipha::mpi_utils::error_printer_if_root() << "upper_dim not supported for this input type IMAGE_DATA" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    create_phat_filtration< dipha::inputs::weighted_cubical_complex >(input_filename, dualize, upper_dim, output_filename);
    break;
  case dipha::file_types::WEIGHTED_BOUNDARY_MATRIX:
    if (upper_dim != std::numeric_limits< int64_t >::max())
    {
      dipha::mpi_utils::error_printer_if_root() << "upper_dim not supported for this input type WEIGHTED_BOUNDARY_MATRIX" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    create_phat_filtration< dipha::inputs::weighted_explicit_complex >(input_filename, dualize, upper_dim, output_filename);
    break;
  case dipha::file_types::DISTANCE_MATRIX:
    create_phat_filtration< dipha::inputs::full_rips_complex >(input_filename, dualize, upper_dim, output_filename);
    break;
  case dipha::file_types::SPARSE_DISTANCE_MATRIX:
    create_phat_filtration< dipha::inputs::sparse_rips_complex >(input_filename, dualize, upper_dim, output_filename);
    break;
  default:
    dipha::mpi_utils::error_printer_if_root() << "Unknown complex type in DIPHA file" << input_filename << std::endl;
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }

  MPI_Finalize();
}
