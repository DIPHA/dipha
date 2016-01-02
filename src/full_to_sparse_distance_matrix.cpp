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
  std::cerr << "Usage: " << "create_phat_filtration [options] upper_value input_filename output_filename" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Options:" << std::endl;
  std::cerr << std::endl;
  std::cerr << "--help    --  prints this screen" << std::endl;
  std::exit(-1);
}

void parse_command_line(int argc, char** argv, double& upper_value, std::string& input_filename, std::string& output_filename)
{
  if (argc < 3)
    print_help_and_exit();

  input_filename = argv[2];
  output_filename = argv[3];

  std::string parameter = std::string(argv[1]);
  size_t pos_last_digit;
  upper_value = std::stod(parameter, &pos_last_digit);
  if (pos_last_digit != parameter.size())
    print_help_and_exit();

}


void create_sparse_representation(const std::string& input_filename, double upper_value, const std::string& output_filename)
{
  dipha::inputs::full_rips_complex complex;
  complex.load_binary(input_filename, 1);

  std::vector< std::vector< std::pair<int64_t, double> > > sparse_matrix;

  int64_t matrix_size = complex.number_of_points();

  sparse_matrix.resize(matrix_size);

  for (int64_t rows = 0; rows < matrix_size; rows++)
  {

    for (int64_t cols = 0; cols < matrix_size; cols++)
    {

      double distance = complex.get_distance(rows, cols);

      if (rows != cols && distance <= 2 * upper_value)
      {
        sparse_matrix[rows].push_back(std::make_pair(cols, distance));
      }
    }
  }

  std::ofstream output_stream(output_filename.c_str(), std::ios_base::binary | std::ios_base::out);
  int64_t magic_number = 8067171840;
  output_stream.write((char*)&magic_number, sizeof(int64_t));
  int64_t file_format_id = dipha::file_types::SPARSE_DISTANCE_MATRIX;
  output_stream.write((char*)&file_format_id, sizeof(int64_t));


  output_stream.write((char*)&matrix_size, sizeof(int64_t));
  for (int64_t cur_row = 0; cur_row < matrix_size; cur_row++)
  {
    int64_t row_size = sparse_matrix[cur_row].size();
    output_stream.write((char*)&row_size, sizeof(int64_t));
  }
  for (int64_t cur_row = 0; cur_row < matrix_size; cur_row++)
  {
    int64_t row_size = sparse_matrix[cur_row].size();
    for (int64_t cur_col_idx = 0; cur_col_idx < row_size; cur_col_idx++)
    {
      output_stream.write((char*)&(sparse_matrix[cur_row][cur_col_idx].first), sizeof(int64_t));
    }
  }
  for (int64_t cur_row = 0; cur_row < matrix_size; cur_row++)
  {
    int64_t row_size = sparse_matrix[cur_row].size();
    for (int64_t cur_col_idx = 0; cur_col_idx < row_size; cur_col_idx++)
    {
      output_stream.write((char*)&(sparse_matrix[cur_row][cur_col_idx].second), sizeof(double));
    }
  }

  output_stream.close();
}

int main(int argc, char** argv)
{

  // mandatory MPI initilization call
  MPI_Init(&argc, &argv);

  std::string input_filename; // name of file that contains the weighted cell complex
  std::string output_filename; // name of file that will contain the PHAT filtration
  double upper_value = std::numeric_limits< double >::max();
  parse_command_line(argc, argv, upper_value, input_filename, output_filename);

  if (dipha::file_types::get_file_type(input_filename) != dipha::file_types::DISTANCE_MATRIX)
  {

    std::cerr << "Requires DIPHA Distance matrix file format for " << input_filename << std::endl;
    std::exit(-1);
  }

  create_sparse_representation(input_filename, upper_value, output_filename);

  MPI_Finalize();
}
