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

#pragma once

#include <dipha/includes.h>

namespace dipha
{
  namespace mpi_utils
  {
    enum message_tag
    {
      MSG_UNREDUCED_COLUMNS, MSG_CO_BOUNDARIES_ANSWERS, MSG_DUALIZED_COLS, MSG_SCATTER_QUERIES, MSG_GATHER_ANSWERS, 
      MSG_SET_GLOBAL_VALUES, MSG_TESTING, MSG_REPORT_INDICES_OF_SPARSE_RIPS, MSG_SPARSE_RIPS_DISTRIBUTE_SIMPLICES,
      MSG_REPORT_MAXIMAL_GLOBAL_INDEX_IN_SPARSE_RIPS, MSG_QUERY_SPARSE_INDICES
    };

    inline void finalize() { MPI_Finalize(); }

    inline int get_rank()
    {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      return rank;
    }

    inline bool is_root() { return get_rank() == 0; }

    struct nullstream : public std::ostream
    {
      nullstream() : std::ostream(nullptr) {};

      template <typename T>
      nullstream& operator<<(T const &) { return *this; }
    };

    inline std::ostream& cout_if_root()
    {
      static nullstream nullstreamer;
      if (is_root())
        return std::cout;
      else
        return nullstreamer;
    }

    inline std::ostream& error_printer_if_root()
    {
      static nullstream nullstreamer;
      if (is_root())
        return std::cout << std::endl << "Error: ";
      else
        return nullstreamer;
    }

    inline int get_num_processes()
    {
      int num_processes;
      MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
      return num_processes;
    }

    template< typename T >
    void receive_vector(std::vector< T >& buf,
                        int source,
                        mpi_utils::message_tag tag)
    {
      int64_t total_msg_size = 0;
      int64_t cur_msg_size = 0;
      MPI_Status status;
      int count;

      do
      {
        MPI_Probe(source, tag, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_BYTE, &count);
        cur_msg_size = count / sizeof(T);
        buf.resize(total_msg_size + cur_msg_size);
        MPI_Recv(buf.data() + total_msg_size, count, MPI_BYTE, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        total_msg_size += cur_msg_size;
      } while (cur_msg_size == globals::MPI_BLOCK_SIZE / sizeof(T));

      buf.shrink_to_fit();

      globals::bytes_received += buf.size() * sizeof(T);
    }

    template< typename T >
    void non_blocking_send_vector(const std::vector< T >& buf,
                                  int target, mpi_utils::message_tag tag,
                                  std::vector< MPI_Request >& requests)
    {
      std::vector< T >& non_const_buf = const_cast<std::vector< T >&>(buf);
      int msg_chunk_size = globals::MPI_BLOCK_SIZE / sizeof(T);
      int64_t num_chunks = buf.size() / msg_chunk_size + 1;
      for (int64_t cur_chunk = 0; cur_chunk < num_chunks; cur_chunk++)
      {
        requests.resize(requests.size() + 1);
        int bytes_to_send = (cur_chunk < num_chunks - 1) ? msg_chunk_size * sizeof(T) : (buf.size() % msg_chunk_size) * sizeof(T);
        MPI_Isend(non_const_buf.data() + cur_chunk * msg_chunk_size, bytes_to_send, 
                  MPI_BYTE, target, tag, MPI_COMM_WORLD, &requests.back());
      }
    }

    template< typename T >
    void all_gather(const T& data, std::vector< T >& gathered_data)
    {
      gathered_data.resize(mpi_utils::get_num_processes());
      MPI_Allgather((void*)&data, sizeof(T), MPI_BYTE, gathered_data.data(), sizeof(T), MPI_BYTE, MPI_COMM_WORLD);
    }

    template< typename T >
    void file_read_at_vector(const MPI_File& file,
                             const MPI_Offset& offset,
                             int64_t size,
                             std::vector< T >& content)
    {
      content.resize(size);
      int chunk_size = globals::MPI_BLOCK_SIZE / sizeof(T);
      int64_t num_chunks = size / chunk_size + 1;
      for (int64_t cur_chunk = 0; cur_chunk < num_chunks; cur_chunk++)
      {
        int bytes_to_read = (cur_chunk < num_chunks - 1) ? chunk_size * sizeof(T) : (size % chunk_size) * sizeof(T);
        MPI_Offset cur_offset = offset + cur_chunk * chunk_size * sizeof(T);
        MPI_File_read_at(file, cur_offset, content.data() + cur_chunk * chunk_size, bytes_to_read, MPI_BYTE, MPI_STATUS_IGNORE);
      }
    }

    template< typename T >
    void file_write_at_vector(const MPI_File& file,
                              const MPI_Offset& offset,
                              const std::vector< T >& content)
    {
      std::vector< T >& non_const_content = const_cast<std::vector< T >&>(content);
      int chunk_size = globals::MPI_BLOCK_SIZE / sizeof(T);
      int64_t num_chunks = content.size() / chunk_size + 1;
      for (int64_t cur_chunk = 0; cur_chunk < num_chunks; cur_chunk++)
      {
        int bytes_to_write = (cur_chunk < num_chunks - 1) ? chunk_size * sizeof(T) : (content.size() % chunk_size) * sizeof(T);
        MPI_Offset cur_offset = offset + cur_chunk * chunk_size * sizeof(T);
        MPI_File_write_at(file, cur_offset, non_const_content.data() + cur_chunk * chunk_size,
                          bytes_to_write, MPI_BYTE, MPI_STATUS_IGNORE);
      }
    }

    template< typename T >
    void file_read_at(const MPI_File& file,
                      const MPI_Offset& offset,
                      T& value)
    {
      MPI_File_read_at(file, offset, &value, sizeof(T), MPI_BYTE, MPI_STATUS_IGNORE);
    }

    template< typename T >
    void file_write_at(const MPI_File& file,
                       const MPI_Offset& offset,
                       const T& value)
    {
      MPI_File_write_at(file, offset, &const_cast<T&>(value), sizeof(T), MPI_BYTE, MPI_STATUS_IGNORE);
    }

    inline MPI_File file_open(const std::string& filename, int access)
    {
      MPI_File file;
      if (MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(filename.c_str()), access, MPI_INFO_NULL, &file) != 0)
      {
        if (mpi_utils::is_root())
        {
          std::cerr << "Error while opening file " << filename << std::endl;
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
      }
      return file;
    }

    inline MPI_File file_open_read_only(const std::string& filename) { return file_open(filename, MPI_MODE_RDONLY); }

    inline MPI_File file_open_write_only(const std::string& filename) { return file_open(filename, MPI_MODE_WRONLY | MPI_MODE_CREATE); }

    inline MPI_File file_open_read_write(const std::string& filename) { return file_open(filename, MPI_MODE_RDWR | MPI_MODE_CREATE); }
  }
}
