/* Copyright CEA-CNRS
 * Labs Involved for :
 *   - CEA  : Maison de la Simulation
 *   - CNRS : LULI, LLR, LPGP, IDRIS
 * Contributors : see file License/Contributors
 *
 * This software, SMILEI, is a PIC (Particles in Cell) program whose purpose
 * is to simulate interactions between matter and extreme intensities light.
 *
 * This software is governed by the CeCILL-B license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL-B
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-B license and that you accept its terms.
 *
 * Modified by: Andrew Sexton, 2023-10-06
 * https://github.com/SmileiPIC/Smilei
 */
#ifndef TRIFORCE_H5_WRAPPER_H
#define TRIFORCE_H5_WRAPPER_H

#include <utility>
#include <vector>
#include <string>
#include <source_location>

#include <mpi.h>
#include <hdf5.h>

#include "error_functions.h"
#include "array2d.h"
#include "vector2d.h"
#include "vector3d.h"


/*************************************/
/********** HDF5 Type Class **********/
template<typename T> struct H5Type { const static hid_t type_id; };
template<> const hid_t H5Type<float>::type_id = H5T_NATIVE_FLOAT;
template<> const hid_t H5Type<double>::type_id = H5T_NATIVE_DOUBLE;
template<> const hid_t H5Type<int>::type_id = H5T_NATIVE_INT;
template<> const hid_t H5Type<unsigned>::type_id = H5T_NATIVE_UINT;
template<> const hid_t H5Type<size_t>::type_id = H5T_NATIVE_ULONG;

/******************************************/
/********** HDF5 Dataspace Class **********/
class H5Space {
public:
  // 1D without selection or chunks
  explicit H5Space(hsize_t size)
  {
    dims_ = {size};
    global_ = size;
    sid_ = H5Screate_simple(1, &size, nullptr);
    if (size <= 0) {
      H5Sselect_none(sid_);
    }
    chunk_.resize(0);
  }

  // 1D
  H5Space(hsize_t size, hsize_t offset, hsize_t npoints, hsize_t chunk=0, bool extendable=false)
  {
    dims_ = {size};
    global_ = size;
    const hsize_t max_dims[] = {H5S_UNLIMITED};
    sid_ = H5Screate_simple(1, &size, extendable ? max_dims : nullptr);

    if (size <= 0 || npoints <= 0) {
      H5Sselect_none(sid_);
    } else {
      hsize_t count = 1;
      H5Sselect_hyperslab(sid_, H5S_SELECT_SET, &offset, nullptr, &count, &npoints);
    }

    if (chunk > 1) {
      chunk_.resize(1, chunk);
    } else {
      chunk_.resize(0);
    }
  }

  // ND
  explicit H5Space(std::vector<hsize_t> size,
                   std::vector<hsize_t> offset={},
                   std::vector<hsize_t> npoints={},
                   std::vector<hsize_t> chunk={},
                   std::vector<bool> extendable={})
  {
    dims_ = size;
    global_ = 1;
    for (auto x: size) {
      global_ *= x;
    }

    if (extendable.empty()) {
      sid_ = H5Screate_simple(static_cast<int>(size.size()), &size[0], nullptr);
    } else {
      std::vector<hsize_t> maxsize(size.size());
      for (size_t i = 0; i < size.size(); ++i) {
        maxsize[i] = extendable[i] ? H5S_UNLIMITED : size[i];
      }
      sid_ = H5Screate_simple(static_cast<int>(size.size()), &size[0], &maxsize[0]);
    }

    if (global_ <= 0) {
      H5Sselect_none(sid_);
    } else if (!offset.empty() || !npoints.empty()) {
      unsigned int i;
      for( i=0; i<npoints.size() && npoints[i]>0; i++ ) {}
      if( i < npoints.size() ) {
        H5Sselect_none( sid_ );
      } else {
        std::vector<hsize_t> count( size.size(), 1 );
        H5Sselect_hyperslab(sid_, H5S_SELECT_SET, &offset[0], nullptr, &count[0], &npoints[0]);
      }
    }
    chunk_ = std::move(chunk);
  }

  ~H5Space()
  {
    H5Sclose(sid_);
  }

  [[nodiscard]] bool valid() const
  {
    return sid_ > 0;
  }

  hid_t sid_;
  hsize_t global_;
  std::vector<hsize_t> dims_;
  std::vector<hsize_t> chunk_;
};

/*************************************/
/********** HDF5 Base Class **********/
class H5Base {
public:
  H5Base()
  : fid_(-1), id_(-1), dcr_(-1), dxpl_(-1)
  {}

  H5Base(const std::string& file, uint access, const MPI_Comm* comm, bool raise)
  {
    // Separate file name and tree inside HDF5 file
    size_t len = file.length();
    size_t idx_h5 = file.find(".h5") + 3;
    size_t idx_hdf5 = file.find(".hdf5") + 5;
    filepath = file;
    std::string group_path{};
    if (idx_h5 < std::string::npos) {
      if (idx_h5 == len || (idx_h5 < len && file[idx_h5] == '/')) {
        filepath = file.substr(0, idx_h5);
        group_path = file.substr(idx_h5);
      }
    }
    else if (idx_hdf5 < std::string::npos) {
      if (idx_hdf5 == len || (idx_hdf5 < len && file[idx_hdf5] == '/')) {
        filepath = file.substr(0, idx_hdf5);
        group_path = file.substr(idx_hdf5);
      }
    }
    // Disable HDF5 error printing stack
    H5E_auto2_t old_func;
    void *old_client_data;
    if(!raise) {
      // Backup default error printing
      H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);
      H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
    }

    // Open or create
    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    if (comm) {
      H5Pset_fapl_mpio(fapl, *comm, MPI_INFO_NULL);
    }
    if (access == H5F_ACC_RDWR) {
      fid_ = H5Fcreate(filepath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
    } else {
      fid_ = H5Fopen(filepath.c_str(), access, fapl);
    }
    H5Pclose(fapl);

    // Check for error
    if (H5Eget_num(H5E_DEFAULT) > 0) {
      fid_ = -1;
      id_ = -1;
      dcr_ = -1;
      dxpl_ = -1;
    } else {
      // Open group from file
      id_ = fid_;
      if (!group_path.empty()) {
        id_ = H5Gopen(fid_, group_path.c_str(), H5P_DEFAULT);
      }
      dxpl_ = H5Pcreate(H5P_DATASET_XFER);
      dcr_ = H5Pcreate(H5P_DATASET_CREATE);
      if (comm) {
        H5Pset_dxpl_mpio(dxpl_, H5FD_MPIO_COLLECTIVE);
        H5Pset_alloc_time(dcr_, H5D_ALLOC_TIME_EARLY);
      }
      H5Pset_fill_time(dcr_, H5D_FILL_TIME_NEVER);
      // Check error
      if (H5Eget_num(H5E_DEFAULT) > 0) {
        id_ = -1;
      }
    }

    if (raise && (fid_ < 0 || id_ < 0)) {
      ERROR("Cannot open file " + filepath, std::source_location::current());
    } else {
      // Restore previous error printing
      H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);
    }
  }

  ~H5Base()
  {
    auto id_type = H5Iget_type(id_);
    if (id_type == H5I_GROUP) {
      H5Gclose(id_);
    } else if (id_type == H5I_DATASET) {
      H5Dclose(id_);
    }

    if (fid_ >= 0) {
      H5Pclose(dxpl_);
      H5Pclose(dcr_);
      auto err = H5Fclose(fid_);
      if (err < 0) {
        H5Eprint(H5E_DEFAULT, nullptr);
        ERROR("Cannot close file " + filepath, std::source_location::current());
      }
    }
  }

  [[nodiscard]] bool valid() const
  {
    return id_ > 0;
  }

//  void flush() const
//  {
//    H5Fflush(id_, H5F_SCOPE_GLOBAL);
//  }

  [[nodiscard]] bool hasGroup(const std::string& grp_name) const {
    return H5Lexists(id_, grp_name.c_str(), H5P_DEFAULT) > 0;
  }

  [[nodiscard]] bool hasDataset(const std::string& dset_name) const {
    return hasGroup(dset_name);
  }

  [[nodiscard]] bool hasAttr(const std::string& attr_name) const {
    return H5Aexists(id_, attr_name.c_str()) > 0;
  }

protected:
  H5Base(hid_t id, hid_t dcr, hid_t dxpl)
  : fid_(-1), id_(id), dcr_(dcr), dxpl_(dxpl)
  {}

  [[nodiscard]] hid_t newGroupId(const std::string& group_name) const
  {
    if (H5Lexists(id_, group_name.c_str(), H5P_DEFAULT) > 0) {
      return H5Oopen(id_, group_name.c_str(), H5P_DEFAULT);
    } else {
      return H5Gcreate(id_, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
  }

  [[nodiscard]] hid_t open(const std::string& name) const
  {
    if (H5Lexists(id_, name.c_str(), H5P_DEFAULT) > 0) {
      return H5Oopen(id_, name.c_str(), H5P_DEFAULT);
    } else {
      ERROR("In HDF5 file, " + name + " does not exist.", std::source_location::current());
    }
    return -1;
  }

  std::string filepath; // Only defined if root '/'
  hid_t fid_; // Only defined if root '/'
  hid_t id_; // Object ID
  hid_t dcr_; // Dataset creation property list
  hid_t dxpl_; // Dataset transfer property list
};

/***************************************/
/********** HDF5 Writer Class **********/
class H5Writer : public H5Base {
public:
  // Open file + location
  explicit H5Writer(const std::string& file, MPI_Comm* comm=nullptr, bool raise=true)
  : H5Base(file, H5F_ACC_RDWR, comm, raise)
  {}

  // Create group inside given location
  H5Writer(H5Writer& loc, const std::string& group_name)
  : H5Base(loc.newGroupId(group_name), loc.dcr_, loc.dxpl_)
  {}

  // Create or open (not write) a dataset inside the given location
  H5Writer(H5Writer& loc, const std::string& name, hid_t type, H5Space& filespace)
  : H5Base(-1, loc.dcr_, loc.dxpl_)
  {
    H5D_layout_t layout = H5Pget_layout(dcr_);
    if (!filespace.chunk_.empty()) {
      H5Pset_chunk(dcr_, static_cast<int>(filespace.chunk_.size()), &(filespace.chunk_[0]));
    }
    if (H5Lexists(loc.id_, name.c_str(), H5P_DEFAULT) == 0) {
      id_ = H5Dcreate(loc.id_, name.c_str(), type, filespace.sid_, H5P_DEFAULT, dcr_, H5P_DEFAULT);
    } else {
      hid_t pid = H5Pcreate(H5P_DATASET_ACCESS);
      id_ = H5Dopen(loc.id_, name.c_str(), pid);
      H5Pclose(pid);
    }
    H5Pset_layout(dcr_, layout);
  }

  // Location is already open
  H5Writer(hid_t id, hid_t dcr, hid_t dxpl)
  : H5Base(id, dcr, dxpl)
  {}

  ~H5Writer() = default;

  H5Writer newGroup(const std::string& group_name)
  {
    return {*this, group_name};
  }

  H5Writer newDataset(const std::string& dset_name, hid_t type, H5Space& filespace)
  {
    return {*this, dset_name, type, filespace};
  }

  // Write string to attribute
  void addAttr(const std::string& attr_name, const std::string& attr_val)
  {
    hid_t atype = H5Tcopy(H5T_C_S1);
    if (attr_val.empty()) {
      H5Tset_size(atype, 1);
    } else {
      H5Tset_size(atype, attr_val.size());
    }
    hid_t sid = H5Screate(H5S_SCALAR);
    hid_t aid = H5Acreate(id_, attr_name.c_str(), atype, sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(aid, atype, attr_val.c_str());
    H5Sclose(sid);
    H5Aclose(aid);
    H5Tclose(atype);
  }

  // Write char* to attribute
  void addAttr(const std::string& attr_name, const char* attr_val)
  {
    addAttr(attr_name, std::string{attr_val});
  }

  // Write any type to attribute
  template<typename T>
  void addAttr(const std::string& attr_name, const T& attr_val)
  {
    auto type = H5Type<T>::type_id;
    hid_t sid = H5Screate(H5S_SCALAR);
    hid_t aid = H5Acreate(id_, attr_name.c_str(), type, sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(aid, type, &attr_val);
    H5Sclose(sid);
    H5Aclose(aid);
  }

  // Write vector of strings to attribute
  void addAttr(const std::string& attr_name, const std::vector<std::string>& attr_val)
  {
    std::string full_string{};
    for (auto& str: attr_val) {
      full_string += str + ";";
    }
    addAttr(attr_name, full_string);
  }

  // Write vector of any type to attribute
  template<typename T>
  void addAttr(const std::string& attr_name, const std::vector<T>& attr_value)
  {
    auto type = H5Type<T>::type_id;
    hsize_t dims = attr_value.size();
    hid_t sid = H5Screate_simple(1, &dims, nullptr);
    hid_t aid = H5Acreate(id_, attr_name.c_str(), type, sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(aid, type, &attr_value[0]);
    H5Aclose(aid);
    H5Sclose(sid);
  }

  template<typename T>
  H5Writer addVector(const std::string& name, const std::vector<vec2<T>>& v, H5Space& filespace, H5Space& memspace)
  {
//    auto type = H5Type<T>::type_id;
//    auto vec_t = H5Tcreate(H5T_COMPOUND, sizeof(vec2<T>));
//    H5Tinsert(vec_t, "x", HOFFSET(vec2<T>, e[0]), type);
//    H5Tinsert(vec_t, "z", HOFFSET(vec2<T>, e[1]), type);
//
//    auto dset_id = H5Dcreate(id_, name.c_str(), vec_t, filespace.sid_, H5P_DEFAULT, dcr_, H5P_DEFAULT);
//    H5Dwrite(dset_id, vec_t, memspace.sid_, filespace.sid_, dxpl_, &v[0]);
//
//    H5Tclose(vec_t);
//    return H5Writer(dset_id, dcr_, dxpl_);
    // Create array type
    size_t dims[1] = {2};
    auto arr_type = H5Tarray_create(H5Type<T>::type_id, 1, dims);
    auto dset_id = H5Dcreate(id_, name.c_str(), arr_type, filespace.sid_, H5P_DEFAULT, dcr_, H5P_DEFAULT);
    // Write dataset
    H5Dwrite(dset_id, arr_type, memspace.sid_, filespace.sid_, dxpl_, &v[0]);
    // Close array type
    H5Tclose(arr_type);

    return H5Writer(dset_id, dcr_, dxpl_);
  }

  template<typename T>
  H5Writer addVector(const std::string& name, const std::vector<vec3<T>>& v, H5Space& filespace, H5Space& memspace)
  {
    // Create array type
    size_t dims[1] = {3};
    auto arr_type = H5Tarray_create(H5Type<T>::type_id, 1, dims);
    auto dset_id = H5Dcreate(id_, name.c_str(), arr_type, filespace.sid_, H5P_DEFAULT, dcr_, H5P_DEFAULT);
    // Write dataset
    H5Dwrite(dset_id, arr_type, memspace.sid_, filespace.sid_, dxpl_, &v[0]);
    // Close array type
    H5Tclose(arr_type);

    return H5Writer(dset_id, dcr_, dxpl_);
  }

  // write a vector of any type to dataset
  template<typename T>
  H5Writer addVector(const std::string& name, const std::vector<T>& v, H5Space& filespace, H5Space& memspace)
  {
    return addVector<T>(name, v[0], filespace, memspace);
  }

  // write a portion of a vector of any type to dataset
  template<typename T>
  H5Writer addVector(const std::string& name, const T& v, H5Space& filespace, H5Space& memspace)
  {
//    // Compression Settings
//    if (deflate > 0) {
//      H5Pset_chunk(dcr_, 1, &dim);
//      H5Pset_deflate(dcr_, std::min(9, deflate));
//    } else {
//      H5Premove_filter(dcr_, H5Z_FILTER_DEFLATE);
//    }

    // Create dataset
    auto type = H5Type<T>::type_id;
    auto dset_id = H5Dcreate(id_, name.c_str(), type, filespace.sid_, H5P_DEFAULT, dcr_, H5P_DEFAULT);
    // Write dataset
    H5Dwrite(dset_id, type, memspace.sid_, filespace.sid_, dxpl_, &v);

    return H5Writer(dset_id, dcr_, dxpl_);
  }

  template<typename T>
  H5Writer addArray(const std::string& name, const std::vector<T>& v, H5Space& filespace, H5Space& memspace, bool independent=false)
  {
    return addArray<T>(name, v[0], filespace, memspace, independent);
  }

  template<typename T>
  H5Writer addArray(const std::string& name, const Array2D<T>& v, H5Space& filespace, H5Space& memspace, bool independent=false)
  {
    return addArray<T>(name, v[0], filespace, memspace, independent);
  }

  template<typename T>
  H5Writer addArray(const std::string& name, const T& v, H5Space& filespace, H5Space& memspace, bool independent=false)
  {
    auto type = H5Type<T>::type_id;
    H5Writer dset = newDataset(name, type, filespace);
    dset.write(v, filespace, memspace, independent);
    return dset;
  }

  // Write to an open dataset
  template<typename T>
  void write(const T& v, const H5Space& filespace, const H5Space& memspace, bool independent=false)
  {
    auto type = H5Type<T>::type_id;
    if (independent) {
      if (memspace.global_ > 0) {
        H5FD_mpio_xfer_t xfer;
        H5Pget_dxpl_mpio(dxpl_, &xfer);
        H5Pset_dxpl_mpio(dxpl_, H5FD_MPIO_INDEPENDENT);
        H5Dwrite(id_, type, memspace.sid_, filespace.sid_, dxpl_, &v);
        H5Pset_dxpl_mpio(dxpl_, xfer);
      }
    } else {
      if (filespace.global_ > 0) {
        H5Dwrite(id_, type, memspace.sid_, filespace.sid_, dxpl_, &v);
      }
    } // end if (independent)
  } // end write()

  // read from an open dataset
  template<typename T>
  void read(T& v, const H5Space& filespace, const H5Space& memspace, bool independent=false)
  {
    auto type = H5Type<T>::type_id;
    if (independent) {
      if (memspace.global_ > 0) {
        H5FD_mpio_xfer_t xfer;
        H5Pget_dxpl_mpio(dxpl_, &xfer);
        H5Pset_dxpl_mpio(dxpl_, H5FD_MPIO_INDEPENDENT);
        H5Dread(id_, type, memspace.sid_, filespace.sid_, dxpl_, &v);
        H5Pset_dxpl_mpio(dxpl_, xfer);
      }
    } else {
      if (filespace.global_ > 0) {
        H5Dread(id_, type, memspace.sid_, filespace.sid_, dxpl_, &v);
      }
    } // end if (independent)
  } // end read()
}; // end class H5Writer

/***************************************/
/********** HDF5 Reader Class **********/
class H5Reader : public H5Base {
public:
  H5Reader() : H5Base() {}

  explicit H5Reader(const std::string& file, const MPI_Comm* comm=nullptr, bool raise=true)
  : H5Base(file, H5F_ACC_RDONLY, comm, raise)
  {}

  H5Reader(hid_t id, hid_t dcr, hid_t dxpl)
  : H5Base(id, dcr, dxpl)
  {}

  H5Reader openGroup(const std::string& group_name)
  {
    return {open(group_name), dcr_, dxpl_};
  }

  H5Reader openDataset(const std::string& dset_name)
  {
    return {open(dset_name), dcr_, dxpl_};
  }

  void openAttr(const std::string& attr_name, std::string& attr_value)
  {
    if (H5Aexists(id_, attr_name.c_str()) > 0) {
      hid_t aid = H5Aopen_name(id_, attr_name.c_str());
      hid_t attr_type = H5Aget_type(aid);
      hsize_t sdim = H5Tget_size(attr_type);
      hid_t mem_type = H5Tcopy(H5T_C_S1);
      H5Tset_size(mem_type, sdim);
      std::vector<char> tmp(sdim);
      H5Aread(aid, mem_type, &tmp[0]);
      attr_value = std::string(tmp.begin(), tmp.end());
      H5Tclose(mem_type);
      H5Tclose(attr_type);
      H5Aclose(aid);
    } else {
      ERROR("Cannot find attribute " + attr_name, std::source_location::current());
    }
  }

  void openAttr(const std::string& attr_name, std::vector<std::string>& attr_val)
  {
    std::string full_string{};
    openAttr(attr_name, full_string);
    size_t back = 0;
    for (size_t i = 0; i < full_string.length(); ++i) {
      if (full_string[i] == ';') {
        attr_val.push_back(full_string.substr(back, i - back));
        back = i + 1;
      }
    }
  }

  template<typename T>
  void openAttr(const std::string& attr_name, std::vector<T>& attr_val)
  {
    if (H5Aexists(id_, attr_name.c_str()) > 0) {
      auto type = H5Type<T>::type_id;
      hid_t aid = H5Aopen(id_, attr_name.c_str(), H5P_DEFAULT);
      hid_t sid = H5Aget_space(aid);
      hssize_t npoints = H5Sget_simple_extent_npoints(sid);
      attr_val.resize(npoints);
      H5Aread(aid, type, &(attr_val[0]));
      H5Sclose(sid);
      H5Aclose(aid);
    } else {
      ERROR("Cannot find attribute " + attr_name, std::source_location::current());
    }
  }

  template<typename T>
  void openAttr(const std::string& attr_name, T& attr_val)
  {
    if (H5Aexists(id_, attr_name.c_str()) > 0) {
      auto type = H5Type<T>::type_id;
      hid_t aid = H5Aopen(id_, attr_name.c_str(), H5P_DEFAULT);
      H5Aread(aid, type, &attr_val);
      H5Aclose(aid);
    } else {
      ERROR("Cannot find attribute " + attr_name, std::source_location::current());
    }
  }

  template<typename T>
  void loadVector(const std::string& vec_name, std::vector<T>& v, hsize_t offset=0, hsize_t npoints=0)
  {
    auto n = getSize(vec_name);
    if (npoints == 0) {
      npoints = n - offset;
    }
    v.resize(npoints);
    loadVector(vec_name, v[0], offset, npoints);
  }

  template<typename T>
  void loadVector(const std::string& vec_name, T& v, hsize_t offset= 0, hsize_t npoints= 0)
  {
    hid_t did = H5Dopen(id_, vec_name.c_str(), H5P_DEFAULT);
    if (did < 0) {
      ERROR("Cannot read dataset " + vec_name, std::source_location::current());
    }
    auto type = H5Type<T>::type_id;
    if (offset != 0 || npoints != 0) {
      hid_t sid = H5Dget_space(did);
      int sdim = H5Sget_simple_extent_ndims(sid);
      if (sdim != 1) {
        ERROR("Expected 1D dataset shape, but got " + std::to_string(sdim) + "D.", std::source_location::current());
      }
      hsize_t dim = H5Sget_simple_extent_npoints(sid);
      if (npoints == 0) {
        npoints = dim - offset;
      }
      hid_t memspace = H5Screate_simple(1, &npoints, nullptr);
      hid_t filespace = H5Screate_simple(1, &dim, nullptr);
      hsize_t o = offset;
      hsize_t c = 1;
      hsize_t n = npoints;
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &o, nullptr, &c, &n);
      H5Dread(did, type, H5S_ALL, H5S_ALL, dxpl_, &v);
      H5Sclose(sid);
      H5Sclose(filespace);
      H5Sclose(memspace);
    } else {
      H5Dread(did, type, H5S_ALL, H5S_ALL, dxpl_, &v);
    }
    H5Dclose(did);
  }

  template<typename T>
  void loadArray(const std::string& arr_name, std::vector<T>& v, H5Space& filespace, H5Space& memspace)
  {
    auto npoints = getSize(arr_name);
    v.resize(npoints);
    loadArray(arr_name, v[0], filespace, memspace);
  }

  template<typename T>
  void loadArray(const std::string& arr_name, Array2D<T>& v, H5Space& filespace, H5Space& memspace)
  {
    // auto npoints = getSize(arr_name);
    // assert(v.size() == npoints);
    loadArray(arr_name, v[0], filespace, memspace);
  }

  template<typename T>
  void loadArray(const std::string& arr_name, T& v, H5Space& filespace, H5Space& memspace)
  {
    auto type = H5Type<T>::type_id;
    hid_t did = open(arr_name);
    if (filespace.global_ > 0) {
      H5Dread(did, type, memspace.sid_, filespace.sid_, dxpl_, &v);
    }
    H5Dclose(did);
  }

  std::vector<hsize_t> getShape(const std::string& name)
  {
    std::vector<hsize_t> shape{};
    if (H5Lexists(id_, name.c_str(), H5P_DEFAULT) > 0) {
      hid_t did = H5Dopen(id_, name.c_str(), H5P_DEFAULT);
      if (did >= 0) {
        hid_t sid = H5Dget_space(did);
        int sdim = H5Sget_simple_extent_ndims(sid);
        shape.resize(sdim);
        H5Sget_simple_extent_dims(sid, &shape[0], nullptr);
        H5Sclose(sid);
        H5Dclose(did);
      }
    }
    return shape;
  }

  hsize_t getSize(const std::string& name)
  {
    auto shape = getShape(name);
    hsize_t n = 1;
    for (auto x: shape) {
      n *= x;
    }
    return n;
  }
};

#endif //TRIFORCE_H5_WRAPPER_H
